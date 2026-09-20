"""Explicit optional wall RHS backends; Python remains the physical reference."""
from dataclasses import dataclass, asdict
from functools import lru_cache
import hashlib
import importlib.util
import json
import os
import tempfile
import inspect
import math
import platform
from pathlib import Path
import sys
import time
import numpy as np
from . import numerical_primitives as numeric


@dataclass(frozen=True)
class WallBackendSettings:
    name: str = 'python'
    profile: bool = False
    disk_cache: bool = False
    cache_directory: str | None = None

    def __post_init__(self):
        if self.name not in ('python','numba'):
            raise ValueError('Wall backend must be python or numba.')
        if type(self.disk_cache) is not bool:
            raise ValueError('Backend disk_cache must be a boolean.')
        if self.disk_cache and self.name!='numba':
            raise ValueError('Disk JIT cache requires the numba backend.')
        if self.cache_directory is not None and (not isinstance(self.cache_directory,str) or not self.cache_directory):
            raise ValueError('Cache directory must be a nonempty path string.')
        if type(self.profile) is not bool:
            raise ValueError('Backend profile must be a boolean.')


def backend_identity(settings=WallBackendSettings()):
    identity=dict(settings=asdict(settings),python=platform.python_version(),
        numpy=np.__version__,platform=platform.platform(),libc=list(platform.libc_ver()),source_sha256=hashlib.sha256(
            Path(__file__).read_bytes()+Path(numeric.__file__).read_bytes()).hexdigest(),
        fastmath=False,parallel=False,float_type='float64')
    if settings.name=='numba':
        try:
            import numba
            import llvmlite
            from llvmlite import binding as llvm
        except ImportError:
            identity.update(numba_available=False)
        else:
            identity.update(numba_available=True,numba=numba.__version__,llvmlite=llvmlite.__version__,
                llvm=list(llvm.llvm_version_info),cpu=str(llvm.get_host_cpu_name()),
                cpu_features=llvm.get_host_cpu_features().flatten(),
                disable_jit=bool(numba.config.DISABLE_JIT),
                compiler_cpu_name=numba.config.CPU_NAME,compiler_cpu_features=numba.config.CPU_FEATURES,
                compiler_opt=str(numba.config.OPT),boundscheck=numba.config.BOUNDSCHECK,
                numba_cache_dir=numba.config.CACHE_DIR)
    return identity


def cache_namespace(identity, source=None):
    """Content + compiler/runtime/CPU addressed; never trust only file mtime."""
    source=Path(numeric.__file__).read_bytes() if source is None else source
    runtime={k:v for k,v in identity.items() if k!='settings'}
    return hashlib.sha256(source+json.dumps(runtime,sort_keys=True).encode()).hexdigest()


@lru_cache(maxsize=8)
def compiled_dispatcher(cache_directory=None, namespace=None):
    import numba
    from numba.extending import register_jitable
    module=numeric
    if cache_directory is not None:
        # A content-addressed copy of the SAME authoritative module avoids
        # Numba's cross-file and timestamp-only invalidation limitations.
        directory=Path(cache_directory)/namespace
        directory.mkdir(parents=True,exist_ok=True)
        name='_dada_numeric_'+namespace
        path=directory/(name+'.py')
        source=Path(numeric.__file__).read_bytes()
        if not path.exists() or path.read_bytes()!=source:
            with tempfile.NamedTemporaryFile(dir=directory,delete=False) as handle:
                handle.write(source);temporary=Path(handle.name)
            temporary.replace(path)
        if name in sys.modules:
            module=sys.modules[name]
        else:
            spec=importlib.util.spec_from_file_location(name,path)
            module=importlib.util.module_from_spec(spec)
            sys.modules[name]=module
            try: spec.loader.exec_module(module)
            except BaseException:
                sys.modules.pop(name,None)
                raise
    # Register compilation support without replacing ordinary Python bindings.
    for name,function in vars(module).items():
        if inspect.isfunction(function) and function.__module__==module.__name__ and name!='wall_kernel':
            register_jitable(function)
    return numba.njit(fastmath=False,parallel=False,cache=cache_directory is not None)(module.wall_kernel)


class PythonWallRHS:
    def __init__(self,wrapper,reason=None,profile=False):
        self.wrapper=wrapper;self.reason=reason;self.profile=profile
        self.calls=0;self.fallbacks=0;self.total_seconds=0.
        self.first_call_seconds=None;self.kinematics_seconds=None;self.dispatch_seconds=None

    def __call__(self,angle,values):
        self.calls+=1
        if self.reason is not None: self.fallbacks+=1
        start=time.perf_counter() if self.profile else 0.
        try: return self.wrapper.derivative(angle,values)
        finally:
            if self.profile: self.total_seconds+=time.perf_counter()-start

    def snapshot(self):
        return dict(actual_backend='python',reason=self.reason,calls=self.calls,
            fallback_calls=self.fallbacks,profile=self.profile,
            total_rhs_seconds=self.total_seconds if self.profile else None)


class WallRHS:
    """One candidate-owned callable and its observable counters, never global."""
    def __init__(self,wrapper,settings=WallBackendSettings()):
        if not isinstance(settings,WallBackendSettings):
            raise TypeError('Expected WallBackendSettings.')
        self.settings=settings
        start=time.perf_counter()
        self.identity=backend_identity(settings)
        self.cache_path=None;self.cache_error=None;self.cache_key=None
        if settings.name=='python':
            self.implementation=PythonWallRHS(wrapper,profile=settings.profile)
        elif not self.identity.get('numba_available'):
            self.implementation=PythonWallRHS(wrapper,'numba_unavailable',settings.profile)
        elif self.identity.get('disable_jit'):
            self.implementation=PythonWallRHS(wrapper,'numba_disabled',settings.profile)
        else:
            if settings.disk_cache:
                self.cache_path=str(Path(settings.cache_directory or
                    str(Path(os.environ.get('XDG_CACHE_HOME',Path.home()/'.cache'))/'dada_solver'/'numba')).expanduser().resolve())
                self.cache_key=cache_namespace(self.identity)
            try:
                try: dispatcher=compiled_dispatcher(self.cache_path,self.cache_key)
                except OSError as error:
                    # A read-only/full cache must not disable the numerical path.
                    self.cache_error=str(error);self.cache_path=None
                    dispatcher=compiled_dispatcher()
                self.implementation=CompiledWallRHS(wrapper,dispatcher,settings.profile)
            except TypeError as error:
                self.implementation=PythonWallRHS(wrapper,'unsupported_model: '+str(error),settings.profile)
        self.preparation_seconds=time.perf_counter()-start

    def __call__(self,angle,values):
        return self.implementation(angle,values)

    def snapshot(self):
        obj=self.implementation
        if isinstance(obj,PythonWallRHS): result=obj.snapshot()
        else:
            result=dict(actual_backend='numba',calls=obj.calls,fallback_calls=obj.fallbacks,
                fallback_reason='unsupported_state' if obj.fallbacks else None,
                first_call_seconds=obj.first_call_seconds,profile=obj.profile,
                dispatcher_cache_path=str(obj.dispatcher.stats.cache_path),
                dispatcher_cache_hits=sum(obj.dispatcher.stats.cache_hits.values()),
                dispatcher_cache_misses=sum(obj.dispatcher.stats.cache_misses.values()),
                kinematics_seconds=obj.kinematics_seconds if obj.profile else None,
                dispatch_seconds=obj.dispatch_seconds if obj.profile else None,
                total_rhs_seconds=obj.total_seconds if obj.profile else None)
        return dict(result,requested_backend=self.settings.name,identity=self.identity,
                    preparation_seconds=self.preparation_seconds,
                    cache_directory=self.cache_path,cache_namespace=self.cache_key,cache_error=self.cache_error)


class CompiledWallRHS:
    """Candidate-owned arrays; unsupported model families are explicitly refused."""
    def __init__(self,wrapper,dispatcher,profile=False):
        from dada_solver.dynamics import ThermodynamicModel
        from dada_solver.fluids import CaloricallyPerfectGas
        from dada_solver.exchangers.air_wall import AirWallMotor,AirWallExchanger
        from dada_solver.exchangers.hardware import TubeHalfLink
        from dada_solver.exchangers.gas_transport import DiluteGasTransport
        from dada_solver.exchangers.gas_correlations import MicrotubeGasModel
        from dada_solver.exchangers.gas_film import MicrotubeGasFilm
        from dada_solver.exchangers.microtube_geometry import MicrotubeBank
        from dada_solver.valves import PassiveCheckValve
        from dada_solver.hydraulics import CompressibleOrifice
        model=wrapper.model
        if (type(wrapper) is not AirWallMotor or type(model) is not ThermodynamicModel or
                type(model.gas) is not CaloricallyPerfectGas or not model.continuous_ideal_diodes):
            raise TypeError('Compiled backend requires the built-in ideal-diode air-wall model.')
        def checked(gas_model):
            if (type(gas_model) is not MicrotubeGasModel or type(gas_model.transport) is not DiluteGasTransport
                    or gas_model.slip is not None):
                raise TypeError('Compiled backend requires built-in dilute transport without an explicit slip model.')
            return gas_model.transport
        for valve in (model.hot_small_valve,model.cold_large_valve):
            if type(valve) is not PassiveCheckValve: raise TypeError('Unsupported valve class.')
        links=[]
        for link in (model.large_hot_link,model.small_cold_link,model.hot_small_valve.flow_model,model.cold_large_valve.flow_model):
            if type(link) is not TubeHalfLink or type(link.bank) is not MicrotubeBank:
                raise TypeError('Compiled backend requires built-in tube links.')
            if type(link._flow_cap) is not CompressibleOrifice or link._flow_cap.pressure_regularization!=0:
                raise TypeError('Compiled backend requires an unregularized sonic cap.')
            tr=checked(link.gas_model);g=link.gas_model;b=link.bank
            if not math.isclose(tr.gas_constant,model.gas.gas_constant,rel_tol=.005):
                raise TypeError('Transport species does not match thermodynamic gas.')
            links.append((b.inner_diameter_m,b.tube_length_m,b.tube_count,b.tube_flow_area_m2,
                link.core_loss_multiplier,link.header_loss_coefficient,link.valve_cda_m2 or 0.,
                g.maximum_mach,g.maximum_relative_pressure_drop,float(g.thermal_entry),tr.minimum_temperature,tr.maximum_temperature,numeric.SPECIES.index(tr.species)))
        walls=[]
        for wall in (wrapper.heat_in,wrapper.heat_out):
            if type(wall) is not AirWallExchanger or type(wall.gas_film) is not MicrotubeGasFilm:
                raise TypeError('Compiled backend requires built-in variable gas films.')
            film=wall.gas_film;g=film.model;tr=checked(g);b=film.bank
            if type(b) is not MicrotubeBank: raise TypeError('Unsupported bank class.')
            walls.append((wall.wall_capacity_j_k,wall.air_inlet_temperature_k,wall._effective_air_conductance,
                film.half_wall_resistance_k_w,b.tube_internal_area_m2,b.inner_diameter_m,b.tube_length_m,
                b.tube_flow_area_m2,g.maximum_mach,g.maximum_relative_pressure_drop,float(g.thermal_entry),
                tr.minimum_temperature,tr.maximum_temperature,numeric.SPECIES.index(tr.species)))
        self.wrapper=wrapper;self.reference=wrapper.derivative
        self.dispatcher=dispatcher;self.profile=profile
        self.kinematics_seconds=0.;self.dispatch_seconds=0.;self.total_seconds=0.
        self.sources=np.array((1,0,3,2),dtype=np.int64)
        self.destinations=np.array((3,2,0,1),dtype=np.int64)
        self.one_way=np.array((False,False,True,True))
        self.ports=np.array(((1,3),(0,2)),dtype=np.int64)
        self.links=np.array(links);self.walls=np.array(walls)
        g=model.gas
        self.gas=np.array((g.gas_constant,g.heat_capacity_cv,g.heat_capacity_cp,g.heat_capacity_ratio,model.angular_speed))
        self.calls=0;self.fallbacks=0
        self.first_call_seconds=None
        for array in (self.links,self.walls,self.gas,self.sources,self.destinations,self.one_way,self.ports): array.flags.writeable=False

    def __call__(self,angle,values):
        full_start=time.perf_counter() if self.profile else 0.
        self.calls+=1
        values=np.asarray(values,dtype=float)
        if values.ndim!=1 or values.size<10:
            self.fallbacks+=1
            return self.reference(angle,values)
        started=time.perf_counter() if self.first_call_seconds is None else None
        model=self.wrapper.model;k=model.kinematics
        kin_start=time.perf_counter() if self.profile else 0.
        provider=getattr(k,'cylinder_volumes_and_derivatives',None)
        if provider is not None:
            vs,vl,ds,dl=provider(angle)
        else:
            vs,vl=k.small_cylinder_volume(angle),k.large_cylinder_volume(angle)
            ds,dl=k.small_cylinder_volume_derivative(angle),k.large_cylinder_volume_derivative(angle)
        if self.profile: self.kinematics_seconds+=time.perf_counter()-kin_start
        volumes=np.array((vs,vl,model.machine_volumes.cold_heat_exchanger,model.machine_volumes.hot_heat_exchanger))
        rates=np.array((model.angular_speed*ds,model.angular_speed*dl))
        dispatch_start=time.perf_counter() if self.profile else 0.
        ok,result=self.dispatcher(np.asarray(values,dtype=float),volumes,rates,self.gas,self.links,self.walls,self.sources,self.destinations,self.one_way,self.ports)
        if self.profile: self.dispatch_seconds+=time.perf_counter()-dispatch_start
        if started is not None: self.first_call_seconds=time.perf_counter()-started
        if not ok:
            self.fallbacks+=1
            result=self.reference(angle,values)
        if self.profile: self.total_seconds+=time.perf_counter()-full_start
        return result
