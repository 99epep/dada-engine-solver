"""Explicit research prototype: compiled air/continuum/laminar wall RHS.

Not a production backend. This deliberately limited numerical port is compared
against the authoritative Python implementation before any adoption. Unsupported
objects/states run that implementation, including all its errors and diagnostics.
Kinematics, breakpoints, integration, periodic convergence and reports stay Python.
No fastmath, parallel execution, disk JIT cache, or approximate state cache.
"""
from contextlib import contextmanager
import math
import time
from unittest.mock import patch
import numpy as np
from numba import njit



@njit
def molar_cp(t, oxygen):
    if oxygen:
        a,b,c,d,e = ((31.32234,-20.23531,57.86644,-36.50624,-.007374) if t<700 else
                     (30.03235,8.772972,-3.988133,.788313,-.741599))
    else:
        a,b,c,d,e = ((28.98641,1.853978,-9.647459,16.63537,.000117) if t<500 else
                     (19.50583,19.88705,-8.598535,1.369784,.527601))
    x=t/1000
    return a+b*x+c*x*x+d*x**3+e/x**2


@njit
def properties(t):
    mu=1.716e-5*(t/273.)**1.5*(273.+111.)/(t+111.)
    k=.0241*(t/273.)**1.5*(273.+194.)/(t+194.)
    cp=(.79*molar_cp(t,False)+.21*molar_cp(t,True))/(.79*.0280134+.21*.0319988)
    return mu,k,cp


@njit
def laminar_diagnostics(flow,p1,p2,t,d,length,area,maximum_mach,maximum_drop,thermal_entry):
    mu,k,cp=properties(t)
    rho=min(p1,p2)/(287.05*t)
    re=abs(flow)*d/(area*mu);pr=cp*mu/k
    ma=abs(flow)/(rho*area)/math.sqrt(cp/(cp-287.05)*287.05*t)
    kn=mu/min(p1,p2)*math.sqrt(math.pi*287.05*t/2)/d
    ratio=max(p1,p2)/min(p1,p2)
    valid=(re<2300 and ma<=maximum_mach and kn<.001 and
           2*(ratio-1)/(ratio+1)<=maximum_drop and .5<=pr<=2000 and length>=.05*re*d)
    gz=re*pr*d/length
    nu=3.66+.0668*gz/(1+.04*gz**(2/3)) if thermal_entry else 3.66
    return valid,nu


@njit
def directed_flow(pin,pout,t,p,r,gamma):
    # p: diameter, length, count, area, multiplier, header K, valve CdA,
    # maximum Mach, maximum pressure drop, thermal-entry flag, Tmin, Tmax.
    if pout>=pin:
        return True,0.
    d,length,count,area,multiplier,header,cda=p[:7]
    if not p[10]<=t<=p[11]:
        return False,0.
    mu,_,_=properties(t)
    kn=mu/((pin+pout)/2)*math.sqrt(math.pi*287.05*t/2)/d
    if kn>=.001:
        return False,0.
    rho=(pin+pout)/(2*r*t)
    nominal=count*math.pi*d**4*(pin-pout)*(pin+pout)/(256*mu*(length/2)*r*t)/multiplier
    linear=(pin-pout)/nominal
    quadratic=header/(4*rho*area**2)
    if cda>0: quadratic+=1/(2*rho*cda**2)
    dp=pin-pout
    flow=2*dp/(linear+math.sqrt(linear**2+4*quadratic*dp))
    # Conservative prototype domain: fall back even for unsupported report-only
    # states. Python owns transition/turbulent/slip treatment and failure strings.
    valid,_=laminar_diagnostics(flow,pin,pout,t,d,length,area,p[7],p[8],p[9])
    if not valid:
        return False,0.
    ratio=pout/pin
    cap_area=min(area,cda) if cda>0 else area
    if ratio<=(2./(gamma+1.))**(gamma/(gamma-1.)):
        factor=math.sqrt(gamma/(r*t))
        factor*=(2./(gamma+1.))**((gamma+1.)/(2.*(gamma-1.)))
        cap=cap_area*pin*factor
    else:
        radicand=2.*gamma/(r*t*(gamma-1.))*(ratio**(2./gamma)-ratio**((gamma+1.)/gamma))
        cap=cap_area*pin*math.sqrt(max(0.,radicand))
    return True,min(flow,cap)


@njit
def kernel(values,volumes,volume_rates,gas,links,walls):
    result=np.zeros(15)
    t=np.empty(4);pressure=np.empty(4)
    r,cv,cp,gamma,omega=gas
    for j in range(4):
        m,u=values[2*j],values[2*j+1]
        if not math.isfinite(m) or not math.isfinite(u) or m<=0 or u<=0 or not math.isfinite(volumes[j]) or volumes[j]<=0:
            return False,result
        t[j]=u/(m*cv)
        pressure[j]=m*r*t[j]/volumes[j]
    flows=np.empty(4)
    # Same transport accumulation order as ThermodynamicModel.assemble_rates.
    sources=(1,0,3,2);destinations=(3,2,0,1)
    for link in range(4):
        a,b=sources[link],destinations[link]
        reverse=pressure[a]<pressure[b]
        if reverse and link>=2:
            flows[link]=0.
            continue
        if reverse: a,b=b,a
        ok,flow=directed_flow(pressure[a],pressure[b],t[a],links[link],r,gamma)
        if not ok: return False,result
        flows[link]=-flow if reverse else flow
        enthalpy=cp*t[a]
        result[2*a]-=flow;result[2*b]+=flow
        result[2*a+1]-=flow*enthalpy;result[2*b+1]+=flow*enthalpy
    for side in range(2):
        j=side+2
        # wall: capacity, air inlet, air conductance, film resistance, internal
        # area, diameter, length, flow area, Mach/drop guards, entry, Tmin/Tmax.
        w=walls[side]
        tw=values[8+side]/w[0]
        if not math.isfinite(tw) or tw<=0 or not w[11]<=t[j]<=w[12]:
            return False,result
        _,k,_=properties(t[j]);conductance=0.
        for port in range(2):
            link=(1 if side==0 else 0) if port==0 else (3 if side==0 else 2)
            a,b=sources[link],destinations[link]
            p1,p2=pressure[a],pressure[b]
            if flows[link]==0: p1,p2=pressure[j],pressure[j]
            ok,nu=laminar_diagnostics(flows[link],p1,p2,t[j],w[5],w[6],w[7],w[8],w[9],w[10])
            if not ok: return False,result
            conductance+=.5*w[4]*nu*k/w[5]
        overall=1/(1/conductance+w[3])
        air_heat=w[2]*(w[1]-tw);gas_heat=overall*(tw-t[j])
        result[2*j+1]+=gas_heat
        result[8+side]=air_heat-gas_heat
        result[10+side]=air_heat;result[12+side]=gas_heat
    result[1]-=pressure[0]*volume_rates[0]
    result[3]-=pressure[1]*volume_rates[1]
    result[14]=pressure[0]*volume_rates[0]+pressure[1]*volume_rates[1]
    return True,result/omega


class PreparedRHS:
    """Candidate-owned arrays; unsupported model families are explicitly refused."""
    def __init__(self,wrapper,reference):
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
            raise TypeError('Prototype requires the built-in ideal-diode air-wall model.')
        def checked(gas_model):
            if (type(gas_model) is not MicrotubeGasModel or type(gas_model.transport) is not DiluteGasTransport
                    or gas_model.transport.species!='air' or gas_model.slip is not None):
                raise TypeError('Prototype supports built-in dilute air without slip only.')
            return gas_model.transport
        for valve in (model.hot_small_valve,model.cold_large_valve):
            if type(valve) is not PassiveCheckValve: raise TypeError('Unsupported valve class.')
        links=[]
        for link in (model.large_hot_link,model.small_cold_link,model.hot_small_valve.flow_model,model.cold_large_valve.flow_model):
            if type(link) is not TubeHalfLink or type(link.bank) is not MicrotubeBank:
                raise TypeError('Prototype requires built-in tube links.')
            if type(link._flow_cap) is not CompressibleOrifice or link._flow_cap.pressure_regularization!=0:
                raise TypeError('Prototype requires an unregularized sonic cap.')
            tr=checked(link.gas_model);g=link.gas_model;b=link.bank
            if not math.isclose(tr.gas_constant,model.gas.gas_constant,rel_tol=.005):
                raise TypeError('Transport species does not match thermodynamic gas.')
            links.append((b.inner_diameter_m,b.tube_length_m,b.tube_count,b.tube_flow_area_m2,
                link.core_loss_multiplier,link.header_loss_coefficient,link.valve_cda_m2 or 0.,
                g.maximum_mach,g.maximum_relative_pressure_drop,float(g.thermal_entry),tr.minimum_temperature,tr.maximum_temperature))
        walls=[]
        for wall in (wrapper.heat_in,wrapper.heat_out):
            if type(wall) is not AirWallExchanger or type(wall.gas_film) is not MicrotubeGasFilm:
                raise TypeError('Prototype requires built-in variable gas films.')
            film=wall.gas_film;g=film.model;tr=checked(g);b=film.bank
            if type(b) is not MicrotubeBank: raise TypeError('Unsupported bank class.')
            walls.append((wall.wall_capacity_j_k,wall.air_inlet_temperature_k,wall._effective_air_conductance,
                film.half_wall_resistance_k_w,b.tube_internal_area_m2,b.inner_diameter_m,b.tube_length_m,
                b.tube_flow_area_m2,g.maximum_mach,g.maximum_relative_pressure_drop,float(g.thermal_entry),
                tr.minimum_temperature,tr.maximum_temperature))
        self.wrapper=wrapper;self.reference=reference
        self.links=np.array(links);self.walls=np.array(walls)
        g=model.gas
        self.gas=np.array((g.gas_constant,g.heat_capacity_cv,g.heat_capacity_cp,g.heat_capacity_ratio,model.angular_speed))
        self.calls=0;self.fallbacks=0
        self.first_call_seconds=None
        for array in (self.links,self.walls,self.gas): array.flags.writeable=False

    def __call__(self,angle,values):
        self.calls+=1
        values=np.asarray(values,dtype=float)
        if values.ndim!=1 or values.size<10:
            self.fallbacks+=1
            return self.reference(self.wrapper,angle,values)
        started=time.perf_counter() if self.first_call_seconds is None else None
        model=self.wrapper.model;k=model.kinematics
        provider=getattr(k,'cylinder_volumes_and_derivatives',None)
        if provider is not None:
            vs,vl,ds,dl=provider(angle)
        else:
            vs,vl=k.small_cylinder_volume(angle),k.large_cylinder_volume(angle)
            ds,dl=k.small_cylinder_volume_derivative(angle),k.large_cylinder_volume_derivative(angle)
        volumes=np.array((vs,vl,model.machine_volumes.cold_heat_exchanger,model.machine_volumes.hot_heat_exchanger))
        rates=np.array((model.angular_speed*ds,model.angular_speed*dl))
        ok,result=kernel(np.asarray(values,dtype=float),volumes,rates,self.gas,self.links,self.walls)
        if started is not None: self.first_call_seconds=time.perf_counter()-started
        if ok: return result
        self.fallbacks+=1
        return self.reference(self.wrapper,angle,values)


@contextmanager
def prototype_backend():
    """Process-local benchmark scope only; never used by persistent campaigns."""
    from dada_solver.exchangers.air_wall import AirWallMotor
    reference=AirWallMotor.derivative
    prepared={};stats=dict(candidates=[],unsupported_calls=0)
    def derivative(wrapper,angle,values):
        key=id(wrapper)
        if key not in prepared:
            try: obj=PreparedRHS(wrapper,reference)
            except TypeError: obj=None
            prepared[key]=(wrapper,obj)  # Retain owner: object IDs cannot be reused.
            stats['candidates'].append(dict(supported=obj is not None))
        obj=prepared[key][1]
        if obj is None:
            stats['unsupported_calls']+=1
            return reference(wrapper,angle,values)
        return obj(angle,values)
    try:
        with patch.object(AirWallMotor,'derivative',derivative):
            yield stats
    finally:
        stats['candidates']=[dict(supported=obj is not None,**(dict(calls=obj.calls,
            fallback_calls=obj.fallbacks,first_call_seconds=obj.first_call_seconds) if obj else {}))
            for _,obj in prepared.values()]
