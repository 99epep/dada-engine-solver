"""Auditable circular-tube gas correlations, without empirical pulse multipliers."""
from dataclasses import dataclass, field
import math
from .gas_transport import DiluteGasTransport, GasTransportModel


class MicrotubeDomainError(ValueError):
    """A required physical closure is outside its declared domain."""


@dataclass(frozen=True)
class GasSurfaceAccommodation:
    momentum_accommodation: float | None = None
    thermal_accommodation: float | None = None
    provenance: str = 'Unknown gas/metal accommodation; no silica value inferred'

    def __post_init__(self):
        for v in (self.momentum_accommodation,self.thermal_accommodation):
            if v is not None and (not math.isfinite(v) or not 0<v<=1):
                raise ValueError('Accommodation must be None or in (0, 1].')
        if not self.provenance: raise ValueError('Accommodation provenance is required.')


@dataclass(frozen=True)
class SecondOrderSlip:
    a1: float
    a2: float
    provenance: str
    mean_free_path_convention: str
    maximum_mean_knudsen: float = .1
    maximum_pressure_ratio: float = 5.

    def __post_init__(self):
        if any(not math.isfinite(v) or v<0 for v in (self.a1,self.a2)) or not self.provenance:
            raise ValueError('Explicit nonnegative slip coefficients and provenance required.')
        if not 0<self.maximum_mean_knudsen<=.31 or self.maximum_pressure_ratio<1:
            raise ValueError('Invalid slip domain.')

    def factor(self,kn_mean,pressure_ratio):
        if not (0<=kn_mean<=self.maximum_mean_knudsen and 1<=pressure_ratio<=self.maximum_pressure_ratio):
            raise MicrotubeDomainError('Slip correlation outside declared Kn/pressure-ratio domain.')
        x=pressure_ratio-1
        ratio_log=2. if x==0 else (2+x)*math.log1p(x)/x
        return 1+8*self.a1*kn_mean+16*self.a2*ratio_log*kn_mean**2


def graur_silica_fit(species):
    """Validation-only coefficients in the paper's VHS Kn convention, P=5.

    These are NOT metal coefficients and cannot be attached to the HS transport
    path without an explicit, independently justified convention conversion.
    """
    a,b,limit={'nitrogen':(11.668,16.626,.291),'argon':(13.218,24.274,.302),
               'helium':(10.812,9.156,.309)}[species]
    return SecondOrderSlip(a/8,b/(16*1.5*math.log(5)),
        f'Ewart et al. CFM2007 Tables 1-2; {species}/fused silica, P=5',
        'graur_vhs',limit,5.)


def knudsen_regime(kn):
    if not math.isfinite(kn) or kn<0: raise ValueError('Kn must be finite and nonnegative.')
    return ('continuum' if kn<.001 else 'slip_onset' if kn<.01 else
            'slip' if kn<=.1 else 'beyond_continuum_model')


def compressible_poiseuille(p1,p2,t,mu,length,diameter,count,gas_constant):
    """Signed isothermal no-slip mass flow; N*pi*D^4*(p1^2-p2^2)/(256*mu*L*R*T)."""
    if any(not math.isfinite(x) or x<=0 for x in (p1,p2,t,mu,length,diameter,count,gas_constant)):
        raise ValueError('Positive finite tube and gas inputs required.')
    return count*math.pi*diameter**4*(p1-p2)*(p1+p2)/(256*mu*length*gas_constant*t)


def laminar_entry_nusselt(graetz):
    """Hausen mean Nu, constant wall temperature, developed velocity profile."""
    if not math.isfinite(graetz) or graetz<0: raise ValueError('Graetz must be nonnegative.')
    return 3.66+.0668*graetz/(1+.04*graetz**(2/3))


def darcy_smooth(reynolds):
    """Reuse the existing Haaland Darcy closure for a hydraulically smooth tube."""
    from .models import _darcy_friction_factor, CorrelationRegime
    if not 4000<=reynolds<=5e6: raise MicrotubeDomainError('Turbulent Reynolds outside 4000-5e6.')
    return _darcy_friction_factor(reynolds,1.,0.,CorrelationRegime.TURBULENT)


def gnielinski(reynolds,prandtl):
    if not .5<=prandtl<=2000: raise MicrotubeDomainError('Gnielinski Pr outside 0.5-2000.')
    f=darcy_smooth(reynolds)
    from .models import _nusselt_number, CorrelationRegime
    return _nusselt_number(reynolds,prandtl,1.,f,CorrelationRegime.TURBULENT)


@dataclass(frozen=True)
class MicrotubeFlowDiagnostics:
    reynolds: float
    prandtl: float
    mach: float
    knudsen: float
    mean_knudsen: float
    graetz: float
    pressure_ratio: float
    compressibility_parameter: float
    thermal_developing: bool
    hydrodynamic_developing: bool
    continuum_regime: str
    slip_regime: str
    compressibility_significant: bool
    correlation_id: str
    nusselt: float | None
    womersley: float
    strouhal: float | None
    viscous_diffusion_time_s: float
    thermal_diffusion_time_s: float
    residence_time_s: float | None
    acoustic_time_s: float
    model_validity: str
    issues: tuple[str,...]


@dataclass(frozen=True)
class MicrotubeGasModel:
    """Production internal-gas settings; None in HardwareInputs selects legacy.

    Domain policy 'report' explicitly permits exploratory extrapolation while
    retaining an INVALID verdict. The default rejects unsupported closures.
    """
    transport: GasTransportModel = field(default_factory=DiluteGasTransport)
    accommodation: GasSurfaceAccommodation = field(default_factory=GasSurfaceAccommodation)
    slip: SecondOrderSlip | None = None
    domain_policy: str = 'reject'
    thermal_entry: bool = True
    maximum_mach: float = .3
    maximum_relative_pressure_drop: float = .2
    boundary_condition: str = 'constant_wall_temperature'

    def __post_init__(self):
        if self.domain_policy not in ('reject','report'): raise ValueError('Unknown domain policy.')
        if self.boundary_condition!='constant_wall_temperature':
            raise ValueError('Only the existing isothermal lumped wall boundary is implemented.')
        if not 0<self.maximum_mach<=.3 or not 0<self.maximum_relative_pressure_drop<=.2:
            raise ValueError('Compressibility guards cannot exceed the declared low-Mach domain.')
        if self.slip is not None and self.slip.mean_free_path_convention!=getattr(self.transport,'mean_free_path_convention',None):
            raise ValueError('Slip and transport mean-free-path conventions differ.')

    def slip_factor(self,p1,p2,t,diameter):
        mean=(p1+p2)/2
        kn=self.transport.mean_free_path(mean,t)/diameter
        if kn<.001: return 1.
        if self.slip is None:
            if self.domain_policy=='reject': raise MicrotubeDomainError('Momentum slip required but gas/surface coefficients are unknown.')
            return 1.
        return self.slip.factor(kn,max(p1,p2)/min(p1,p2))

    def diagnose(self,bank,flow,p1,p2,t,*,frequency=0.,length=None):
        length=bank.tube_length_m if length is None else length
        if not all(math.isfinite(v) and v>0 for v in (p1,p2,t,length)) or not math.isfinite(flow) or frequency<0:
            raise ValueError('Invalid instantaneous state.')
        tr=self.transport; mu=tr.viscosity(t); conductivity=tr.conductivity(t);cp=tr.cp(t)
        r=tr.gas_constant;rho=min(p1,p2)/(r*t);pm=(p1+p2)/2
        area=bank.tube_count*math.pi*bank.inner_diameter_m**2/4;d=bank.inner_diameter_m
        u=abs(flow)/(rho*area);re=abs(flow)*d/(area*mu);pr=cp*mu/conductivity
        speed=math.sqrt(cp/(cp-r)*r*t);ma=u/speed
        kn=tr.mean_free_path(min(p1,p2),t)/d;knm=tr.mean_free_path(pm,t)/d
        gz=re*pr*d/length;ratio=max(p1,p2)/min(p1,p2);compressibility=2*(ratio-1)/(ratio+1)
        issues=[]
        if ma>self.maximum_mach: issues.append('high_mach')
        if kn>=.001: issues.append('thermal_slip_not_implemented')
        if kn>.1: issues.append('beyond_continuum_model')
        if compressibility>self.maximum_relative_pressure_drop: issues.append('large_relative_pressure_drop')
        if not .5<=pr<=2000: issues.append('prandtl_outside_domain')
        if re<2300:
            nu=laminar_entry_nusselt(gz) if self.thermal_entry else 3.66
            correlation='hausen_constant_wall' if self.thermal_entry else 'fully_developed_3_66_screening'
            if re==0: correlation='stagnant_radial_screening'
            if length<.05*re*d: issues.append('hydrodynamic_entry_unresolved')
        elif re>=4000 and re<=5e6 and .5<=pr<=2000:
            nu=gnielinski(re,pr);correlation='gnielinski_smooth'
            if length<10*d: issues.append('turbulent_entry_unresolved')
        else:
            nu=None;correlation='unavailable_transition_or_out_of_range'
            issues.append('reynolds_outside_correlation_domain')
        viscous=(d/2)**2*rho/mu;thermal=viscous*pr
        return MicrotubeFlowDiagnostics(re,pr,ma,kn,knm,gz,ratio,compressibility,
            length<.05*re*pr*d,length<.05*re*d,knudsen_regime(kn),
            'negligible' if knm<.001 else 'second_order' if self.slip else 'unknown_accommodation',
            ma>self.maximum_mach or compressibility>self.maximum_relative_pressure_drop,
            correlation,nu,math.sqrt(2*math.pi*frequency*viscous),
            frequency*length/u if u else None,viscous,thermal,length/u if u else None,
            length/speed,'invalid' if issues else 'valid',tuple(issues))

    def require(self,diagnostics):
        if diagnostics.issues and self.domain_policy=='reject':
            raise MicrotubeDomainError('; '.join(diagnostics.issues))
        if diagnostics.nusselt is None:
            raise MicrotubeDomainError('No thermal closure for this Reynolds/Pr domain.')
