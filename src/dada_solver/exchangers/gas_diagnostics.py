"""Time- and heat-weighted instantaneous microtube domain reports."""
import math
import numpy as np
from dada_solver.state import ThermodynamicState


def cycle_microtube_diagnostics(wrapper, angles, trajectory):
    """Report actual variable-film domains; no static-UA trajectory replay.

    Half-tube heat shares are allocated in proportion to their film conductance.
    Fractions use trapezoidal physical-time weighting, never sample counts.
    Hydraulic port temperatures are selected by the actual signed flow.
    """
    if not any(getattr(w,'requires_flow_context',False) for w in (wrapper.heat_in,wrapper.heat_out)):
        return None
    samples={name:[] for name in ('Hi.inlet','Hi.outlet','Ho.inlet','Ho.outlet')}
    weights={name:[] for name in samples}
    hydraulics={name:[] for name in samples}
    for angle,values in zip(angles,trajectory.T):
        contexts=wrapper.flow_contexts(float(angle),values)
        temperatures=ThermodynamicState.from_array(values[:8]).temperatures(wrapper.model.gas)
        thermal=wrapper.thermal_rates(float(angle),values)
        for side,index,wall,context,rates,ports in zip(('Hi','Ho'),(2,3),(wrapper.heat_in,wrapper.heat_out),
                contexts,thermal,(((0,2),(2,1)),((1,3),(3,0)))):
            if not getattr(wall,'requires_flow_context',False): continue
            film=wall.gas_film
            _,diagnostics=film.evaluate(temperatures[index],values[index+6]/wall.wall_capacity_j_k,context)
            total_nu=sum(d.nusselt for d in diagnostics)
            for label,diagnostic,passage,pair in zip(('inlet','outlet'),diagnostics,context['passages'],ports):
                name=f'{side}.{label}'
                samples[name].append(diagnostic)
                weights[name].append(abs(rates['gas_heat_w'])*diagnostic.nusselt/total_nu)
                flow,p1,p2=passage
                upstream=pair[0] if flow>=0 else pair[1]
                hydraulic = film.model.diagnose(film.bank,flow,p1,p2,temperatures[upstream],
                    frequency=context['frequency_hz'],length=film.bank.tube_length_m/2)
                hydraulics[name].append(hydraulic)
    time=np.asarray(angles)/wrapper.model.angular_speed
    integrate=lambda y: float(np.trapz(np.asarray(y,dtype=float),time))
    duration=time[-1]-time[0]; result={}; failures=set()
    fields=('reynolds','prandtl','mach','knudsen','mean_knudsen','graetz','pressure_ratio',
            'compressibility_parameter','womersley','strouhal','viscous_diffusion_time_s',
            'thermal_diffusion_time_s','residence_time_s','acoustic_time_s')
    for name,rows in samples.items():
        if not rows: continue
        heat=np.asarray(weights[name]); total_heat=integrate(heat)
        ranges={}
        for key in fields:
            values=[getattr(r,key) for r in rows if getattr(r,key) is not None]
            ranges[key]=dict(minimum=min(values),maximum=max(values)) if values else None
        domains={}
        for field in ('correlation_id','continuum_regime','slip_regime','model_validity',
                      'thermal_developing','hydrodynamic_developing','compressibility_significant'):
            domains[field]={}
            for value in sorted(set(getattr(r,field) for r in rows),key=str):
                indicator=np.array([getattr(r,field)==value for r in rows])
                domains[field][str(value)]=dict(time_fraction=integrate(indicator)/duration,
                    absolute_heat_fraction=integrate(indicator*heat)/total_heat if total_heat else None)
        for row in rows:
            failures.update(name+':'+issue for issue in row.issues)
        for row in hydraulics[name]:
            failures.update(name+':hydraulic_'+issue for issue in row.issues
                            if issue in ('high_mach','beyond_continuum_model','reynolds_outside_correlation_domain'))
        hydraulic_ranges={key:dict(minimum=min(getattr(r,key) for r in hydraulics[name]),
                                   maximum=max(getattr(r,key) for r in hydraulics[name]))
                          for key in ('reynolds','mach','knudsen','pressure_ratio')}
        result[name]=dict(ranges=ranges,domains=domains,hydraulic_upstream_ranges=hydraulic_ranges,
                         absolute_gas_heat_J=total_heat)
    return dict(passages=result,model_validity='invalid' if failures else 'valid',
        failed_criteria=sorted(failures),time_weighting='trapezoidal physical time',
        heat_weighting='absolute wall-to-gas heat, split by instantaneous half-film conductance',
        limitations=['quasi_steady_pulse_response_unvalidated','stagnant_radial_heat_is_lumped_screening',
                     'no_axial_wall_or_gas_temperature_resolution','external_air_film_unchanged'])
