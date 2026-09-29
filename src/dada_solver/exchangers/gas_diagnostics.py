"""Time- and heat-weighted instantaneous microtube domain reports."""
import math
import numpy as np
from scipy.integrate import trapezoid
from dada_solver.state import ThermodynamicState


def cycle_microtube_diagnostics(wrapper, angles, trajectory, *, replay=None):
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
    if replay is not None:
        replay.require(wrapper=wrapper, angles=angles, trajectory=trajectory)
        for sample in replay.samples:
            for side,wall_index,wall in zip(('Hi','Ho'),(0,1),sample.walls):
                if not wall.film_diagnostics: continue
                total_nu=sum(d.nusselt for d in wall.film_diagnostics)
                for port,(label,diagnostic) in enumerate(zip(('inlet','outlet'),wall.film_diagnostics)):
                    name=f'{side}.{label}'
                    samples[name].append(diagnostic)
                    weights[name].append(abs(wall.gas_heat_w)*diagnostic.nusselt/total_nu)
                    hydraulics[name].append(sample.hydraulic_diagnostics[2*wall_index+port])
    else:
        for angle,values in zip(angles,trajectory.T):
            contexts=wrapper.flow_contexts(float(angle),values)
            from dada_solver.fluids import CaloricallyPerfectGas
            volumes=None if type(wrapper.model.gas) is CaloricallyPerfectGas else wrapper.model.volumes(float(angle))
            temperatures=ThermodynamicState.from_array(values[:8]).temperatures(wrapper.model.gas, volumes)
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
    integrate=lambda y: float(trapezoid(np.asarray(y,dtype=float),time))
    duration=time[-1]-time[0]; result={}; failures=set()
    transition_any=np.zeros(len(time),dtype=bool)
    transition_heat=0.; all_heat=0.
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
        for field in ('flow_regime','correlation_id','continuum_regime','slip_regime','model_validity',
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
        transition=np.array([r.flow_regime=='transition' or h.flow_regime=='transition'
                             for r,h in zip(rows,hydraulics[name])])
        transition_re=[r.reynolds for r in (*rows,*hydraulics[name]) if r.flow_regime=='transition']
        transition_any |= transition
        heat_in_transition=integrate(transition*heat)
        transition_heat+=heat_in_transition; all_heat+=total_heat
        result[name]=dict(ranges=ranges,domains=domains,hydraulic_upstream_ranges=hydraulic_ranges,
            absolute_gas_heat_J=total_heat, transition_model_used=bool(np.any(transition)),
            transition_time_fraction=integrate(transition)/duration,
            transition_absolute_heat_fraction=heat_in_transition/total_heat if total_heat else None,
            transition_reynolds_range=dict(minimum=min(transition_re),maximum=max(transition_re)) if transition_re else None)
    return dict(passages=result,model_validity='invalid' if failures else 'valid',
        transition_model_used=bool(np.any(transition_any)),
        confidence='transition_uncertainty' if np.any(transition_any) else 'standard_correlation_scope',
        transition_time_fraction=integrate(transition_any)/duration,
        transition_absolute_heat_fraction=transition_heat/all_heat if all_heat else None,
        failed_criteria=sorted(failures),time_weighting='trapezoidal physical time',
        heat_weighting='absolute wall-to-gas heat, split by instantaneous half-film conductance',
        limitations=['quasi_steady_pulse_response_unvalidated','stagnant_radial_heat_is_lumped_screening',
                     'no_axial_wall_or_gas_temperature_resolution','external_air_film_unchanged'])
