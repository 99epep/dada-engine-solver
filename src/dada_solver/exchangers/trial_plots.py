"""Plots and sampled passive-valve event diagnostics for air/wall motor trials."""
from pathlib import Path
import math
import numpy as np
from dada_solver.state import ThermodynamicState


def save_trial_plots(wrapper, angles, trajectory, ports, prefix, *, external_losses_excluded=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    temperatures=[];pressures=[];volumes=[];heats=[]
    for angle,state in zip(angles,trajectory.T):
        gas=ThermodynamicState.from_array(state[:8])
        volume=wrapper.model.volumes(float(angle))
        t,p=gas.temperatures_and_pressures(wrapper.model.gas,volume)
        temperatures.append(t);pressures.append(p)
        volumes.append([volume.small_cylinder,volume.large_cylinder])
        hi,ho=wrapper.thermal_rates(float(angle),state)
        heats.append([hi['gas_heat_w'],ho['gas_heat_w'],hi['air_heat_w'],ho['air_heat_w']])
    temperatures=np.asarray(temperatures);pressures=np.asarray(pressures)
    volumes=np.asarray(volumes);heats=np.asarray(heats)
    degrees=np.degrees(angles)
    events=[]
    for label,difference in [('H_i_to_L',pressures[:,2]-pressures[:,1]),('H_o_to_S',pressures[:,3]-pressures[:,0])]:
        for i in range(len(angles)-1):
            if angles[i+1]<=angles[i] or (difference[i]>0)==(difference[i+1]>0):continue
            fraction=-difference[i]/(difference[i+1]-difference[i])
            angle=float(angles[i]+fraction*(angles[i+1]-angles[i]))
            volume=wrapper.model.volumes(angle)
            events.append(dict(valve=label,transition='opening' if difference[i+1]>0 else 'closing',
                cycle_progress_degrees=math.degrees(angle),
                lambda_small=wrapper.model.machine_volumes.small_cylinder.closure_fraction(volume.small_cylinder),
                lambda_large=wrapper.model.machine_volumes.large_cylinder.closure_fraction(volume.large_cylinder),
                method='linear_interpolation_of_sampled_pressure_zero_crossing'))
    fig,ax=plt.subplots(3,2,figsize=(13,11))
    for i,label in enumerate(['S','L','H_i','H_o']):
        ax[0,0].plot(degrees,temperatures[:,i]-273.15,label=label)
    ax[0,0].plot(degrees,trajectory[8]/wrapper.heat_in.wall_capacity_j_k-273.15,'--',label='H_i wall')
    ax[0,0].plot(degrees,trajectory[9]/wrapper.heat_out.wall_capacity_j_k-273.15,'--',label='H_o wall')
    ax[0,0].axhline(wrapper.heat_in.air_inlet_temperature_k-273.15,color='red',ls=':',label='Hot air inlet')
    ax[0,0].axhline(wrapper.heat_out.air_inlet_temperature_k-273.15,color='blue',ls=':',label='Cold air inlet')
    ax[0,0].set_ylabel('Temperature (deg C)')
    ax[0,1].plot(degrees,wrapper.heat_in.air_inlet_temperature_k-temperatures[:,2],label='Hot inlet air - H_i gas')
    ax[0,1].plot(degrees,temperatures[:,3]-wrapper.heat_out.air_inlet_temperature_k,label='H_o gas - cold inlet air')
    ax[0,1].set_ylabel('Inlet-air / working-gas difference (K)')
    for i,label in enumerate(['S','L','H_i','H_o']):
        ax[1,0].plot(degrees,pressures[:,i]/1e5,label=label)
    for i,label in enumerate(['S','L']):
        ax[1,1].plot(degrees,volumes[:,i]*1000,label=label)
    ax[1,0].set_ylabel('Pressure (bar absolute)')
    ax[1,1].set_ylabel('Cylinder volume (L)')
    for i,label in enumerate(['S to H_i','H_i to L','L to H_o','H_o to S']):
        ax[2,0].plot(degrees,ports[:,i]*1000,label=label)
    ax[2,0].set_ylabel('Signed port mass flow (g/s)')
    for i,label in enumerate(['H_i gas','H_o gas','Hot external air','Cold external air']):
        ax[2,1].plot(degrees,heats[:,i],label=label)
    ax[2,1].set_ylabel('Heat received by gas or wall (W)')
    for a in ax.flat:
        a.grid(alpha=.3);a.legend(fontsize=8)
        a.set_xlabel('Motor cycle progress (degrees)')
        for event_angle in sorted({round(e['cycle_progress_degrees'],6) for e in events}):
            a.axvline(event_angle,color='0.4',ls='--',lw=.55,alpha=.7,label='_nolegend_')
    fig.suptitle('Coupled four-bar trial: screening assumptions' + ('; external aerodynamic losses excluded' if external_losses_excluded else ''))
    fig.tight_layout();fig.savefig(Path(str(prefix)+'_graphs.png'),dpi=150);plt.close(fig)
    pv,panels=plt.subplots(1,2,figsize=(11,4.5))
    for i,(panel,label) in enumerate(zip(panels,['S','L'])):
        panel.plot(volumes[:,i]*1000,pressures[:,i]/1e5)
        work=float(np.trapz(pressures[:,i],volumes[:,i]))
        panel.set_title(f'{label}: signed gas work ~ {work:.2f} J/cycle')
        for phase in (45,165,285):
            x0=np.interp(phase,degrees,volumes[:,i])*1000
            x1=np.interp(phase+10,degrees,volumes[:,i])*1000
            y0=np.interp(phase,degrees,pressures[:,i])/1e5
            y1=np.interp(phase+10,degrees,pressures[:,i])/1e5
            panel.annotate('',xy=(x1,y1),xytext=(x0,y0),arrowprops=dict(arrowstyle='->',color='black',lw=1))
        for event in events:
            phase=event['cycle_progress_degrees']
            panel.plot(np.interp(phase,degrees,volumes[:,i])*1000,
                       np.interp(phase,degrees,pressures[:,i])/1e5,'.',color='black',ms=5)
        panel.set_xlabel('Cylinder volume (L)');panel.set_ylabel('Pressure (bar absolute)');panel.grid(alpha=.3)
    pv.suptitle('Pressure-volume cycles: arrows show time direction; dots mark valve events')
    pv.tight_layout();pv.savefig(Path(str(prefix)+'_pv.png'),dpi=150);plt.close(pv)
    return events
