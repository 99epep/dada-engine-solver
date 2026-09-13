"""Plot the completed K2 harmonic phase scan and export its compact table."""
import csv
import json
import os
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def main():
    os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'dada_solver_matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report=json.loads((ROOT/'outputs/motor_harmonic_k2_fine_scan.json').read_text())
    rows=[r for r in report['results'] if r['status']=='converged']
    fields=['phase_degrees','indicated_thermal_efficiency','indicated_power_w','heat_input_w','elapsed_seconds']
    with (ROOT/'outputs/motor_harmonic_k2_fine_scan.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        writer.writerows({k:r[k] for k in fields} for r in rows)
    fig,axes=plt.subplots(2,1,figsize=(9,7),sharex=True)
    phase=[r['phase_degrees'] for r in rows]
    axes[0].plot(phase,[100*r['indicated_thermal_efficiency'] for r in rows],'o-')
    axes[1].plot(phase,[r['indicated_power_w'] for r in rows],'o-')
    coupler=json.loads((ROOT/'outputs/compact_motor_coupler_k2.json').read_text())['result']
    axes[0].axhline(100*coupler['indicated_thermal_efficiency'],color='tab:orange',ls='--',label='Finite-rod coupler candidate, K2')
    axes[1].axhline(coupler['indicated_power_w'],color='tab:orange',ls='--')
    axes[0].set_ylabel('Indicated thermal efficiency (%)');axes[0].legend()
    axes[1].set_ylabel('Indicated gas power (W)');axes[1].set_xlabel('Small-cylinder phase offset (deg, solver convention)')
    for ax in axes: ax.grid(alpha=.25)
    fig.suptitle('Ideal harmonic motion — K2 hardware, 0.5° phase scan')
    fig.tight_layout();fig.savefig(ROOT/'outputs/motor_harmonic_k2_fine_scan.png',dpi=140)

if __name__=='__main__':main()
