"""Produce transparent Doty hydraulic errors and Graur fit-reproduction checks."""
from dataclasses import asdict
import json
from pathlib import Path
from dada_solver.exchangers.doty import validate_gas_hydraulics
from dada_solver.exchangers.gas_correlations import graur_silica_fit

ROOT=Path(__file__).resolve().parents[1]
rows=[]
for gas in ('nitrogen','helium'):
    rows.extend(validate_gas_hydraulics(ROOT/f'examples/data/doty_1991_{gas}_reference.csv'))
graur=[]
for gas,a,b,sa,sb,tmac,stmac in [('nitrogen',11.668,16.626,.967,4.059,.908,.041),
                                ('argon',13.218,24.274,.799,3.698,.871,.017),
                                ('helium',10.812,9.156,.371,1.500,.914,.009)]:
    model=graur_silica_fit(gas)
    for kn in (.01,.03,.1,.2):
        expected=1+a*kn+b*kn*kn;predicted=model.factor(kn,5)
        graur.append(dict(gas=gas,surface='fused_silica',knudsen_convention='graur_vhs',
            mean_knudsen=kn,pressure_ratio=5,measured=None,published_fit=expected,predicted=predicted,
            absolute_error_against_fit=abs(predicted-expected),relative_error_against_fit=(predicted-expected)/expected,
            coefficient_uncertainty_envelope=sa*kn+sb*kn*kn,
            inside_experimental_uncertainty=None,reference_kind='published_regression_not_raw_measurement',
            model_validity='within_published_fit_range_only',source='HAL_microtubes.pdf Table 1 and Eq. 3,5',
            accommodation_reference=dict(momentum=tmac,uncertainty=stmac,applicable_to_DADA_metal=False),
            slip_coefficients=asdict(model)))
report=dict(doty=rows,graur=graur,calibration='No fitted multiplier; no metal TMAC inferred',
            uncertainty_note='Doty output uncertainty only; missing geometric/temperature uncertainties are not invented.')
(ROOT/'outputs/microtube_gas_validation.json').write_text(json.dumps(report,indent=2)+'\n')
for row in rows:
    print(row['fluid'],row['row'],row['model_id'],round(row['predicted']),round(100*row['relative_error'],2),row['inside_experimental_uncertainty'])
