"""Audit paired temperature changes without assuming a fitted loss mechanism."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
rows = []
with (ROOT/'examples/data/doty_1991_nitrogen_reference.csv').open() as stream:
    for row in csv.DictReader(stream):
        t1, t2, t3, t4 = (float(row[key]) for key in ('T1_K','T2_K','T3_K','T4_K'))
        hot, cold = t1-t2, t4-t3
        # Explicit approximate nitrogen cp; Eq. 19 uses mean stream temperature difference.
        cp = 1040.0
        mean_difference = (t1+t2-t3-t4)/2
        reduced_ua = float(row['mass_flow_kg_s'])*cp*hot/mean_difference
        rows.append(dict(mass_flow_kg_s=row['mass_flow_kg_s'],
                         hot_temperature_drop_k=hot, cold_temperature_rise_k=cold,
                         assumed_nitrogen_cp_j_kg_k=cp,
                         reduced_apparent_ua_w_k=reduced_ua,
                         reported_ua_w_k=float(row['UA_W_K']),
                         within_reported_ua_uncertainty=abs(reduced_ua-float(row['UA_W_K'])) <= float(row['UA_uncertainty_W_K']),
                         heat_imbalance_fraction_of_hot=(hot-cold)/hot,
                         assumption='equal_mass_flows_and_equal_constant_cp'))
output = ROOT/'outputs/doty_energy_balance_audit.csv'
output.parent.mkdir(exist_ok=True)
with output.open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(output)
