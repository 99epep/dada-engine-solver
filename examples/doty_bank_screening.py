"""Write illustrative Doty bank scenarios; these are not motor predictions.

Run from the checkout: PYTHONPATH=src python3 examples/doty_bank_screening.py
"""
import csv
from pathlib import Path

from dada_solver.exchangers.doty import load_references, scale_bank

ROOT = Path(__file__).resolve().parents[1]
reference = load_references(ROOT / "examples/data/doty_1991_nitrogen_reference.csv")[1]
rows = []
for banks in (1, 4, 16):
    for duty in (1.0, 0.5, 0.25):
        for exponent in (0.0, 0.5):
            result = scale_bank(
                reference, parallel_banks=banks,
                total_mass_flow_kg_s=banks * reference.mass_flow_kg_s / duty,
                ua_flow_exponent=exponent, ua_condition_factor=1.0,
                viscosity_ratio=1.0, density_ratio=1.0,
                additional_pressure_drop_pa=0.0,
            )
            rows.append(dict(rectangular_pulse_duty=duty, ua_flow_exponent=exponent,
                             **result))
output = ROOT / "outputs/doty_bank_screening.csv"
output.parent.mkdir(exist_ok=True)
with output.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(output)
