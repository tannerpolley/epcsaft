"""Reproduce #214 and the common-pressure scalar difference with the installed wheel."""
import csv
import hashlib
import json
import math
from pathlib import Path

from epcsaft import Mixture, Parameters

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'data/packets/held-cameretti-sadowski-2008/3/parameters/parameters.json'
model = Mixture(Parameters.from_json(INPUT, components=['water', 'sodium-cation', 'chloride-anion']))
T, P, MW = 298.15, 100000.0, 0.01801528
reference = model.reference(T, P, solvent_fractions=[1, 0, 0],
                            ions=('sodium-cation', 'chloride-anion'), stoichiometry=(1, 1))
pure = model.state(T, P=P, x=[1, 0, 0], phase='liquid')
rows = []
for m in (0.001, 1e-6, 1e-8, 1e-10, 1e-12):
    total = 1 / MW + 2 * m
    x = [1 / MW / total, m / total, m / total]
    finite = model.state(T, P=P, x=x, phase='liquid')
    mu = finite.residual_chemical_potential_over_rt[0] - pure.residual_chemical_potential_over_rt[0]
    log_rho = math.log1p((finite.molar_density - pure.molar_density) / pure.molar_density)
    log_activity = math.fsum((-math.log1p(2 * m * MW), mu, log_rho))
    rows.append(dict(molality_mol_kg=m, osmotic_native=reference.activity(m).osmotic_coefficient,
                     osmotic_common_pressure=-log_activity / (2 * m * MW),
                     finite_pressure_residual_Pa=finite.pressure-P,
                     pure_pressure_residual_Pa=pure.pressure-P,
                     residual_chemical_potential_difference_over_RT=mu, log_density_ratio=log_rho))
with (HERE / 'results/diagnosis.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=rows[0])
    writer.writeheader()
    writer.writerows(rows)
wheel, = (ROOT / 'build/environment-wheel').glob('*.whl')
(HERE / 'results/diagnosis-identity.json').write_text(json.dumps({
    'wheel_sha256': hashlib.sha256(wheel.read_bytes()).hexdigest(),
    'input_path': str(INPUT.relative_to(ROOT)), 'input_sha256': hashlib.sha256(INPUT.read_bytes()).hexdigest(),
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'source_base': '3f5ee19f', 'temperature_K': T, 'pressure_Pa': P,
}, indent=2) + '\n')
print(rows[2])
