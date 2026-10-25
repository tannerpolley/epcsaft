"""Compare the installed scalar activity with retained chloride values and the Debye limit."""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

from epcsaft import Mixture, Parameters, Stop, provenance

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
OUTPUT = HERE / 'results/current-main'
OUTPUT.mkdir(parents=True, exist_ok=True)
PRIOR = ROOT / 'analyses/2026-divalent-chloride-activity'
T, P, MW = 298.15, 100000.0, 0.01801528
SALTS = {'LiCl': ('lithium-cation', 1), 'NaCl': ('sodium-cation', 1),
         'MgCl2': ('magnesium-cation', 2), 'CaCl2': ('calcium-cation', 2)}
models = {(mode, salt): Mixture(Parameters.from_json(
    ROOT / f'data/packets/held-cameretti-sadowski-2008/{mode[-1]}/parameters/parameters.json',
    components=['water', cation, 'chloride-anion']))
    for mode in ('packet3', 'packet4') for salt, (cation, _) in SALTS.items()}
mixed_input = ROOT / 'data/packets/figiel-2025-reference-electrolytes/4/parameters/parameters.json'
mixed_mapping = json.loads(mixed_input.read_text())
# Packet 4 has no water/methanol cross edge; declared to keep that model (#569), not a physical-absence claim.
mixed_mapping['topology']['no_cross_association'] = [{
    'component_id_a': 'water', 'component_id_b': 'methanol',
    'provenance': {'source_id': mixed_mapping['document_id'], 'locator': 'packet 4: no water/methanol edge'}}]
mixed = Mixture(Parameters.from_mapping(mixed_mapping, components=['water', 'methanol', 'sodium-cation', 'chloride-anion']))


def declarations():
    return models, mixed


if '--declare-only' in sys.argv:
    raise SystemExit(0)
retained = list(csv.DictReader((PRIOR / 'results/model-curves.csv').open()))
prior_pure_pressure = {(r['mode'], r['salt']): P + float(r['pure_pressure_residual_Pa']) for r in
                       csv.DictReader((PRIOR / 'results/solvent-activity-conditioning.csv').open())}
current = {(r['mode'], r['salt'], float(r['molality_mol_kg'])): r for r in
           csv.DictReader((PRIOR / 'results/current-main/model-curves.csv').open())}
rows = []
for mode in ('packet3', 'packet4'):
    for salt, (cation, nu_minus) in SALTS.items():
        path = ROOT / f'data/packets/held-cameretti-sadowski-2008/{mode[-1]}/parameters/parameters.json'
        model = models[mode, salt]
        reference = model.reference(T, P, solvent_fractions=[1, 0, 0],
                                    ions=(cation, 'chloride-anion'), stoichiometry=(1, nu_minus))
        pure = model.state(T, P=P, x=[1, 0, 0], phase='liquid')
        for old in (r for r in retained if r['mode'] == mode and r['salt'] == salt):
            m = float(old['molality_mol_kg'])
            s = (1 + nu_minus) * m * MW
            total = 1 / MW + (1 + nu_minus) * m
            x = [1 / MW / total, m / total, nu_minus * m / total]
            state = model.state(T, P=P, x=x, phase='liquid')
            scalar = -math.fsum((-math.log1p(s),
                state.residual_chemical_potential_over_rt[0] - pure.residual_chemical_potential_over_rt[0],
                math.log1p((state.molar_density - pure.molar_density) / pure.molar_density))) / s
            activity = reference.activity(m)
            assert isinstance(activity.osmotic_coefficient, float)
            solvent = reference.solvent_activity(x=x)
            assert abs(activity.osmotic_coefficient - solvent.osmotic_coefficient) <= 1e-12
            old_phi = float(old['osmotic_coefficient'])
            delta = activity.osmotic_coefficient - old_phi
            # The removed pressure error belongs to the historical roots, not the current roots.
            prior_pressure = P + float(old['pressure_residual_Pa'])
            prior_pure = prior_pure_pressure[mode, salt]
            pressure_term = -math.log1p((prior_pressure - prior_pure) / prior_pure) / s
            # Retain every discrepancy; qualification checks the accepted criterion afterward.
            now = current[mode, salt, m]
            assert abs(state.molar_density - float(now['molar_density_mol_m3'])) <= 1e-10
            assert abs(activity.mean_ionic_molality_scale - float(now['gamma_pm_m'])) <= 1e-12
            rows.append(dict(mode=mode, salt=salt, molality_mol_kg=m, osmotic_retained=old_phi,
                             osmotic_common_pressure=scalar, osmotic_native=activity.osmotic_coefficient,
                             gamma_retained=float(old['gamma_pm_m']), gamma_current=activity.mean_ionic_molality_scale,
                             density_retained=float(old['molar_density_mol_m3']), density_current=state.molar_density,
                             osmotic_change=delta, pressure_term=pressure_term,
                             rounded_six_changed=round(scalar, 6) != round(old_phi, 6)))
        print(mode, salt, 'resolved rows checked', flush=True)
with (OUTPUT / 'resolved-comparison.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
trace = []
path = ROOT / 'data/packets/held-cameretti-sadowski-2008/3/parameters/parameters.json'
model = models['packet3', 'NaCl']
reference = model.reference(T, P, solvent_fractions=[1, 0, 0],
                            ions=('sodium-cation', 'chloride-anion'), stoichiometry=(1, 1))
pure = model.state(T, P=P, x=[1, 0, 0], phase='liquid')
e, eps, kb, na = 1.602176634e-19, 8.8541878128e-12 * 78.09, 1.380649e-23, 6.02214076e23
A = e * e / (8 * math.pi * eps * kb * T) * math.sqrt(2 * pure.molar_density * MW * na * e * e / (eps * kb * T))
for m in (1e-8, 1e-9, 1e-10, 1e-12):
    activity = reference.activity(m)
    phi = activity.osmotic_coefficient
    deb = 1 - A * math.sqrt(m) / 3
    trace.append(dict(molality_mol_kg=m, osmotic_coefficient=phi if isinstance(phi, float) else '',
                      precision='reported' if isinstance(phi, float) else str(phi),
                      debye_limit=deb, gamma_pm_m=activity.mean_ionic_molality_scale))
    assert math.isfinite(activity.mean_ionic_molality_scale)
    if m == 1e-8:
        assert isinstance(phi, float) and abs(phi - deb) <= 1e-4
    else:
        assert phi == Stop.converged_to_floor
with (OUTPUT / 'trace-verification.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=trace[0], lineterminator='\n')
    writer.writeheader()
    writer.writerows(trace)
mixed_reference = mixed.reference(T, P, solvent_fractions=[0.5, 0.5, 0, 0],
                                  ions=('sodium-cation', 'chloride-anion'), stoichiometry=(1, 1))
assert mixed_reference.activity(0.001).osmotic_coefficient is None
maximum_change = max(abs(r['osmotic_change']) for r in rows)
maximum_identity = max(abs(r['osmotic_change'] - r['pressure_term']) for r in rows)
(OUTPUT / 'verification-identity.json').write_text(json.dumps({
    'engine_provenance': list(provenance()),
    'source_base': '502e8099051c00bf169198c970f3525bbcfac7fd',
    'prior_delivery': 'c5ba43b1376777884620b3ad9aeea938dbbba116',
    'diagnosis_base': '3f5ee19f',
    'runtime_source_sha256': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in (
        'engine/native/core/numerical_policy.hpp', 'engine/native/properties/helmholtz_model.hpp',
        'engine/native/properties/infinite_dilution_reference.cpp', 'engine/src/epcsaft/_core.pyi')},
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'retained_values_path': str((PRIOR / 'results/model-curves.csv').relative_to(ROOT)),
    'retained_values_sha256': hashlib.sha256((PRIOR / 'results/model-curves.csv').read_bytes()).hexdigest(),
    'current_values_path': str((PRIOR / 'results/current-main/model-curves.csv').relative_to(ROOT)),
    'current_values_sha256': hashlib.sha256((PRIOR / 'results/current-main/model-curves.csv').read_bytes()).hexdigest(),
    'input_sha256': {mode: hashlib.sha256((ROOT / f'data/packets/held-cameretti-sadowski-2008/{mode[-1]}/parameters/parameters.json').read_bytes()).hexdigest() for mode in ('packet3', 'packet4')},
    'temperature_K': T, 'pressure_Pa': P, 'resolved_row_count': len(rows),
    'maximum_resolved_change': maximum_change, 'maximum_pressure_corrected_identity_error': maximum_identity,
    'rounded_six_changed': sum(r['rounded_six_changed'] for r in rows),
    'trace_debye_absolute_error': abs(float(trace[0]['osmotic_coefficient']) - trace[0]['debye_limit']),
    'mixed_solvent_osmotic': 'absent', 'mixed_input_path': str(mixed_input.relative_to(ROOT)),
    'mixed_input_sha256': hashlib.sha256(mixed_input.read_bytes()).hexdigest(),
    'debye_A_ln': A, 'falsifier_criteria': {'trace_absolute': 1e-4, 'resolved_absolute': 1e-6, 'pressure_corrected_absolute': 5e-9},
}, indent=2) + '\n')
assert maximum_change <= 1e-6 and maximum_identity <= 5e-9
print('maximum scalar change', maximum_change)
print('maximum pressure-corrected identity error', maximum_identity)
print('trace', trace, flush=True)
