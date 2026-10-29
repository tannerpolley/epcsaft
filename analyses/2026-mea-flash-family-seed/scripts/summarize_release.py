"""Check retained v1.0.0 states and compare them with main-8298c602; no Engine calls."""
import ast
import csv
import hashlib
import json
import math
from pathlib import Path

import sympy

from compare import difference

ROOT = Path(__file__).resolve().parents[3]
RESULTS = Path(__file__).resolve().parents[1] / 'results'
OUT = RESULTS / 'v1.0.0'
TOLERANCE = 1e-10  # Existing declared row and frozen state-comparison criteria.


def main():
    current = [json.loads(s) for s in (OUT / 'family.jsonl').read_text().splitlines()]
    previous = [json.loads(s) for s in (RESULTS / 'main-8298c602/family.jsonl').read_text().splitlines()]
    points = [tuple(map(float, s.split())) for s in
              (ROOT / 'analyses/2026-mea-flash-family-seed/scripts/family-points.txt').read_text().splitlines() if s.strip()]
    provenance = json.loads((OUT / 'provenance.json').read_text())
    assert hashlib.sha256((OUT / 'family.jsonl').read_bytes()).hexdigest() == provenance['family_sha256']
    for name, expected in provenance['source_sha256'].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    assert len(current) == len(previous) == len(points) == 73
    tree = ast.parse((ROOT / 'engine/tests/equilibrium/test_trace_reliability.py').read_text())
    constants = {n.targets[0].id: ast.literal_eval(n.value) for n in tree.body
                 if isinstance(n, ast.Assign) and isinstance(n.targets[0], ast.Name)
                 and n.targets[0].id in ('_MEA', '_MEA_REACTIONS')}
    species = constants['_MEA']
    reactions = sympy.Matrix([[sympy.Rational(r.get(s, 0)) for s in species]
                             for r in constants['_MEA_REACTIONS']])
    invariants = [[float(v) for v in column] for column in reactions.nullspace()]
    comparison, checks = [], []
    def accepted(row):
        return row['success'] and row.get('met', False)
    for index, (new, old, point) in enumerate(zip(current, previous, points, strict=True)):
        assert new['index'] == old['index'] == index
        assert (new['loading'], new['T'], new['P']) == (old['loading'], old['T'], old['P']) == point
        assert new['engine_provenance'] == provenance['engine_provenance']
        both = accepted(new) and accepted(old)
        delta, where = difference(new, old) if both else (None, None)
        comparison.append(dict(index=index, loading=new['loading'], T_K=new['T'], P_Pa=new['P'],
                               previous_converged=accepted(old), release_converged=accepted(new),
                               state_difference=delta, worst_quantity=where,
                               previous_message=old['message'], release_message=new['message']))
        if not accepted(new):
            continue
        n, x, rho = (new[k] for k in ('amounts', 'mole_fractions', 'molar_densities'))
        assert len(n) == len(rho) == 2 and len(x) == 2 * len(species)
        assert all(math.isfinite(v) for v in n + x + rho + [new['residual']])
        assert min(n + rho) > 0 and min(x) >= 0 and max(x) <= 1
        normalization = max(abs(math.fsum(x[j:j + len(species)]) - 1) for j in (0, len(species)))
        total = [math.fsum(n[j] * x[j * len(species) + k] for j in range(2)) for k in range(len(species))]
        feed = [new['loading'] * 9.7292821, 9.7292821, 76.9751925] + [0.] * (len(species) - 3)
        # Each nullspace vector is an exact reaction-invariant linear combination of species amounts.
        errors = [abs(math.fsum(v * (a - b) for v, a, b in zip(vector, total, feed)))
                  for vector in invariants]
        scaled = [error / max(1., math.fsum(abs(v * b) for v, b in zip(vector, feed)))
                  for vector, error in zip(invariants, errors)]
        assert new['residual'] <= TOLERANCE and normalization <= TOLERANCE and max(scaled) <= TOLERANCE
        checks.append(dict(index=index, residual=new['residual'], phase_fraction_sum_error=normalization,
                           max_reaction_invariant_error_mol=max(errors), max_scaled_invariant_error=max(scaled)))
    common = [r for r in comparison if r['state_difference'] is not None]
    worst = max(common, key=lambda r: r['state_difference'])
    assert worst['state_difference'] <= TOLERANCE
    summary = dict(total=len(current), converged=sum(map(accepted, current)),
                   previous_converged=sum(map(accepted, previous)), common=len(common),
                   lost=[r['index'] for r in comparison if r['previous_converged'] and not r['release_converged']],
                   recovered=[r['index'] for r in comparison if not r['previous_converged'] and r['release_converged']],
                   failed=[r for r in comparison if not r['release_converged']],
                   row_tolerance=TOLERANCE, state_criterion=TOLERANCE,
                   state_metric='max |a-b|/max(|a|,|b|,1e-14), over phase amounts, mole fractions and molar densities',
                   max_state_difference=worst['state_difference'], worst_state_index=worst['index'],
                   worst_quantity=worst['worst_quantity'],
                   largest_converged_residual=max(r['residual'] for r in checks),
                   largest_phase_fraction_sum_error=max(r['phase_fraction_sum_error'] for r in checks),
                   reaction_invariant_vectors=invariants,
                   largest_reaction_invariant_error_mol=max(r['max_reaction_invariant_error_mol'] for r in checks),
                   largest_scaled_invariant_error=max(r['max_scaled_invariant_error'] for r in checks),
                   invariant_scale='max(1 mol, sum_i |invariant_i * feed_i|); exact nullspace of declared reaction rows',
                   physical_array_checks='finite values; positive phase amounts/densities; fractions in [0,1]; normalized phases',
                   key_flash_work={k: current[0][k] for k in ('iterations', 'evaluations', 'hessian_calls')})
    for name, rows in (('comparison.csv', comparison), ('state-checks.csv', checks)):
        with (OUT / name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
