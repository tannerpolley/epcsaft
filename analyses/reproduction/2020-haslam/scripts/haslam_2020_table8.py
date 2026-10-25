"""Declare native ePC-SAFT comparisons; retain refusals without changing packet physics."""
import csv
import hashlib
import json
import math
import time

import epcsaft
from epcsaft.parameters import parameter_document
from check_haslam_2020_source_packet import ROOT, SOURCE_HASHES, validate

REPO = ROOT.parents[2]
PROPERTIES = {'osmotic_coefficient': 'osmotic_coefficient', 'gamma_pm_m': 'mean_ionic_molality_scale'}
CATIONS = {'Li': 'lithium-cation', 'Na': 'sodium-cation', 'K': 'potassium-cation'}
ANIONS = {'Cl': 'chloride-anion', 'Br': 'bromide-anion', 'I': 'iodide-anion'}


def write_csv(name, rows):
    with (ROOT / 'results/data' / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(declare_only=False):
    start = time.perf_counter()
    source = validate()
    rows, diagnostics = [], {}
    for salt in dict.fromkeys(r['salt'] for r in source):
        reference, cause, record = None, '', {}
        if salt == 'LiI':
            cause = 'No source-backed Li+/I- pair in selected Figiel packets'
        else:
            cation = next(c for c in CATIONS if salt.startswith(c))
            components = ['water', CATIONS[cation], ANIONS[salt[len(cation):]]]
            packet = 'figiel-2025-' + ('aqueous-iodides/3' if salt.endswith('I') else 'reference-electrolytes/4')
            path = REPO / 'data/packets' / packet / 'parameters/parameters.json'
            record = {'packet': packet, 'component_ids': components,
                      'packet_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
            try:
                document = parameter_document(REPO / 'data/packets', [(packet, components)])
                record['resolved_document_sha256'] = hashlib.sha256(json.dumps(document, sort_keys=True).encode()).hexdigest()
                record['model_families'] = document['model_families']
                parameters = epcsaft.Parameters.from_mapping(document)
                if declare_only:
                    continue
                reference = epcsaft.Mixture(parameters).reference(
                    298.15, 101000.0, solvent_fractions=(1, 0, 0), ions=components[1:], stoichiometry=(1, 1))
                record.update(terminal_molality_mol_kg=reference.molality, convergence_error=reference.convergence_error,
                              log_fugacity=reference.log_fugacity, pure_solvent_log_fugacity=reference.pure_solvent_log_fugacity)
                if reference.convergence_error is not None and (not math.isfinite(reference.convergence_error) or reference.convergence_error > 5e-5):
                    cause = 'Native reference-limit diagnostic exceeds fixed tolerance'
            except Exception as exc:
                if declare_only:
                    raise
                cause = f'{type(exc).__name__}: {exc}'
        record.update(cause=cause, effective_activity_solver_options=None)
        diagnostics[salt] = record
        for point in (r for r in source if r['salt'] == salt):
            values, point_cause = {}, cause
            if not cause:
                try:
                    activity = reference.activity(float(point['molality_mol_kg']))
                    record['effective_activity_solver_options'] = activity.solver_options
                    values = {prop: getattr(activity, field) for prop, field in PROPERTIES.items()}
                except Exception as exc:
                    point_cause = f'{type(exc).__name__}: {exc}'
            for prop in PROPERTIES:
                value = values.get(prop)
                valid = isinstance(value, float) and math.isfinite(value) and value > 0
                residual = value - float(point[prop]) if valid else None
                rows.append(dict(salt=salt, observable=prop, row_id=point['row_id'],
                                 molality_mol_kg=point['molality_mol_kg'], literature_value=point[prop],
                                 model_value=value if valid else None, residual=residual,
                                 absolute_residual=abs(residual) if valid else None,
                                 relative_abs_error_percent=100 * abs(residual) / float(point[prop]) if valid else None,
                                 printed_rounding_agreement=abs(residual) <= 0.0005 if valid else None,
                                 status='evaluated' if valid else 'not_evaluated',
                                 cause='' if valid else point_cause or f'Native property unavailable: {value}',
                                 evidence='physical comparison', haslam_selection='not independently established'))
    if declare_only:
        return 0
    with (ROOT / 'inputs/source/haslam-2020-table8.csv').open() as stream:
        printed = list(csv.DictReader(stream))
    summaries = []
    for entry in (r for r in printed if r['salt'] != 'average'):
        group = [r for r in rows if (r['salt'], r['observable']) == (entry['salt'], entry['observable'])]
        valid = [r for r in group if r['status'] == 'evaluated']
        summaries.append(dict(salt=entry['salt'], observable=entry['observable'], source_points=len(group),
                              evaluated_points=len(valid), unevaluated_points=len(group) - len(valid),
                              epcsaft_aad_percent=sum(r['relative_abs_error_percent'] for r in valid) / len(group) if len(valid) == len(group) else None,
                              haslam_printed_points=entry['number_of_points'], haslam_printed_aad_percent=entry['reported_aad_percent'],
                              haslam_status='not recomputed: model absent from the Engine',
                              selection_note=entry['selection_note']))
    result = dict(engine_provenance=epcsaft.provenance(), temperature_K=298.15, pressure_Pa=101000,
                  packet_pressure_Pa=100000, pressure_extrapolation_percent=1, pressure_error_bound=None,
                  printed_rounding_absolute_tolerance=0.0005, reference_convergence_tolerance=5e-5,
                  effective_tolerances=[dict(name=t.name, value=t.value, source=t.source) for t in epcsaft.tolerances()],
                  source_hashes=SOURCE_HASHES, salts=diagnostics,
                  aad_definition='100/n * sum(abs((model-reference)/reference)); complete grid only',
                  per_point_density_residual_and_stop='not exposed by bound activity result; no second density solve',
                  haslam_model='not recomputed: model absent from the Engine',
                  run_wall_seconds=time.perf_counter() - start)
    write_csv('haslam-2020-table8.csv', rows)
    write_csv('haslam-2020-table8-aad-comparison.csv', summaries)
    (ROOT / 'results/data/haslam-2020-table8.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(f"Retained {len(rows)} property records; {sum(r['status'] == 'evaluated' for r in rows)} evaluated")
    for salt, record in diagnostics.items():
        print(salt, record['cause'] or 'reference available')
    return 0 if any(r['status'] == 'evaluated' for r in rows) else 1


def declarations():
    return main(declare_only=True)


if __name__ == '__main__':
    raise SystemExit(main())
