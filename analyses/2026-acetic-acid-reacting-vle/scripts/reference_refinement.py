"""Declare a native finite-dilution diagnostic and read native pressure re-solves."""
import csv
import hashlib
import json
import subprocess
from pathlib import Path

from run import CONFIG, ROOT, T, model, problem, solve


def main():
    identity = json.loads((ROOT / "results/executable.json").read_text())
    evidence = Path(identity["wheel"]).parent.parent
    binary = evidence / "reference_ladder"
    mix = model("primary")
    states = json.loads((ROOT / "results/states.json").read_text())
    center = next(r for r in states if r["variant"] == "primary" and not r["tighter"])
    declaration = problem(center["x_water"])._native(mix)
    doc = mix.parameters._document
    filename = evidence / "reference-declarations.txt"
    with filename.open("w") as out:
        out.write(f'{len(doc.component_ids)} {len(doc.site_component)} {len(doc.entries)}\n')
        out.write(' '.join(map(json.dumps, doc.component_ids)) + '\n')
        for comp, mult in zip(doc.site_component, doc.site_multiplicity):
            out.write(f'{comp} {mult}\n')
        out.write(' '.join(json.dumps(getattr(doc, name)) for name in ('base_choice', 'association_choice', 'electrolyte_choice', 'ionic_dispersion_choice', 'permittivity_choice', 'polar_choice')) + '\n')
        for e in doc.entries:
            out.write(f'{json.dumps(e.identity)} {json.dumps(e.field)} {len(e.index)} ' + ' '.join(map(str, e.index)) + f' {e.magnitude:.17g} {json.dumps(e.unit)}\n')
        out.write(f'{T} {CONFIG["water_molar_mass_kg_per_mol"]} {len(declaration.reactions)} {len(CONFIG["reference_molalities_mol_kg"])}\n')
        out.write(' '.join(map(str, declaration.neutral_reference.solute_molality_ratios)) + '\n')
        out.write(' '.join(map(str, CONFIG['reference_molalities_mol_kg'])) + '\n')
        for reaction, reference in zip(declaration.reactions, center['reaction_references']):
            out.write(f'{reference["source_pressure_Pa"]} {reference["contraction"]:.17g}\n')
            out.write(' '.join(map(str, reaction.stoichiometry)) + '\n')
    completed = subprocess.run([str(binary), str(filename)], check=True, capture_output=True, text=True)
    (evidence / "reference-ladder.csv").write_text(completed.stdout)
    raw = list(csv.DictReader(completed.stdout.splitlines()))
    references, shifts = [], [[], [], []]
    for reaction in range(len(declaration.reactions)):
        group = [r for r in raw if int(r['reaction_index']) == reaction]
        pairs = [r for r in group if r['terminal_pair']]
        for row, value in zip(shifts, pairs):
            row.append(float(value['log_shift']))
        references.append(dict(reaction_index=reaction, source_pressure_Pa=float(group[0]['source_pressure_Pa']),
            molalities_mol_kg=[float(r['molality_mol_kg']) for r in group],
            finite_contractions=[float(r['finite_contraction']) for r in group],
            terminal_pair_contractions=[float(r['terminal_pair']) for r in pairs],
            analytic_contraction=float(group[0]['analytic_contraction'])))
    comparisons, failures = [], []
    for number in (9, 10, 11):
        centers = [r for r in states if r['measurement_row'] == number and r['variant'] == 'primary' and not r['tighter']]
        center = centers[len(centers) // 2]
        pressures = []
        for rung, shift in enumerate(shifts):
            r = solve(center['x_water'], 'primary', reference_log_shifts=shift)
            pressures.append(r['P_Pa'])
            if not r.get('checks_met'):
                failures.append(dict(measurement_row=number, rung=rung, result=r))
        comparisons.append(dict(row=number, pressures_Pa=pressures,
            pressure_minus_analytic_Pa=[p - center['P_Pa'] for p in pressures],
            consecutive_pressure_changes_Pa=[b - a for a, b in zip(pressures, pressures[1:])],
            maximum_absolute_pressure_effect_Pa=max(abs(p - center['P_Pa']) for p in pressures)))
    output = dict(equation='reference_limit_terminal_pair', references=references, comparisons=comparisons,
        attempted=len(shifts) * len(comparisons), failures=failures,
        native_diagnostic=dict(binary=str(binary), sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
            source_sha256=hashlib.sha256((ROOT / 'scripts/reference_ladder.cpp').read_bytes()).hexdigest(),
            loaded_library_sha256={k: v for k, v in identity['loaded_sha256'].items() if Path(k).name == 'libepcsaft_native.so'},
            build_command_file=str(evidence / 'reference-build-command.json')),
        meaning='Numerical refinement comparison on a fixed neutral dilution path, not a rigorous truncation bound or an uncertainty for supplied parameters.')
    (ROOT / 'results/reference-refinement.json').write_text(json.dumps(output, indent=2) + '\n')
    assert not failures
    assert max(r['maximum_absolute_pressure_effect_Pa'] for r in comparisons) < CONFIG['numerical_pressure_criterion_Pa']
    print(json.dumps(comparisons, indent=2))


if __name__ == '__main__':
    main()
