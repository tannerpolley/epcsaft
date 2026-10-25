"""Read retained cation-mixing results into the notebook; no native calculation."""
import json

from present import D, HERE, RESULTS, read, table


def main():
    for name in ('acceptance', 'endpoint-parameters', 'endpoint-comparisons', 'yao-singular-values', 'trace-divalent'):
        rows = read(f'cation-mixing/{name}.csv')
        table(f'cation-mixing/{name}.md', list(rows[0]), [list(row.values()) for row in rows])
    out = RESULTS / 'cation-mixing'
    f, c, ref = [json.loads(p.read_text()) for p in (out/'native-fit.json', out/'checks.json', HERE/'inputs/cation-mixing-reference.json')]
    def paragraph(name, text):
        (out/f'{name}.md').write_text(text+'\n')
    paragraph('conditions', f"The opt-in solvent-mass term is $\\Delta A^r/(RT)=2\\theta n_{{Li}}n_{{Mg}}/W+\\psi n_{{Li}}n_{{Mg}}n_{{Cl}}/W^2$, with $W=M_w n_w$, θ={ref['theta_kg_per_mol']} kg/mol and ψ={ref['psi_kg2_per_mol2']} kg²/mol². At {D['conditions']['temperature_K']} K and {D['conditions']['pressure_Pa']:g} Pa, the native five-coordinate Li fit passes the owner-amended calibration limits and unchanged transfer and trace guards.")
    acceptance = read('cation-mixing/acceptance.csv')
    closest = max(acceptance, key=lambda r: abs(float(r['value']))/float(r['limit']))
    paragraph('margin', f"The largest limit ratio is {abs(float(closest['value']))/float(closest['limit']):.6f}, from {closest['check']}: {float(closest['value']):.11g} against {closest['limit']}.")
    paragraph('convergence', f"The [native fit record](results/cation-mixing/native-fit.json) reports cost {f['final_cost']:.17g} and {f['iterations']} iterations: `{f['message']}` The t07-p10 reference took {ref['iterations']} iterations: `{ref['message']}` No optimizer tolerance was tightened. The Born diameter remains at its {f['physical'][3]:g} Å lower bound; covariance is withheld, with five-coordinate Jacobian rank {f['covariance']['rank']} and projected gradient norm {f['covariance']['projected_gradient_norm']:.6g}.")
    paragraph('reproduction', f"The lead's §10/O1 ruling applies absolute {read('cation-mixing/reproduction.csv')[0]['limit']} reproduction to the prediction-error metrics against observations. Their largest difference is {c['calibration_reproduction_max_abs']:.6g}. F2/F3/F4 use their unchanged physical guards; their reported differences are optimizer-endpoint differences, with no reproduction verdict assigned. At equal parameters, observed-error metrics differ by at most {c['equal_parameter_error_metric_max_abs']:.6g}, model-comparison metrics by {c['equal_parameter_guard_max_abs']:.6g}, and Li Jacobian columns match exactly.")
    paragraph('equal-parameter', f"The [equal-parameter record](results/cation-mixing/equal-parameter.csv) separates native base values, integrated values and the analytic D7 increment: the increment error is at most {c['equal_parameter_increment_max_abs']:.6g}. The original joint-point file's last 22 raw predictions exclude D7 and are not direct targets. After adding D7, the largest raw cold/seeded observation-route difference is {c['equal_parameter_all_96_max_abs']:.6g}; this is distinct from the physical metric comparison.")
    raw = read('cation-mixing/falsifiers.csv')
    analytic = max(abs(float(r['value'])) for r in raw if r['quantity'].endswith('error'))
    finest = min(float(r['step']) for r in raw if r['quantity'].startswith('gradient-'))
    fd = max(abs(float(r['value'])) for r in raw if r['quantity'].startswith(('gradient-', 'hessian-')) and float(r['step']) == finest)
    higher = max(abs(float(r['value'])) for r in raw if r['quantity'] in ('E3','E4','E5'))
    recovery = read('cation-mixing/manufactured-recovery.csv')
    paragraph('falsifiers', f"The [native falsifiers](results/cation-mixing/falsifiers.csv) use the two declared unequal mixtures. The largest analytic/invariant discrepancy is {analytic:.6g}. Central differences show quadratic refinement until cancellation; the finest first/second-action discrepancy is {fd:.6g}. Exact E3/E4/E5 discrepancies are at most {higher:.6g}, with separate refinement records retained. Starting from zero retains nonzero coefficient sensitivities. Native manufactured recovery gives "+'; '.join(f"{r['coefficient']}={r['recovered']} for target {r['target']}" for r in recovery)+'.')
    song = json.loads((RESULTS/'trace-divalent-song.json').read_text())
    paragraph('song', f"The eight divalent-background predictions use LiCl={read('cation-mixing/trace-divalent.csv')[0]['m_LiCl_mol_kg']} mol/kg at the declared background ionic strengths. Song's published θ={song['theta_Li_Mg']} kg/mol and ψ={song['psi_Li_Mg_Cl']} kg²/mol² remain a report-only native Pitzer comparator with shared Yao lineage.")
    deriv, occupied = read('cation-mixing/neutral-derivatives.csv'), read('cation-mixing/neutral-occupied.csv')
    residual = next(float(r['A']) for r in deriv if r['quantity'] == 'residual')
    paragraph('neutral-components', f"Removing the neutral-component count admits the complete `/2` water/six-ion model with zero-amount dodecane from `gross-sadowski-2001/3` ([inputs](scripts/neutral_components_inputs.py), [native comparison](scripts/neutral_components.cpp)). At the comparison state the residual potential is $A^r/(RT)$={residual!r} in both component selections; all {len(deriv)} [comparisons](results/cation-mixing/neutral-derivatives.csv) of its value, common-species first and second amount derivatives and zero-sum composition actions differ by at most {max(float(r['absolute_difference']) for r in deriv):g}. With an occupied neutral, the [isolated contribution](results/cation-mixing/neutral-occupied.csv) is F={float(occupied[0]['F_mol']):.17g} mol against the analytic {float(occupied[0]['F_reference_mol']):g} mol (largest relative error {max(float(r['relative_error']) for r in occupied):.3g}) for every neutral amount and molar mass tried, so the term uses the declared water mass only; its molar value remains F/N.")
    paragraph('preservation', f"Absent a declaration, retained corpus values, key-solve states and work counts, and 149-row MEA F1 controls are unchanged. The full 1,467-realization #575 comparison on frozen `/1` is bitwise identical: {c['corpora']['compared_fields']['replay-575']:,} replay fields and {c['corpora']['compared_fields']['corpora']:,} corpus table fields agree, as do {c['corpora_json']['numeric_fields']:,} corpus numeric JSON fields. The only inherited amended #575 failure remains r1408; all {c['existing_packet_bytes']['files']} existing packet files are byte-identical.")



if __name__ == "__main__":
    main()
