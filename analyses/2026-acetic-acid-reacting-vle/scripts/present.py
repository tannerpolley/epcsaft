"""Prepare notebook displays from retained tables; never import the Engine."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def read(name):
    return list(csv.DictReader((ROOT / "results" / name).open()))


def table(name, headers, rows):
    text = "| " + " | ".join(headers) + " |\n| " + " | ".join("---" for _ in headers) + " |\n"
    text += "".join("| " + " | ".join(map(str, row)) + " |\n" for row in rows)
    (ROOT / "figures" / name).write_text(text)


def main():
    out = ROOT / "figures"
    out.mkdir(exist_ok=True)
    s = json.loads((ROOT / "results/summary.json").read_text())
    c = json.loads((ROOT / "inputs/reactions.json").read_text())
    p, control, interval = read("pressures.csv"), read("controls.csv"), read("uncertainty.csv")
    def f(row, key):
        return float(row[key])
    table("pressures.md", ["ThermoML row", "Analytical water fraction", "Measured P ± expanded U / kPa", "Calculated P / kPa", "Signed error / kPa", "Signed relative error / %"],
          [[r["row"], f'{f(r,"x_water"):.4f}', f'{f(r,"P_measured_kPa"):.3f} ± {f(r,"U95_P_kPa"):.3f}', f'{f(r,"P_calculated_kPa"):.6f}', f'{f(r,"signed_error_kPa"):+.6f}', f'{f(r,"signed_relative_error_percent"):+.5f}'] for r in p])
    table("controls.md", ["Row", "No dissociation P / kPa", "Dissociation effect / Pa", "Original water-only dielectric effect of dissociation / Pa", "Acid-permittivity effect / Pa", "Water-ionization omission, max over input interval / Pa"],
          [[r["row"], f'{f(r,"P_no_dissociation_kPa"):.6f}', f'{f(r,"dissociation_effect_Pa"):+.6f}', f'{f(r,"original_dissociation_effect_Pa"):+.6f}', f'{f(r,"acid_permittivity_effect_Pa"):+.6f}', f'{f(r,"autoprotolysis_grid_max_abs_Pa"):.2e}'] for r in control])
    table("uncertainty.md", ["Row", "Input water-fraction interval", "Calculated P interval / kPa", "Measured P interval / kPa", "Interval separation / kPa", "Disjoint?"],
          [[r["row"], f'[{f(r,"x_low"):.4f}, {f(r,"x_high"):.4f}]', f'[{f(r,"P_calculated_low_kPa"):.6f}, {f(r,"P_calculated_high_kPa"):.6f}]', f'[{f(r,"P_measured_low_kPa"):.3f}, {f(r,"P_measured_high_kPa"):.3f}]', f'{f(r,"separation_kPa"):.6f}', r["disjoint"]] for r in interval])
    table("numerical.md", ["Row", "Tighter-solve change / Pa", "Maximum reference-refinement change / Pa", "Native local numerical pressure estimate, max over interval / Pa", "Exact dP/dx, sampled range / Pa"],
          [[r["row"], f'{f(r,"tighter_solve_abs_effect_Pa"):.2e}', f'{f(r,"reference_refinement_max_abs_Pa"):.2e}', f'{f(r,"native_numerical_pressure_estimate_Pa"):.2e}', f'[{f(q,"min_dP_dx_Pa"):.2f}, {f(q,"max_dP_dx_Pa"):.2f}]'] for r, q in zip(control, interval)])
    (out / "outcome.md").write_text(
        f'At {s["temperature_K"]:.2f} K, the fixed published-parameter calculation overpredicts all {s["measurements"]} measured total pressures. '
        f'The mean signed error is **{s["mean_signed_error_kPa"]:+.6f} kPa** and the root-mean-square error is **{s["RMSE_kPa"]:.6f} kPa**. '
        f'The mean absolute relative error is **{s["mean_absolute_relative_error_percent"]:.5f}%**, with a maximum of **{s["maximum_absolute_relative_error_percent"]:.5f}%**. '
        f'All {s["interval_disjoint_count"]} conditional calculated intervals are disjoint from the measurement intervals.\n')
    (out / "conditions.md").write_text(
        f'The analytical feed is one mole of water plus acetic acid at **{c["temperature_K"]:.2f} K**. '
        f'The acid constant is **$K_{{a,m}}={c["Ka"]:.3e}$**, with a standard molality of **{c["standard_molality_mol_per_kg"]:g} mol kg⁻¹** and adopted source pressure **{c["acid_reference_pressure_Pa"]:g} Pa**. '
        f'The water-ionization constant is **$pK_w={c["pKw"]:.2f}$** at source pressure **{c["water_reference_pressure_Pa"]:g} Pa**. '
        f'The adopted acetic-acid relative permittivity is **{c["acid_permittivity"]["value"]:.9f}**.\n')
    (out / "numerical-summary.md").write_text(
        f'All **{s["attempted"]} native bubble solves** passed the requested numerical checks, including the separate reference comparisons; **{s["numerical_failures"]} numerical failures** are retained in `results/failures.json`. '
        f'The maximum row residual is {s["max_row_residual"]:.2e}; maximum absolute charge and inventory errors are {s["max_charge_error_mol"]:.2e} and {s["max_inventory_error_mol"]:.2e} mol. '
        f'The maximum relative pressure-closure error is {s["max_relative_pressure_closure"]:.2e}; maximum neutral log-fugacity mismatch and liquid log-reaction residual are {s["max_neutral_log_fugacity_difference"]:.2e} and {s["max_reaction_log_residual"]:.2e}.\n\n'
        f'The ordinary row tolerance is {c["row_tolerance"]:.0e}, repeated at {c["tighter_row_tolerance"]:.0e}. '
        f'The checks require charge error ≤ {c["charge_tolerance_mol"]:.0e} mol, inventory error ≤ {c["inventory_tolerance_mol"]:.0e} mol, log-fugacity and log-reaction residuals ≤ {c["log_fugacity_tolerance"]:.0e}, and relative pressure closure ≤ {c["relative_pressure_tolerance"]:.0e}. '
        f'The pressure-refinement criterion and conservative interval expansion are **{c["numerical_pressure_criterion_Pa"]:g} Pa**. '
        f'Each composition interval retains **{c["interval_grid_points"]} native solves** with exact first derivatives.\n')
    reference = json.loads((ROOT / "results/reference-refinement.json").read_text())
    (out / "reference-conditions.md").write_text(
        'The independent neutral dilution ladder uses molalities **' + ', '.join(f'{m:.0e}' for m in reference["references"][0]["molalities_mol_kg"]) +
        ' mol kg⁻¹**. Successive square-root terminal-pair extrapolations perturb both reaction references in native bubble re-solves. '
        'The resulting pressure differences quantify this refinement check; they are not a rigorous truncation enclosure.\n')
    original_signals = [abs(f(r, "original_dissociation_effect_Pa")) for r in control]
    uncertainties = [1000 * f(r, "U95_P_kPa") for r in p]
    (out / "speciation-limit.md").write_text(
        f'With the original dielectric declaration, the pressure signal from dissociation is only about {min(original_signals):.0f}–{max(original_signals):.0f} Pa, '
        f'against {min(uncertainties):.0f}–{max(uncertainties):.0f} Pa measurement uncertainty.\n')
    x = [f(r, "x_water") for r in p]
    fig, ax = plt.subplots(figsize=(7.4, 4.6), layout="constrained")
    ax.errorbar(x, [f(r, "P_measured_kPa") for r in p],
        xerr=[(f(q, "x_high") - f(q, "x_low")) / 2 for q in interval],
        yerr=[f(r, "U95_P_kPa") for r in p], fmt="o", color="#0072B2", capsize=4,
        label="Measured, compiler expanded uncertainty")
    center = [f(r, "P_calculated_kPa") for r in p]
    ax.errorbar(x, center,
        yerr=[[a - f(q, "P_calculated_low_kPa") for a, q in zip(center, interval)],
              [f(q, "P_calculated_high_kPa") - a for a, q in zip(center, interval)]],
        fmt="s", color="#D55E00", capsize=4, label="Calculated, composition-interval effect")
    ax.set(xlabel="Analytical liquid water mole fraction", ylabel="Total equilibrium pressure / kPa",
           title="Calculated pressures exceed the three measurement intervals")
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax.grid(alpha=.18)
    fig.savefig(out / "pressure-comparison.png", dpi=180)
    plt.close(fig)
    from decompose import main as decompose
    decompose()
    print("Prepared notebook tables, retained-value prose and pressure figure.")


if __name__ == "__main__":
    main()
