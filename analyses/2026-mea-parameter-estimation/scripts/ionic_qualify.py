"""Build the ionic-electrostatic decision records from the runtime packet and the Born-form comparison.

Run compare.py first; this script only reads its CSV outputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from io import StringIO
from pathlib import Path

from ionic_compare import APPENDIX_F, D4, DECISIVE_DELTA_CHI_SQUARE, U_CARBAMATE, U_LOADING, U_PRESSURE, U_TEMPERATURE, U_VAPOR, packet_document
from epcsaft import Parameters

HERE = Path(__file__).resolve().parents[1] / "results/analysis/ionic-electrostatics"
ROOT = HERE.parents[4]
CAMPAIGN = ROOT / "analyses/2026-mea-parameter-estimation"
PACKET = ROOT / "data/packets/mea-co2-h2o-nine-species-estimation/5/parameters/parameters.json"
SEARCH_LOG = CAMPAIGN / "inputs/source/library/manifests/source_search_log.csv"

IONS = {
    "protonated-monoethanolamine": "MEAH+",
    "carbamate-anion": "MEACOO-",
    "bicarbonate-anion": "HCO3-",
    "carbonate-anion": "CO3^2-",
    "hydronium-cation": "H3O+",
    "hydroxide-anion": "OH-",
}
ROLES = {
    "charge_number": "formal charge used by electrostatic terms",
    "segment_diameter": "PC-SAFT segment sigma used by the molecular reference term",
    "packing_diameter": "hard-sphere packing diameter; 0.88 times segment sigma",
    "debye_huckel_diameter": "closest-approach diameter used by the Debye-Huckel term",
    "born_diameter": "radial length used by the Born charging term",
    "dispersion_energy_over_k": "PC-SAFT dispersion energy divided by Boltzmann constant",
}
BORN_BLIND = "density is Born-blind here: the runtime permittivity has no density dependence, so the Born term adds no pressure"
PERMITTIVITY_ONLY = "permittivity tests the shared bulk-permittivity input, not the Born form"
EVIDENCE_DECISIONS = {
    "SEARCH-P2-001": ("fulltext_review_pending", f"{PERMITTIVITY_ONLY}; Table 1 static values rise with loading (71 to 93 at 293 K), opposite to any ionic decrement, and come from 0.9-2.45 GHz fits"),
    "SEARCH-P2-002": ("lead_identity_resolution_pending", PERMITTIVITY_ONLY),
    "SEARCH-P2-003": ("lead_identity_resolution_pending", PERMITTIVITY_ONLY),
    "SEARCH-P2-004": ("lead_identity_resolution_pending", PERMITTIVITY_ONLY),
    "SEARCH-P3-002": ("negative_search_no_extraction", "no MEAH+ salt activity or osmotic data found; repeated 2026-09-24"),
    "SEARCH-P5-001": ("extracted_and_independently_verified", BORN_BLIND),
    "SEARCH-P5-002": ("extracted_and_independently_verified", BORN_BLIND),
    "SEARCH-P5-003": ("blocked_no_lawful_numeric_extraction", BORN_BLIND),
    "SEARCH-P5-004": ("negative_search_no_extraction", BORN_BLIND),
    "SEARCH-P5-005": ("boundary_only_no_extraction", BORN_BLIND),
    "SEARCH-P5-006": ("context_only_no_regression_target", BORN_BLIND),
}
# Searched 2026-09-24 (repository, MEA-Thermodynamics, Zotero library, open literature).
ADDED_EVIDENCE = [
    {"source": "Hilliard (2008) dissertation, Appendix D.4; Boettinger, Maiwald and Hasse (2008) doi:10.1016/j.fluid.2007.09.017",
     "observable": "water and MEA partial pressures over 30 wt% MEA at 313.15/333.15 K, loading 0.15-0.49, at NMR speciation measured at the same T and MEA strength",
     "model_comparison_eligible": True,
     "disposition": "admitted: fixed-speciation solvent fugacity; no reaction constant or CO2 pressure enters; 25 states, 50 rows; "
                    "#161 weights from Hilliard's stated uncertainties and Appendix D.4 measured temperature and MEA molality"},
    {"source": "Hilliard (2008) dissertation, Appendix F.2 (printed p. 933), 13C NMR speciation of 7 m MEA at 40 and 60 C",
     "observable": "carbamate and bicarbonate+carbonate per NMR amine total at the NMR loading, second speciation source for the same Hilliard states",
     "model_comparison_eligible": True,
     "disposition": "holdout: fitted candidates scored without refitting; internal-standard scale cancels in the ratios; "
                    "the tabulated loading breaks the carbon balance at high loading and is not used"},
    {"source": "Lee and Lee (1998) doi:10.1021/je970263h, frozen in MEA-Thermodynamics PR #122",
     "observable": "vapor pressure over aqueous MEAH+Cl-, 298.15-333.15 K",
     "model_comparison_eligible": False,
     "disposition": "#161 1(d) stopped by owner decision (2026-09-29) and superseded by MEA-Thermodynamics #121; "
                    "chloride.py fits are retained as diagnostics only, with no ranking claim"},
    {"source": "Bates and Pinching (1951) doi:10.6028/jres.046.039",
     "observable": "EMF of MEA/MEAH+Cl- buffers, 273-323 K",
     "model_comparison_eligible": False,
     "disposition": "ionic strength at most 0.08 mol/kg; the Born forms differ too little there to rank"},
    {"source": "Agieienko and Buchner (2024) doi:10.1039/d3cp05027k",
     "observable": "static permittivity of unloaded MEA-water, 298.15 K",
     "model_comparison_eligible": False,
     "disposition": f"no ions; {PERMITTIVITY_ONLY}"},
]
MEA_BORN_FORM = ("Packet v5 here states original Born (c_shell = c_dielectric = 0). The MEA-Thermodynamics application does "
                 "not: its selected record omits both coefficients and its shared_evaluation.py fills (1, 1), so the application "
                 "evaluates SSM+DS (MEA-Thermodynamics PR #125, check N1). The MEA Born-form decision belongs to "
                 "MEA-Thermodynamics #121 decision 21.")
RELEVANT_FAMILIES = {"Born_solvation_and_dielectric", "ion_interaction_parameters", "MEAH_MEACOO_segment_parameters"}
RELEVANT_EVIDENCE_SHA256 = "342b4e6474086184496dfb0661464911e59a55636c26f2105aef167b1dff1247"
FORMS = {"original-born": (0.0, 0.0), "ssm": (1.0, 0.0), "ssm-plus-ds": (1.0, 1.0)}
HILLIARD = "Hilliard (2008) dissertation, Zotero DNQUPRMT, PDF LB7LUETY"
STATED_UNCERTAINTIES = [
    {"quantity": "vapor-phase concentration (FT-IR)", "value": "+-2.0 % relative, expanded (coverage factor not stated; used as standard)",
     "used_in_ln_sigma": U_VAPOR, "locator": f"{HILLIARD}, section 2.3.4, printed p. 29 (PDF p. 86), citing Goff (2005)"},
    {"quantity": "FT-IR analyzer accuracy", "value": "+-3.0 % of full scale", "used_in_ln_sigma": None,
     "locator": f"{HILLIARD}, section 2.3.4, printed p. 28 (PDF p. 85); full scale per species not stated, superseded by the 2.0 % concentration uncertainty"},
    {"quantity": "reactor pressure", "value": "+-0.08 % of 0-8 bar full scale", "used_in_ln_sigma": round(U_PRESSURE, 5),
     "locator": f"{HILLIARD}, section 2.3.4, printed p. 27 (PDF p. 84); reactor gas at atmospheric pressure"},
    {"quantity": "reactor temperature", "value": "controlled within +-0.1 C (thermometer +-0.01 C)", "used_in_ln_sigma": f"model d ln f/dT x {U_TEMPERATURE} K",
     "locator": f"{HILLIARD}, section 2.3.4, printed pp. 27-28 (PDF pp. 84-85)"},
    {"quantity": "CO2 loading", "value": "+-2.0 % relative standard", "used_in_ln_sigma": f"model d ln f/d loading x {U_LOADING} loading",
     "locator": f"{HILLIARD}, section 2.3.4, printed p. 29 (PDF p. 86)"},
    {"quantity": "apparatus check: pure water and MEA vapor pressure", "value": "AARD 4.40 % (water, vs Kell 1984) and 7.31 % (MEA, vs DIPPR); low-temperature vapor-pressure analysis uncertainty +-2.0 %",
     "used_in_ln_sigma": None, "locator": f"{HILLIARD}, section 2.4, printed pp. 34-37 (PDF pp. 91-94)"},
    {"quantity": "Hilliard's own regression weight for P_MEA", "value": "sigma 2.0 % (DRS default)", "used_in_ln_sigma": None,
     "locator": f"{HILLIARD}, Table 13.4-2, printed p. 440 (PDF p. 497)"},
    {"quantity": "NMR amine total against the dioxane internal standard (context only)", "value": "+-3.0 % (1-5 wt% dioxane) to +-5.0 % (<1 wt%) relative; "
     "no per-species carbamate uncertainty is stated, and Hilliard's own regression uses sigma_x 1.0 % (Table 13.4-2, p. 440)", "used_in_ln_sigma": None,
     "locator": f"{HILLIARD}, section 3.11, printed pp. 62-63 (PDF pp. 119-120)"},
    {"quantity": "carbamate speciation (project assumption, not Hilliard)", "value": f"+-{U_CARBAMATE:.0%} relative, as the #163 carbamate scales",
     "used_in_ln_sigma": f"model half-spread at carbamate x(1 +- {U_CARBAMATE}); dominates MEA rows above loading 0.35", "locator": "#163 assumption"},
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(name: str) -> list[dict]:
    with (HERE / name).open() as stream:
        return list(csv.DictReader(stream))


def calibrated_ranking() -> dict:
    """#161: original Born minus SSM+DS chi-square in every calibrated case (compare.py calibrated).

    Each difference is divided by max(1, reduced chi-square of the better form), so weights that prove too small
    cannot inflate it. The ranking is decisive only if every case favours the same form by at least
    DECISIVE_DELTA_CHI_SQUARE."""
    fits = read_csv("solvent-fugacity-refit-calibrated.csv")
    scores = read_csv("solvent-fugacity-holdout.csv")

    def case(label: str, rows: list[dict], chi: str, reduced: str) -> dict:
        by = {r["model"]: r for r in rows}
        born, ds = (by[m] for m in ("original-born", "ssm-plus-ds"))
        best = min((born, ds), key=lambda r: float(r[chi]))
        delta = float(born[chi]) - float(ds[chi])
        return {"case": label, "chi_square": {m: round(float(r[chi]), 4) for m, r in by.items()},
                "reduced_chi_square": {m: round(float(r[reduced]), 4) for m, r in by.items()},
                "original_born_minus_ssm_plus_ds": round(delta, 4),
                "effective_difference": round(delta / max(1.0, float(best[reduced])), 4),
                "difference_at_observed_scatter": round(delta / float(best[reduced]), 4)}

    cases = []
    for state, bias, scale in [(s, b, c) for s, b in (("nominal", "none"), ("measured", "none"), ("measured", "fitted"))
                               for c in ("0.95", "1.0", "1.05")] + [("measured-no-carbamate-sigma", "none", "1.0")]:
        rows = [r for r in fits if (r["state"], r["loading_bias"], r["carbamate_scale"]) == (state, bias, scale)]
        cases.append({**case(f"refit {state} state, loading bias {bias}, carbamate x{scale}", rows, "chi_square", "reduced_chi_square"),
                      "decides": state == "measured",
                      "carbamate_bias_per_loading": {r["model"]: round(float(r["carbamate_bias_per_loading"]), 4) for r in rows},
                      "k_ion_water": {r["model"]: round(float(r["k_ion_water"]), 4) for r in rows},
                      "k_ion_mea": {r["model"]: round(float(r["k_ion_mea"]), 4) for r in rows},
                      "max_gradient_inf_norm": round(max(float(r["multi_start_gradient_inf_norm"]) for r in rows), 5),
                      "max_multi_start_chi_square_spread": float(f'{max(float(r["multi_start_chi_square_spread"]) for r in rows):.3g}')})
    for scale in ("0.95", "1.0", "1.05"):
        rows = [r for r in scores if r["fit_carbamate_scale"] == scale]
        cases.append({**case(f"Appendix F holdout, frozen fit from carbamate x{scale}", rows, "chi_square_frozen", "reduced_chi_square_frozen"),
                      "decides": True})
        cases.append({**case(f"Appendix F holdout, neutral offsets re-profiled, fit from carbamate x{scale}", rows,
                             "chi_square_offset_profiled", "reduced_chi_square_offset_profiled"), "decides": False})
    deciding = [c["effective_difference"] for c in cases if c["decides"]]
    same_sign = all(d > 0 for d in deciding) or all(d < 0 for d in deciding)
    return {"rule": f"decisive when every deciding case favours the same form by an effective chi-square difference >= {DECISIVE_DELTA_CHI_SQUARE}",
            "decisive": same_sign and min(abs(d) for d in deciding) >= DECISIVE_DELTA_CHI_SQUARE,
            "favoured_form_in_every_deciding_case": ("ssm-plus-ds" if deciding[0] > 0 else "original-born") if same_sign else None,
            "smallest_effective_difference": min(abs(d) for d in deciding),
            "cases": cases,
            "stated_uncertainties": STATED_UNCERTAINTIES}


def build() -> tuple[dict, str, str]:
    wheels = list((ROOT / "build/environment-wheel").glob("epcsaft-*.whl"))
    if len(wheels) != 1:
        raise SystemExit(f"expected exactly one Engine wheel, found {len(wheels)}")
    wheel = wheels[0]
    verify = subprocess.run(
        [str(ROOT / ".venv/bin/python"), str(ROOT / "tools/verify-installed-wheel.py"), "--wheel", str(wheel)],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout.strip().replace(str(ROOT), ".")
    packet = json.loads(PACKET.read_text())
    Parameters.from_mapping(packet_document())
    electrolyte = next(f for f in packet["model_families"] if f["family_id"] == "model/electrolyte")
    permittivity = next(f for f in packet["model_families"] if f["family_id"] == "model/relative_permittivity")
    runtime = {pair: name for name, pair in FORMS.items()}[(electrolyte["c_shell"], electrolyte["c_dielectric"])]

    decisions = []
    for component in packet["components"]:
        if component["component_id"] not in IONS:
            continue
        records = {r["family"]: r for r in component["coefficients"]} | {"charge_number": component["fixed"]["charge_number"]}
        for family, role in ROLES.items():
            record = records[family]
            fixed = family == "charge_number"
            decisions.append({
                "species": IONS[component["component_id"]], "component_id": component["component_id"], "coordinate": family,
                "role": role, "value": record["value"]["magnitude"], "unit": record["value"]["unit"],
                "status": "established" if fixed else "provisional", "qualification": record["qualification"],
                "candidate_domain_id": record["candidate_domain_id"], "provenance": record["provenance"],
                "decision": "retain fixed species structure" if fixed else
                "retain as provisional; not identified by the solvent-fugacity comparison; do not fit to reactive pressure",
            })

    with SEARCH_LOG.open(newline="") as stream:
        relevant = [row for row in csv.DictReader(stream) if row["parameter_family"] in RELEVANT_FAMILIES]
    relevant_sha256 = hashlib.sha256(json.dumps(relevant, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if relevant_sha256 != RELEVANT_EVIDENCE_SHA256:
        raise SystemExit("ionic evidence content changed; reassess source qualification and admission")
    by_id = {row["search_id"]: row for row in relevant}
    if set(by_id) != set(EVIDENCE_DECISIONS):
        raise SystemExit("ionic evidence inventory changed; classify every relevant search record before regeneration")
    evidence = []
    for search_id, (expected, reason) in EVIDENCE_DECISIONS.items():
        if by_id[search_id]["extraction_decision"] != expected:
            raise SystemExit(f"{search_id} extraction decision changed; reassess evidence admission")
        evidence.append({**by_id[search_id], "model_comparison_eligible": False, "exclusion_reason": reason})

    fits = [{k: (v if k == "model" else float(v)) for k, v in r.items()} for r in read_csv("solvent-fugacity-refit.csv")]
    trends = read_csv("solvent-fugacity-trends.csv")
    osmotic_values = [{"model": r["model"], "molality_mol_kg": float(r["molality_mol_kg"]), "coefficient": float(r["value"])} for r in read_csv("born-sensitivity.csv") if r["quantity"] == "osmotic_coefficient"]

    def metric(name: str) -> dict:
        own = [f for f in fits if f["model"] == name]
        water = [float(t["slope_per_loading"]) for t in trends
                 if t["model"] == name and t["component"] == "water" and t["temperature_K"] == "313.15"]
        return {
            "rows": int(own[0]["rows"]), "fitted_parameters": int(own[0]["parameters"]),
            "chi_square_by_carbamate_scale": {str(f["carbamate_scale"]): round(f["chi_square"], 2) for f in own},
            "k_ion_water_range": [round(min(f["k_ion_water"] for f in own), 3), round(max(f["k_ion_water"] for f in own), 3)],
            "k_ion_mea_range": [round(min(f["k_ion_mea"] for f in own), 3), round(max(f["k_ion_mea"] for f in own), 3)],
            "unfitted_water_slope_313K_range": [round(min(water), 3), round(max(water), 3)],
        }

    metrics = {name: metric(name) for name in FORMS}
    calibrated = calibrated_ranking()
    record = {
        "schema_version": 3,
        "issue": 117,
        "ranking_issue": 161,
        "rediagnosis_issue": 145,
        "decision": "packet_v5_ssm_rejected_on_physical_coordinates_original_born_vs_ssm_plus_ds_not_decisive_original_born_retained",
        "claim": (
            "On packet v5 at measured speciation, refits weighted by Hilliard's stated measurement uncertainties plus the project's +-5 % "
            "carbamate assumption give SSM+DS a lower chi-square than "
            "original Born in every case: three uniform carbamate scales, measured and nominal states, with and without a loading-linear "
            "carbamate bias, and on the Hilliard Appendix F NMR holdout without refitting. The difference is 0.3-1.3 chi-square units, "
            "below the decisive threshold of 10 even when the weights are shrunk to the observed scatter, so the order is consistent "
            "but not decisive. SSM stays rejected on physical grounds: implausible ion-solvent unlike energies and negative packet-coordinate "
            "aqueous osmotic coefficients at k_ij=0. No packet coordinates are changed. " + MEA_BORN_FORM
        ),
        "diagnosed_cause_of_remaining_negative": (
            "With the two ion-solvent k_ij fitted, the admissible Born forms differ in these partial pressures by less than their "
            "measurement and speciation uncertainty. At Hilliard's stated uncertainties the reduced chi-square is 0.30-0.54: water rows "
            "carry 2.2 % against about 0.8 % observed scatter, and MEA rows at high loading are dominated by speciation and loading "
            "uncertainty. The fitted k_ij absorb the unfitted Born-form slope difference. Partial pressures at matched speciation "
            "do not separate the forms while k_ij are free, also without the carbamate term in the weights; data that fix the ion-solvent "
            "interaction independently can. "
            "Density cannot help: it is Born-blind in this model."
        ),
        "needed_to_separate_original_born_from_ssm_plus_ds": [
            "MEA-Thermodynamics #121 decision 21: both forms fitted on the full MEA calibration data with the exact-derivative fitter, 40-80 C pass rule.",
            "MEAH+ salt osmotic or activity data that fix the ion-solvent k_ij independently of the Born form.",
            "Solvent partial pressures and NMR speciation measured on the same samples, with reported uncertainty.",
        ],
        "forbidden_comparison": "reactive pressure alone",
        "calibrated_ranking": calibrated,
        "osmotic_coefficient_check": {
            "basis": "aqueous 1:1 MEAH+/MEACOO- at 298.15 K with packet ion-solvent k_ij=0",
            "values": osmotic_values,
            "source": "born-sensitivity.csv",
        },
        "independent_eligible_observation_count": sum(1 for e in ADDED_EVIDENCE if e["model_comparison_eligible"]),
        "runtime_model_configuration": {
            "parameter_document": str(PACKET.relative_to(ROOT)),
            "born_form": runtime,
            "c_shell": electrolyte["c_shell"],
            "c_dielectric": electrolyte["c_dielectric"],
            "bulk_permittivity": permittivity["choice"],
            "matches_selection": runtime == "original-born",
        },
        "model_decisions": [
            {"model": "original-born", "configuration": {"c_shell": 0.0, "c_dielectric": 0.0}, "status": "candidate",
             "metric": metrics["original-born"],
             "disposition": "packet v5 baseline in this repository; 0.3-1.3 chi-square above SSM+DS in every calibrated case, not decisively worse"},
            {"model": "ssm", "configuration": {"c_shell": 1.0, "c_dielectric": 0.0}, "status": "rejected",
             "metric": metrics["ssm"],
             "disposition": (
                 "rejected on physical grounds despite the lowest chi-square: fitted ion-solvent k_ij values are implausible, "
                 "and packet-coordinate aqueous osmotic coefficients at k_ij=0 are negative at 1 mol/kg and above; see the retained diagnostics"
             )},
            {"model": "ssm-plus-ds", "configuration": {"c_shell": 1.0, "c_dielectric": 1.0}, "status": "candidate",
             "metric": metrics["ssm-plus-ds"],
             "disposition": "admissible alternative; lower chi-square than original Born in every calibrated case by 0.3-1.3, not decisive"},
        ],
        "fitted_candidate_coordinates": {
            "scope": "one shared k_ij for MEAH+ and MEACOO- with water and with MEA, per Born form; candidates only, not promoted",
            "fits": fits,
        },
        "coordinate_decisions": decisions,
        "evidence_inventory": evidence,
        "added_evidence": ADDED_EVIDENCE,
        "evidence_admission": {
            "solvent_fugacity_at_measured_speciation": "admitted",
            "direct_meah_salt_activity_or_osmotic": "diagnostic_only (Lee and Lee 1998; superseded by MEA-Thermodynamics #121)",
            "volumetric": "cannot discriminate: Born-blind",
            "static_permittivity": "cannot discriminate: shared bulk-permittivity input",
            "reactive_pressure": "forbidden_for_ionic_identification",
        },
        "domain": "30 wt% MEA, 313.15 and 333.15 K, 6.6-22 kPa, CO2 loading 0.15-0.49 mol/mol",
        "source_authority": {
            "parameter_document_sha256": sha256(PACKET),
            "search_log": str(SEARCH_LOG.relative_to(ROOT)),
            "search_log_sha256": sha256(SEARCH_LOG),
            "relevant_evidence_sha256": relevant_sha256,
            "comparison_outputs_sha256": {n: sha256(HERE / n) for n in (
                "born-sensitivity.csv", "solvent-fugacity.csv", "solvent-fugacity-trends.csv", "solvent-fugacity-refit.csv",
                "solvent-fugacity-refit-v4.csv", "coordinate-sensitivity.csv", "solvent-fugacity-weights.csv",
                "solvent-fugacity-refit-calibrated.csv", "solvent-fugacity-holdout.csv")},
            "hilliard_transcriptions_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in (D4, APPENDIX_F)},
            "born_model_semantics": "engine/docs/science/full-eos.md section 7",
        },
        "installed_wheel_check": {"wheel": str(wheel.relative_to(ROOT)), "wheel_sha256": sha256(wheel), "verification": verify},
        "falsification_conditions": [
            "A holdout at matched speciation separates original Born from SSM+DS by an effective chi-square difference of at least 10.",
            "MEAH+ salt activity or osmotic data with a qualified counterion rank the forms.",
            "A refit of the ionic coordinates lets SSM reproduce the flat water partial pressure with physical k_ij.",
        ],
    }

    table = StringIO()
    fields = ["species", "component_id", "coordinate", "value", "unit", "status", "qualification",
              "candidate_domain_id", "source_id", "source_domain_id", "source_locator", "decision"]
    writer = csv.DictWriter(table, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in decisions:
        writer.writerow({**{k: row[k] for k in fields if k in row}, "source_id": row["provenance"]["source_id"],
                         "source_domain_id": row["provenance"]["domain_id"], "source_locator": row["provenance"]["locator"]})

    colors = {"established": "#2f855a", "provisional": "#d69e2e"}
    labels = ["charge", "segment σ", "packing d", "DH d", "Born r", "ε/k"]
    cell_w, cell_h, left, top = 120, 42, 120, 70
    svg = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="410" viewBox="0 0 900 410">',
        '<rect width="900" height="410" fill="white"/>',
        '<style>text{font-family:system-ui,sans-serif;fill:#1a202c}.small{font-size:12px}.label{font-size:14px}.title{font-size:20px;font-weight:600}</style>',
        '<text x="20" y="30" class="title">Ionic-coordinate status</text>',
        '<text x="20" y="50" class="small">Green: established species structure. Amber: inherited value, not identified by independent MEA data.</text>',
    ]
    for column, label in enumerate(labels):
        svg.append(f'<text x="{left + column * cell_w + cell_w / 2}" y="{top - 10}" text-anchor="middle" class="small">{label}</text>')
    for row_index, species in enumerate(IONS.values()):
        y = top + row_index * cell_h
        svg.append(f'<text x="{left - 12}" y="{y + 27}" text-anchor="end" class="label">{species}</text>')
        for column, family in enumerate(ROLES):
            status = next(d["status"] for d in decisions if d["species"] == species and d["coordinate"] == family)
            x = left + column * cell_w
            svg.append(f'<rect x="{x + 2}" y="{y + 2}" width="{cell_w - 4}" height="{cell_h - 4}" rx="4" fill="{colors[status]}"/>')
            svg.append(f'<text x="{x + cell_w / 2}" y="{y + 27}" text-anchor="middle" class="small" style="fill:white">{status}</text>')
    chi = next(c for c in calibrated["cases"] if c["case"] == "refit measured state, loading bias none, carbamate x1.0")["chi_square"]
    svg.extend([
        f'<text x="20" y="350" class="label">Hilliard-weighted refit (50 rows, carbamate x1.0), chi-square: original Born {chi["original-born"]:.1f}, '
        f'SSM {chi["ssm"]:.1f}, SSM+DS {chi["ssm-plus-ds"]:.1f}</text>',
        '<text x="20" y="372" class="small">SSM rejected on physical grounds. SSM+DS lower in every case by 0.3-1.3: not decisive. Packet v5: original Born; MEA application: SSM+DS (MEA PR #125).</text>',
        '</svg>',
    ])
    return record, table.getvalue(), "\n".join(svg) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    record, table, svg = build()
    rendered = {"decision-table.json": json.dumps(record, indent=2, sort_keys=True) + "\n",
                "decision-table.csv": table, "status-summary.svg": svg}
    if args.check:
        stale = [name for name, content in rendered.items() if (HERE / name).read_text() != content]
        if stale:
            raise SystemExit(f"stale generated evidence: {', '.join(stale)}")
        print("ionic qualification evidence is current")
        return
    for name, content in rendered.items():
        (HERE / name).write_text(content)


def declarations():
    return Parameters.from_mapping(packet_document())


if __name__ == "__main__":
    main()
