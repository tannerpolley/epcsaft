"""Validate the source-only Ascani 2022 Case Study 2 tracer contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CSV = ROOT / "inputs/source" / "ascani-2022-case-study-2-tracer.csv"
DEFAULT_METADATA = ROOT / "inputs/source" / "ascani-2022-case-study-2-tracer.yaml"

EXPECTED_COLUMNS = [
    "case_id",
    "role",
    "temperature_K",
    "pressure_Pa",
    "feed_w_water_formula",
    "feed_w_butanol_formula",
    "feed_w_nacl_formula",
    "feed_w_kcl_formula",
    "feed_x_water_explicit",
    "feed_x_butanol_explicit",
    "feed_x_na_explicit",
    "feed_x_k_explicit",
    "feed_x_cl_explicit",
    "organic_x_water_formula",
    "organic_x_butanol_formula",
    "organic_x_nacl_formula",
    "organic_x_kcl_formula",
    "aqueous_x_water_formula",
    "aqueous_x_butanol_formula",
    "aqueous_x_nacl_formula",
    "aqueous_x_kcl_formula",
    "organic_x_water_explicit",
    "organic_x_butanol_explicit",
    "organic_x_na_explicit",
    "organic_x_k_explicit",
    "organic_x_cl_explicit",
    "aqueous_x_water_explicit",
    "aqueous_x_butanol_explicit",
    "aqueous_x_na_explicit",
    "aqueous_x_k_explicit",
    "aqueous_x_cl_explicit",
    "ln_f_water_bar",
    "ln_f_butanol_bar",
    "ln_f_kcl_pair_bar",
    "ln_f_nacl_pair_bar",
]
EXPECTED_SOURCE_HASHES = {
    "main_markdown": "c0b73c10aa1ce9830e29f34aa3c1d1af4b889971959c3245a2deb7efdd979cd6",
    "supporting_information_markdown": "a6a61508cbaae805f2e360686785318953747a66c27d9102936e77be7f472c03",
}
EXPECTED_AUTHORITY_PDF_HASHES = {
    "main_markdown": "3c446bfa5274f8c4f0963f23c76a5cf7a87496db94dc29896369f4946b53e6e6",
    "supporting_information_markdown": "33f43e99557130757d0893d78772cbd59c05da69ff2e85377de242c04b8af4d0",
}
EXPECTED_MIGRATION_BINDING = {
    "decision": "D-024",
    "gate_commit": "4527864ffcd37f5e9a524500dfd99d5a34c85672",
    "gate_tree": "6a4f9711787c333c146566be8947df2c60f1fc68",
}
EXPECTED_D026_BINDING = {
    "decision": "D-026",
    "gate_commit": "3a4ef0a0c6b98c43405d3cafc1ac4f5f87afa68d",
    "gate_tree": "9307c3f79581b6e0479d4ac2468932b2a68e5f5b",
    "selection": "ascani_case_study_2_source_complete_pending_distinct_provider_bundle",
}
EXPECTED_PARAMETER_SOURCE_HASHES = {
    "ascani_held_2021": {
        "pdf_sha256": "9ab259a8dfb27a052fcf49782e6ab75132140d94c5f8f695215a2703f1d010ab",
        "permanent_lab_pdf_sha256": "b26bb94b84e676c15a5ff40f878f0e0d1488c6994990c3c3e7d968630f2fa464",
        "markdown_sha256": "948a50f7d447a5d3622b3bd401932d4dc84a3b715ae93dd0ffbb01fe64866c68",
    },
    "nann_held_sadowski_2013": {
        "pdf_sha256": "aeb3562f76ac9b3a8ee11779696933aead55a5bd0e6e6bb67bd82f38bf691313",
        "mathpix_sha256": "de0de1510c0520d4c900421ddab3edf84e613c798840617f5d07095e427d4426",
    },
    "held_et_al_2014": {
        "pdf_sha256": "dea9aa05e2ee8eb1c675873fa9d8737943312484874d844f01b45613da261acf",
        "markdown_sha256": "b8b1e46bf870224de5de68b5989f9cb377d17445d87109a5462a94f1efaafbda",
    },
}
DERIVATION_ABS = Decimal("1e-27")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decimals(row: dict[str, str], names: list[str]) -> list[Decimal]:
    return [Decimal(row[name]) for name in names]


def explicit_phase(formula: list[Decimal]) -> list[Decimal]:
    water, butanol, nacl, kcl = formula
    amounts = [water, butanol, nacl, kcl, nacl + kcl]
    total = sum(amounts)
    return [amount / total for amount in amounts]


def check(
    csv_path: Path, metadata_path: Path, verify_local_sources: bool = False
) -> dict[str, object]:
    with csv_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise ValueError(f"unexpected CSV columns: {reader.fieldnames}")
        rows = list(reader)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    if len(rows) != 1:
        raise ValueError(f"expected one frozen tracer row, found {len(rows)}")
    row = rows[0]
    if metadata["migration_binding"] != EXPECTED_MIGRATION_BINDING:
        raise ValueError("unexpected D-024 Migration binding")
    if metadata["d026_selection_binding"] != EXPECTED_D026_BINDING:
        raise ValueError("unexpected D-026 Migration binding")
    if metadata["citation"]["doi"] != "10.1021/acs.jced.1c00866":
        raise ValueError("unexpected Ascani DOI")
    source_hashes = {
        name: source["sha256"] for name, source in metadata["primary_sources"].items()
    }
    if source_hashes != EXPECTED_SOURCE_HASHES:
        raise ValueError("primary Markdown source hashes changed")
    authority_pdf_hashes = {
        name: source["authority_pdf_sha256"]
        for name, source in metadata["primary_sources"].items()
    }
    if authority_pdf_hashes != EXPECTED_AUTHORITY_PDF_HASHES:
        raise ValueError("primary PDF source hashes changed")
    parameter_source_hashes = {
        name: {
            key: source[key]
            for key in (
                "pdf_sha256",
                "permanent_lab_pdf_sha256",
                "markdown_sha256",
                "mathpix_sha256",
            )
            if key in source
        }
        for name, source in metadata["bound_parameter_sources"].items()
    }
    if parameter_source_hashes != EXPECTED_PARAMETER_SOURCE_HASHES:
        raise ValueError("bound parameter source hashes changed")
    csv_hash = sha256(csv_path)
    if metadata["retained_files"]["tracer_csv"]["sha256"] != csv_hash:
        raise ValueError("retained tracer CSV hash does not match metadata")

    if (
        row["case_id"] != "ascani2022-case-study-2"
        or row["role"] != "secondary_external_held2_tracer"
    ):
        raise ValueError("unexpected tracer identity")
    if Decimal(row["temperature_K"]) != Decimal("298.15") or Decimal(
        row["pressure_Pa"]
    ) != Decimal("100000"):
        raise ValueError("unexpected tracer state")

    feed_mass = decimals(
        row,
        [
            "feed_w_water_formula",
            "feed_w_butanol_formula",
            "feed_w_nacl_formula",
            "feed_w_kcl_formula",
        ],
    )
    if feed_mass != [
        Decimal(value) for value in ("0.8094", "0.1728", "0.0054", "0.0124")
    ] or sum(feed_mass) != Decimal(1):
        raise ValueError("formula-basis feed changed")
    explicit_feed = decimals(
        row,
        [
            "feed_x_water_explicit",
            "feed_x_butanol_explicit",
            "feed_x_na_explicit",
            "feed_x_k_explicit",
            "feed_x_cl_explicit",
        ],
    )
    if explicit_feed != [
        Decimal("0.9403742328474496"),
        Decimal("0.04879524350242891"),
        Decimal("0.0019339333595015447"),
        Decimal("0.0034813284655591703"),
        Decimal("0.005415261825060715"),
    ]:
        raise ValueError("explicit species feed changed")
    if (
        abs(sum(explicit_feed) - Decimal(1)) > Decimal("1e-15")
        or explicit_feed[2] + explicit_feed[3] != explicit_feed[4]
    ):
        raise ValueError(
            "explicit feed does not preserve normalization and charge stoichiometry"
        )

    formula_names = ["water", "butanol", "nacl", "kcl"]
    organic_formula = decimals(
        row, [f"organic_x_{name}_formula" for name in formula_names]
    )
    aqueous_formula = decimals(
        row, [f"aqueous_x_{name}_formula" for name in formula_names]
    )
    if organic_formula != [
        Decimal(value) for value in ("0.4426", "0.5570", "0.0000415", "0.000420")
    ]:
        raise ValueError("paper organic formula composition changed")
    if aqueous_formula != [
        Decimal(value) for value in ("0.9627", "0.0122", "0.0076", "0.0174")
    ]:
        raise ValueError("paper aqueous formula composition changed")
    explicit_names = ["water", "butanol", "na", "k", "cl"]
    for phase, formula in (("organic", organic_formula), ("aqueous", aqueous_formula)):
        observed = decimals(
            row, [f"{phase}_x_{name}_explicit" for name in explicit_names]
        )
        expected = explicit_phase(formula)
        if any(
            abs(left - right) > DERIVATION_ABS
            for left, right in zip(observed, expected, strict=True)
        ):
            raise ValueError(f"{phase} formula-to-explicit transformation changed")
        if (
            abs(sum(observed) - Decimal(1)) > DERIVATION_ABS
            or abs(observed[2] + observed[3] - observed[4]) > DERIVATION_ABS
        ):
            raise ValueError(
                f"{phase} explicit composition is not normalized and charge balanced"
            )

    fugacity = decimals(
        row,
        [
            "ln_f_water_bar",
            "ln_f_butanol_bar",
            "ln_f_kcl_pair_bar",
            "ln_f_nacl_pair_bar",
        ],
    )
    if fugacity != [
        Decimal(value) for value in ("-3.521", "-5.088", "-206.733", "-224.891")
    ]:
        raise ValueError("Table 5 log-fugacity values changed")

    if metadata["source_classification"]["model_comparison_allowance"] is not None:
        raise ValueError("source packet must not invent a model cutoff")
    if len(metadata["missing_upstream_provenance"]) != 3:
        raise ValueError("missing-upstream provenance list changed")
    archive = metadata["lab_archive_provenance"]
    negative_space = " ".join(archive["negative_space"])
    for token in (
        "Do not reuse",
        "seeds",
        "controller",
        "tolerances",
        "cannot validate",
    ):
        if token not in negative_space:
            raise ValueError(f"archive negative-space boundary missing: {token}")
    if (
        archive["historical_witnesses"]["strong_negative_witness_min_tpd"]
        != "-0.09607343786579076"
    ):
        raise ValueError("historical negative witness changed")

    source_screen = metadata["d026_source_completeness_screen"]
    if (
        source_screen["case_source_status"]
        != "SOURCE_COMPLETE_FOR_BOUNDED_ASCANI_EPCSAFT_ADVANCED_MODEL"
        or source_screen["unresolved_ascani_source_records"]
        or source_screen["provider_status"] != "PROVIDER_NOT_YET_CAPABLE"
        or source_screen["residual_source_need"]
        != "NONE_FOR_BOUNDED_298_15_K_ASCANI_MODEL"
    ):
        raise ValueError("D-026 Ascani source screen changed")
    if source_screen["hybrid_figiel_status"] != (
        "NOT_SOURCE_COMPLETE_FOR_HYBRID_FIGIEL_SSM_DS_1_BUTANOL"
    ):
        raise ValueError("Figiel hybrid source boundary changed")
    if (
        "Ascani 2022 SI Eq. S4"
        not in source_screen["water_butanol_l_ij_sign_convention"]["target_authority"]
    ):
        raise ValueError("water/1-butanol l_ij target convention changed")
    if source_screen["model_output"] != "not_run":
        raise ValueError("D-026 source screen must remain model-free")

    support = metadata["supporting_numerical_literature"]
    if (
        support["sha256"]
        != "501d9bdb5dfd89cf5584e050cfcd9cc6580fac5395f809cb91d16d81d1bd1f38"
    ):
        raise ValueError("Belov-Aristova source hash changed")
    if support["line_ranges"] != {
        "broad_concentration_scales_conditioning": "23-27",
        "exact_derivatives_and_automatic_differentiation": "121-132",
        "dimensionless_objective": "134-138",
        "nonnegative_variables_phase_sums_material_balances": "141-187",
        "inventory_normalization_and_review_checks": "301-307",
    }:
        raise ValueError("Belov-Aristova line locators changed")
    if (
        support["classification"]
        != "general numerical-conditioning support only; not electrolyte-LLE or HELD2 evidence"
    ):
        raise ValueError("Belov-Aristova classification changed")

    challenge = metadata["future_installed_challenge"]
    if set(challenge["five_decision_layers"]) != {
        "artifact_input",
        "solver",
        "numerical",
        "physical",
        "predictive",
    } or set(challenge["five_decision_layers"].values()) != {"not_run"}:
        raise ValueError("five pre-model decision layers changed")
    if challenge["globality_certificate"] != "not_guaranteed":
        raise ValueError("globality classification changed")

    verified_sources: list[str] = []
    if verify_local_sources:
        for name, source in metadata["primary_sources"].items():
            for prefix, hash_key in (
                ("", "sha256"),
                ("authority_pdf_", "authority_pdf_sha256"),
            ):
                path = Path(source[f"{prefix}execution_locator"])
                if not path.is_file() or sha256(path) != source[hash_key]:
                    raise ValueError(
                        f"local source verification failed: {name}:{prefix}"
                    )
                verified_sources.append(f"{name}:{prefix or 'markdown'}")
        for name, source in metadata["bound_parameter_sources"].items():
            for prefix in ("pdf_", "permanent_lab_pdf_", "markdown_", "mathpix_"):
                hash_key = f"{prefix}sha256"
                if hash_key not in source:
                    continue
                path = Path(source[f"{prefix}execution_locator"])
                if not path.is_file() or sha256(path) != source[hash_key]:
                    raise ValueError(
                        f"local source verification failed: {name}:{prefix}"
                    )
                verified_sources.append(f"{name}:{prefix}")

    return {
        "status": "source_complete_for_bounded_ascani_epcsaft_advanced_model",
        "case_id": row["case_id"],
        "rows": len(rows),
        "csv_sha256": csv_hash,
        "metadata_sha256": sha256(metadata_path),
        "primary_source_sha256": source_hashes,
        "primary_authority_pdf_sha256": authority_pdf_hashes,
        "parameter_source_sha256": parameter_source_hashes,
        "missing_upstream_sources": len(metadata["missing_upstream_provenance"]),
        "d026_unresolved_ascani_source_records": len(
            source_screen["unresolved_ascani_source_records"]
        ),
        "provider_status": source_screen["provider_status"],
        "local_sources_verified": verified_sources,
        "model_output": "not_run",
        "globality_certificate": challenge["globality_certificate"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate the frozen source-only Ascani 2022 HELD2 tracer contract."
    )
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--verify-local-sources", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            check(args.csv, args.metadata, args.verify_local_sources),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
