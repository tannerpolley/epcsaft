from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CAMPAIGN = Path(__file__).resolve().parents[1]
LIBRARY = CAMPAIGN / "inputs/source" / "library"
BASELINE = (
    ROOT / "analyses/2026-mea-grepe-regression/inputs/source/parameters.json"
)
PACKET = CAMPAIGN / "inputs/data/mea-co2-h2o-nine-species-estimation/1"
PARAMETERS = PACKET / "parameters/parameters.json"
SOURCE_COMMIT = "412f85b4a34d713aa3b00985f6bef044c1f12327"
SOURCE_DOMAIN = "mea-provider-source-qualified-313-15-k"
CANDIDATE_DOMAIN = "mea-candidate-293-15-to-393-15-k"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def build_parameters() -> str:
    source_hash = sha256(BASELINE)
    candidate = {
        "candidate_domain_id": CANDIDATE_DOMAIN,
        "qualification": "candidate_extrapolation",
        "source_sha256": f"sha256:{source_hash}",
    }
    value = json.loads(BASELINE.read_text(encoding="utf-8"))
    value["document_id"] = "mea-co2-h2o-nine-species-estimation-candidate-v1"
    value["sources"].extend(
        [
            {
                "source_id": "tong-2012-aroua-1999",
                "citation": "Tong et al. (2012), Table 5 K5 MEA, via Aroua et al. (1999).",
                "use_basis": "R4 carbamate hydrolysis source-centred correlation.",
            },
            {
                "source_id": "bates-pinching-1951",
                "citation": "Bates and Pinching (1951), equation 7.",
                "use_basis": "R5 monoethanolammonium dissociation source-centred correlation.",
            },
        ]
    )
    value["domains"] = [
        {
            "domain_id": "mea-diagnostic-293-15-to-393-15-k",
            "kind": "fit-range",
            "temperature_min": {"magnitude": 293.15, "unit": "kelvin"},
            "temperature_max": {"magnitude": 393.15, "unit": "kelvin"},
            "pressure_min": {"magnitude": 10.0, "unit": "pascal"},
            "pressure_max": {"magnitude": 300000.0, "unit": "pascal"},
        },
        {
            "domain_id": SOURCE_DOMAIN,
            "kind": "reported-conditions",
            "temperature_min": {"magnitude": 313.15, "unit": "kelvin"},
            "temperature_max": {"magnitude": 313.15, "unit": "kelvin"},
            "pressure_min": {"magnitude": 100000.0, "unit": "pascal"},
            "pressure_max": {"magnitude": 100000.0, "unit": "pascal"},
        },
        {
            "domain_id": CANDIDATE_DOMAIN,
            "kind": "fit-range",
            "temperature_min": {"magnitude": 293.15, "unit": "kelvin"},
            "temperature_max": {"magnitude": 393.15, "unit": "kelvin"},
            "pressure_min": {"magnitude": 1.0, "unit": "pascal"},
            "pressure_max": {"magnitude": 10000000.0, "unit": "pascal"},
        },
    ]
    value["model_families"] = [
        family for family in value["model_families"] if family["kind"] != "polar"
    ]
    for family in value["model_families"]:
        family["provenance"]["domain_id"] = SOURCE_DOMAIN
        family.update(candidate)
    for component in value["components"]:
        for fixed in component["fixed"].values():
            fixed.update(candidate)
        for coefficient in component["coefficients"]:
            coefficient.update(candidate)
    for pair in value["pairs"]:
        for coefficient in pair["coefficients"]:
            coefficient.update(candidate)
    for coefficient in value["model_coefficients"]:
        coefficient.update(candidate)
    for correlation in value["correlations"]:
        correlation.update(candidate)
    for record_kind in ("presets", "sites", "edges"):
        for record in value["topology"][record_kind]:
            record.update(candidate)
    # The source copy has no CO2/MEA cross edge; new documents declare it (#569), not a physical-absence claim.
    value["topology"]["no_cross_association"] = [{
        "component_id_a": "carbon-dioxide", "component_id_b": "monoethanolamine",
        "provenance": {"source_id": "mea-retained-phase-two-artifact", "locator": "baseline copy: no CO2/MEA edge"}}]
    value["reaction_correlations"] = [
        {
            "reaction_id": "R4",
            "kind": "ln-k-a-plus-b-over-t",
            "source": {
                "source_id": "tong-2012-aroua-1999",
                "locator": "Tong 2012 Table 5 K5 MEA; Aroua 1999 infinite-dilution relation",
                "domain_id": SOURCE_DOMAIN,
            },
            "candidate_domain_id": CANDIDATE_DOMAIN,
            "qualification": "candidate_extrapolation",
            "source_sha256": "sha256:78af3d67d58ff610dcc7c34fd55a75763470b7cd5bac0f3e39c256c7d7ffd61c",
            "coefficients": [
                {
                    "name": "a",
                    "identity": "reaction:R4:correlation:a",
                    "value": {"magnitude": 2.151, "unit": "dimensionless"},
                    "optimizer_domain": {"lower": -20.0, "upper": 20.0},
                },
                {
                    "name": "b_k",
                    "identity": "reaction:R4:correlation:b_k",
                    "value": {"magnitude": -1545.3, "unit": "kelvin"},
                    "optimizer_domain": {"lower": -10000.0, "upper": 5000.0},
                },
            ],
        },
        {
            "reaction_id": "R5",
            "kind": "negative-log10-temperature-polynomial",
            "source": {
                "source_id": "bates-pinching-1951",
                "locator": "Equation 7, p. 349",
                "domain_id": SOURCE_DOMAIN,
            },
            "candidate_domain_id": CANDIDATE_DOMAIN,
            "qualification": "candidate_extrapolation",
            "source_sha256": "sha256:aff8621efd41dbce106bb189fd104840a1513d9d163f4f57dfb1abbf53a30db1",
            "coefficients": [
                {
                    "name": "a_k",
                    "identity": "reaction:R5:correlation:a_k",
                    "value": {"magnitude": 2677.91, "unit": "kelvin"},
                    "optimizer_domain": {"lower": 500.0, "upper": 6000.0},
                },
                {
                    "name": "b",
                    "identity": "reaction:R5:correlation:b",
                    "value": {"magnitude": 0.3869, "unit": "dimensionless"},
                    "optimizer_domain": {"lower": -10.0, "upper": 10.0},
                },
                {
                    "name": "c_per_k",
                    "identity": "reaction:R5:correlation:c_per_k",
                    "value": {"magnitude": 0.0004277, "unit": "1 / kelvin"},
                    "optimizer_domain": {"lower": -0.01, "upper": 0.01},
                },
            ],
        },
    ]
    write_json(PARAMETERS, value)
    return source_hash


def build_admission() -> None:
    cai = rows(
        LIBRARY / "observations/vapor_liquid_equilibrium/Cai_1996_MEA_H2O_VLE.csv"
    )
    kiepe = rows(LIBRARY / "observations/physical_co2_water/kiepe_2002_table_1.csv")
    vle = rows(
        LIBRARY / "observations/vapor_liquid_equilibrium/Canonical_VLE_Observations.csv"
    )
    speciation = rows(LIBRARY / "manifests/speciation_target_membership.csv")
    pressure_candidates = Counter(
        row["source_key"]
        for row in vle
        if row["source_key"] in {"Hilliard2008", "Jou1995", "Xu2011"}
        and (row["source_key"] != "Hilliard2008" or row["MEA_weight_fraction"] == "0.3")
        and (row["source_key"] != "Jou1995" or row["active_view_member"] == "yes")
        and (row["source_key"] != "Xu2011" or row["active_view_member"] == "yes")
    )
    speciation_targets = Counter(
        row["source_key"] for row in speciation if row["target_eligible"] == "yes"
    )
    bottinger_states = {
        row["state_id"]
        for row in speciation
        if row["source_key"] == "Bottinger2008" and row["mea_mass_fraction"] == "0.3"
    }
    counts = {
        "cai_neutral_rows": len(cai),
        "kiepe_physical_rows": sum(row["included"].lower() == "true" for row in kiepe),
        "hilliard_30wt_pressure_rows": pressure_candidates["Hilliard2008"],
        "jou_pressure_rows": pressure_candidates["Jou1995"],
        "xu_untouched_pressure_rows": pressure_candidates["Xu2011"],
        "bottinger_eligible_targets": speciation_targets["Bottinger2008"],
        "bottinger_30wt_states": len(bottinger_states),
        "matin_eligible_targets": speciation_targets["Matin2012"],
        "jakobsen_untouched_targets": speciation_targets["Jakobsen2005"],
        "raw_calorimetry_rows": 158,
    }
    expected = {
        "cai_neutral_rows": 29,
        "kiepe_physical_rows": 39,
        "hilliard_30wt_pressure_rows": 31,
        "jou_pressure_rows": 48,
        "xu_untouched_pressure_rows": 18,
        "bottinger_eligible_targets": 55,
        "bottinger_30wt_states": 39,
        "matin_eligible_targets": 76,
        "jakobsen_untouched_targets": 67,
        "raw_calorimetry_rows": 158,
    }
    if counts != expected:
        raise RuntimeError(f"admission count drift: {counts}")
    write_json(
        CAMPAIGN / "inputs/source/admission.json",
        {
            "schema_version": 1,
            "identity": "mea-pressure-speciation-admission-v1",
            "source_authority_commit": SOURCE_COMMIT,
            "supersedes": "old candidate lifecycle and active-v1 row-number cutoffs",
            "selection_rule": "source fields, 30 wt% reactive membership, and frozen target eligibility",
            "counts": counts,
            "pressure_scale": {
                "value": 0.05,
                "unit": "natural_log",
                "role": "diagnostic_model_discrepancy_not_measurement_uncertainty",
            },
            "speciation_scale": {
                "value": 0.01,
                "unit": "absolute_mole_fraction",
                "role": "diagnostic_model_discrepancy_not_measurement_uncertainty",
            },
            "heat_status": "unsupported_source_state_current_engine_scope",
        },
    )

    groups = sorted(
        {
            (
                row["source_key"],
                row["temperature_canonical_C"] or row["temperature_reported_C"],
            )
            for row in vle
            if row["source_key"] in {"Hilliard2008", "Jou1995"}
            and (
                row["source_key"] != "Hilliard2008"
                or row["MEA_weight_fraction"] == "0.3"
            )
            and (row["source_key"] != "Jou1995" or row["active_view_member"] == "yes")
        }
        | {
            (row["source_key"], row["temperature_C"])
            for row in speciation
            if row["source_key"] in {"Bottinger2008", "Matin2012"}
            and row["target_eligible"] == "yes"
        }
    )
    write_json(
        CAMPAIGN / "inputs/source/folds.json",
        {
            "schema_version": 1,
            "grouping_rule": "source-temperature-series-indivisible",
            "fold_count": 5,
            "groups": [
                {
                    "group_id": f"{source}|{temperature}C",
                    "source": source,
                    "temperature_c": float(temperature),
                    "fold": index % 5,
                }
                for index, (source, temperature) in enumerate(groups)
            ],
        },
    )


def build_packet_manifest(source_hash: str) -> None:
    provenance = PACKET / "provenance/source-origin.toml"
    provenance.parent.mkdir(parents=True, exist_ok=True)
    provenance.write_text(
        "\n".join(
            [
                'provenance_kind = "source-record-provenance-v1"',
                "schema_version = 1",
                "",
                "[[records]]",
                'retained_path = "parameters/parameters.json"',
                'origin_repository = "tannerpolley/MEA-Thermodynamics"',
                f'origin_commit = "{SOURCE_COMMIT}"',
                f'origin_tree = "{SOURCE_COMMIT}"',
                'origin_path = "data/reference/epcsaft_bundles/mea-co2-h2o-nine-species-regression-input/1"',
                f'sha256 = "{sha256(PARAMETERS)}"',
                "",
            ]
        ),
        encoding="utf-8",
    )
    manifest = PACKET / "packet.toml"
    manifest.write_text(
        "\n".join(
            [
                "schema_version = 2",
                'packet_id = "mea-co2-h2o-nine-species-estimation"',
                'packet_version = "1"',
                'title = "Nine-species CO2-MEA-H2O estimation candidate"',
                'description = "Complete reusable nine-species ePC-SAFT baseline with fail-closed source and candidate model domains and typed R4/R5 correlations."',
                'license_or_use_basis = "Project-authorized numeric parameter records and source metadata from the pinned MEA-Thermodynamics authority commit; no publisher full text is included."',
                'tolerance_rationale = "The packet asserts no fit tolerance. Source values are immutable inputs; the consuming analysis owns diagnostic scales, scoring, and candidate status."',
                "",
                "[[sources]]",
                'source_id = "mea-thermodynamics-authority"',
                f'citation = "tannerpolley/MEA-Thermodynamics commit {SOURCE_COMMIT}."',
                'use_basis = "Complete nine-species baseline, component and model records."',
                "",
                "[[sources]]",
                'source_id = "reaction-primary-records"',
                'citation = "Tong et al. (2012)/Aroua et al. (1999) R4 and Bates and Pinching (1951) R5 records."',
                'use_basis = "Source-centred R4/R5 correlation coefficients and original source-qualified ranges."',
                "",
                "[[domains]]",
                f'domain_id = "{SOURCE_DOMAIN}"',
                'description = "Original installed Provider source-transfer executable result at 313.15 K and 100 kPa; not the candidate support claim."',
                "",
                "[[domains]]",
                f'domain_id = "{CANDIDATE_DOMAIN}"',
                'description = "Candidate evaluation domain 293.15-393.15 K and 1 Pa-10 MPa, qualified candidate_extrapolation pending campaign evidence."',
                "",
                "[[files]]",
                'path = "parameters/parameters.json"',
                'role = "parameters"',
                'media_type = "application/json"',
                f'sha256 = "{sha256(PARAMETERS)}"',
                'source_refs = ["mea-thermodynamics-authority", "reaction-primary-records"]',
                f'domain_refs = ["{SOURCE_DOMAIN}", "{CANDIDATE_DOMAIN}"]',
                'units_basis = "Every coefficient has an explicit canonical unit; reaction coefficients retain dimensionless, K, and K^-1 units."',
                "",
                "[[files]]",
                'path = "provenance/source-origin.toml"',
                'role = "provenance"',
                'media_type = "application/toml"',
                f'sha256 = "{sha256(provenance)}"',
                'source_refs = ["mea-thermodynamics-authority"]',
                f'domain_refs = ["{SOURCE_DOMAIN}", "{CANDIDATE_DOMAIN}"]',
                f'units_basis = "Repository identity, pinned commit, source payload hash {source_hash}, and retained path are non-numeric provenance."',
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    source_hash = build_parameters()
    build_admission()
    build_packet_manifest(source_hash)
    print(json.dumps({"parameters_sha256": sha256(PARAMETERS), "admission": "clean"}))


if __name__ == "__main__":
    main()
