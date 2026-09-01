from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "analyses/2026-mea-parameter-estimation"
PACKET = ROOT / "data" / "packets" / "mea-co2-h2o-nine-species-estimation" / "1"
MATERIALIZED = CAMPAIGN / "inputs/data/mea-co2-h2o-nine-species-estimation/1"
SELECTED = CAMPAIGN / "inputs/neutral-selected-parameters.json"
ACTIVE = CAMPAIGN / "inputs/active-bundle.json"
SELECTION = CAMPAIGN / "results/selection-ladder.json"
FAILURE = CAMPAIGN / "results/failure-analysis.json"
LOCKS = ROOT / "data/packet-locks.json"
GROSS_PARAMETERS = (
    ROOT / "data" / "packets" / "gross-sadowski-2001" / "1" / "parameters/single.csv"
)
CAI_OBSERVATIONS = (
    CAMPAIGN
    / "inputs/source/library/observations/vapor_liquid_equilibrium/Cai_1996_MEA_H2O_VLE.csv"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def packet_fingerprint(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix().encode()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(path.stat().st_size.to_bytes(8, "big"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def status(identity: str) -> tuple[str, str]:
    if identity == "pair/carbon-dioxide/monoethanolamine/k_ij":
        return (
            "missing",
            "zero is an explicit placeholder; no neutral CO2-MEA data identify it",
        )
    if identity == "pair/monoethanolamine/water/k_ij":
        return "fitted", "Cai 1996 grouped 2B selection and complete 25-row refit"
    provisional = (
        "born_diameter",
        "debye_huckel_diameter",
        "packing_diameter",
        "suppression",
        "relative_permittivity",
        "ionic_region",
    )
    ionic = (
        "protonated-monoethanolamine",
        "carbamate-anion",
        "bicarbonate-anion",
        "carbonate-anion",
        "hydronium-cation",
        "hydroxide-anion",
    )
    if any(item in identity for item in provisional):
        return (
            "provisional",
            "not independently identified by MEA electrolyte-sensitive evidence",
        )
    if identity.startswith("pair/") and any(item in identity for item in ionic):
        return (
            "provisional",
            "inherited ionic interaction; pressure did not uniquely identify it",
        )
    if identity.startswith("component/") and any(item in identity for item in ionic):
        if identity.endswith(("charge_number", "molar_mass")):
            return "fixed", "verified species structure"
        return "provisional", "inherited ionic molecular value; not a hydrated diameter"
    if identity.startswith("reaction:"):
        return (
            "fixed",
            "source-audited baseline correlation; R4/R5 above 323 K are extrapolated",
        )
    if identity.startswith("association/carbon-dioxide/"):
        return "candidate", "Schick/Pabsch induced CO2-water association baseline"
    return "fixed", "source or selected neutral structure"


def coordinate(
    identity: str, value: Any, unit: str, record: dict[str, Any]
) -> dict[str, Any]:
    state, rationale = status(identity)
    provenance = record.get("provenance") or record.get("source") or {}
    return {
        "identity": identity,
        "value": value,
        "unit": unit,
        "status": state,
        "rationale": rationale,
        "candidate_domain_id": record.get("candidate_domain_id"),
        "qualification": record.get("qualification"),
        "provenance": provenance,
    }


def coordinates(parameters: dict[str, Any]) -> list[dict[str, Any]]:
    values: list[dict[str, Any]] = []
    for component in parameters["components"]:
        component_id = component["component_id"]
        for family, record in component["fixed"].items():
            value = record["value"]
            values.append(
                coordinate(
                    f"component/{component_id}/{family}",
                    value["magnitude"],
                    value["unit"],
                    record,
                )
            )
        for record in component["coefficients"]:
            value = record["value"]
            values.append(
                coordinate(
                    record["identity"], value["magnitude"], value["unit"], record
                )
            )
    for pair in parameters["pairs"]:
        for record in pair["coefficients"]:
            value = record["value"]
            values.append(
                coordinate(
                    record["identity"], value["magnitude"], value["unit"], record
                )
            )
    for record in parameters["model_coefficients"]:
        value = record["value"]
        values.append(
            coordinate(record["identity"], value["magnitude"], value["unit"], record)
        )
    for family in parameters["model_families"]:
        values.append(
            coordinate(family["family_id"], family["choice"], "choice", family)
        )
    for site in parameters["topology"]["sites"]:
        identity = f"association/{site['component_id']}/{site['site_id']}/multiplicity"
        values.append(coordinate(identity, site["multiplicity"], "count", site))
    for edge in parameters["topology"]["edges"]:
        for field in ("energy_over_k", "volume"):
            if field not in edge:
                continue
            record = edge[field]
            value = record["value"]
            values.append(
                coordinate(record["identity"], value["magnitude"], value["unit"], edge)
            )
    for reaction in parameters["reaction_correlations"]:
        for record in reaction["coefficients"]:
            values.append(
                coordinate(
                    record["identity"],
                    record["value"]["magnitude"],
                    record["value"]["unit"],
                    reaction,
                )
            )
    identities = [item["identity"] for item in values]
    if len(identities) != len(set(identities)):
        raise RuntimeError("active bundle contains duplicate coordinate identities")
    return sorted(values, key=lambda item: item["identity"])


def promote_parameters() -> dict[str, Any]:
    selection = json.loads(SELECTION.read_text())
    expected = selection["neutral_association"][
        "selected_parameters_sha256"
    ].removeprefix("sha256:")
    if (
        sha256(SELECTED) != expected
        or selection["neutral_association"]["selected_scheme"] != "2b"
    ):
        raise RuntimeError("retained neutral selection identity changed")
    parameters = json.loads(SELECTED.read_text())
    co2 = next(
        item
        for item in parameters["components"]
        if item["component_id"] == "carbon-dioxide"
    )
    segment_count = next(
        item for item in co2["coefficients"] if item["family"] == "segment_count"
    )
    segment_count["value"]["magnitude"] = 2.0729
    segment_count["source_sha256"] = f"sha256:{sha256(GROSS_PARAMETERS)}"
    segment_count["provenance"] = {
        "source_id": "gross-sadowski-2001",
        "domain_id": "carbon-dioxide-fit-range",
        "locator": "gross-sadowski-2001/1/parameters/single.csv:carbon-dioxide-segment-count",
    }
    domains = {item["domain_id"]: item for item in parameters["domains"]}
    if "carbon-dioxide-fit-range" not in domains:
        parameters["domains"].append(
            {
                "domain_id": "carbon-dioxide-fit-range",
                "kind": "fit-range",
                "temperature_min": {"magnitude": 216.0, "unit": "kelvin"},
                "temperature_max": {"magnitude": 304.0, "unit": "kelvin"},
            }
        )
    cai_domain = domains["cai-1996-neutral-refit"]
    cai_domain.clear()
    cai_domain.update(
        {
            "domain_id": "cai-1996-neutral-refit",
            "kind": "fit-range",
            "temperature_min": {"magnitude": 362.81, "unit": "kelvin"},
            "temperature_max": {"magnitude": 431.85, "unit": "kelvin"},
            "pressure_min": {"magnitude": 66660.0, "unit": "pascal"},
            "pressure_max": {"magnitude": 101330.0, "unit": "pascal"},
        }
    )
    mea_water = next(
        record
        for pair in parameters["pairs"]
        for record in pair["coefficients"]
        if record["identity"] == "pair/monoethanolamine/water/k_ij"
    )
    mea_water["source_sha256"] = f"sha256:{sha256(CAI_OBSERVATIONS)}"
    mea_water["provenance"]["domain_id"] = "mea-diagnostic-293-15-to-393-15-k"
    mea_water["provenance"]["locator"] = (
        "grouped pressure-level one-SE selection and complete 25-row refit; "
        "observation bounds are domain cai-1996-neutral-refit"
    )
    if not any(
        item["source_id"] == "gross-sadowski-2001" for item in parameters["sources"]
    ):
        parameters["sources"].append(
            {
                "source_id": "gross-sadowski-2001",
                "citation": "Gross and Sadowski (2001), Table 2, carbon dioxide row.",
                "use_basis": "CO2 m=2.0729; 2.079 is retained only as a retired historical calculation value.",
            }
        )
    cai_source = next(
        item
        for item in parameters["sources"]
        if item["source_id"] == "cai-1996-neutral-refit"
    )
    cai_source["use_basis"] = (
        "Provisional project fit of MEA-water k_ij to 25 admitted Cai (1996) rows; "
        "not an accepted reactive parameter set."
    )
    runtime = PACKET / "parameters/parameters.json"
    write_json(runtime, parameters)
    return parameters


def update_packet() -> str:
    runtime = PACKET / "parameters/parameters.json"
    provenance = PACKET / "provenance/source-origin.toml"
    provenance.write_text(
        f'''provenance_kind = "derived-record-provenance-v1"
schema_version = 1
transformations = ["2B neutral-family selection", "Cai (1996) MEA-water fit", "Gross-Sadowski (2001) CO2 segment-count correction"]

[[records]]
retained_path = "parameters/parameters.json"
sha256 = "{sha256(runtime)}"

[[sources]]
source_id = "neutral-selected-parameters"
role = "2B-selected base parameter document and MEA-water fit"
retained_path = "analyses/2026-mea-parameter-estimation/inputs/neutral-selected-parameters.json"
sha256 = "{sha256(SELECTED)}"

[[sources]]
role = "CO2 segment-count source"
source_id = "gross-sadowski-2001"
packet_id = "gross-sadowski-2001"
packet_version = "1"
packet_member = "parameters/single.csv"
sha256 = "{sha256(GROSS_PARAMETERS)}"

[[sources]]
role = "MEA-water fit observations"
source_id = "cai-1996-neutral-refit"
retained_path = "analyses/2026-mea-parameter-estimation/inputs/source/library/observations/vapor_liquid_equilibrium/Cai_1996_MEA_H2O_VLE.csv"
sha256 = "{sha256(CAI_OBSERVATIONS)}"
''',
        encoding="utf-8",
    )
    manifest = PACKET / "packet.toml"
    text = manifest.read_text()
    text = re.sub(
        r'(path = "parameters/parameters.json"\nrole = "parameters"\nmedia_type = "application/json"\nsha256 = ")[0-9a-f]{64}',
        rf"\g<1>{sha256(runtime)}",
        text,
    )
    text = text.replace(
        'source_refs = ["mea-thermodynamics-authority", "reaction-primary-records"]',
        'source_refs = ["mea-thermodynamics-authority", "reaction-primary-records", "gross-sadowski-2001", "cai-1996-neutral-refit"]',
    ).replace(
        'source_refs = ["mea-thermodynamics-authority"]',
        'source_refs = ["mea-thermodynamics-authority", "gross-sadowski-2001", "cai-1996-neutral-refit"]',
    )
    text = re.sub(
        r'(path = "provenance/source-origin.toml"\nrole = "provenance"\nmedia_type = "application/toml"\nsha256 = ")[0-9a-f]{64}',
        rf"\g<1>{sha256(provenance)}",
        text,
    )
    manifest.write_text(text, encoding="utf-8")
    for relative in (
        Path("packet.toml"),
        Path("parameters/parameters.json"),
        Path("provenance/source-origin.toml"),
    ):
        target = MATERIALIZED / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((PACKET / relative).read_bytes())
    fingerprint = packet_fingerprint(PACKET)
    locks = json.loads(LOCKS.read_text())
    lock = next(
        item
        for item in locks["packets"]
        if item["packet_id"] == "mea-co2-h2o-nine-species-estimation"
    )
    lock["packet_fingerprint"] = fingerprint
    write_json(
        LOCKS,
        {
            "schema_version": locks["schema_version"],
            "repository": locks["repository"],
            "packet_root": locks["packet_root"],
            "packets": [
                {
                    "packet_id": item["packet_id"],
                    "packet_version": item["packet_version"],
                    "packet_path": item["packet_path"],
                    "packet_fingerprint": item["packet_fingerprint"],
                    "materialization": item["materialization"],
                }
                for item in locks["packets"]
            ],
        },
    )
    write_json(
        CAMPAIGN / "inputs/data/mea-co2-h2o-nine-species-estimation/materialization-v1.json",
        {
            "schema_version": 3,
            "repository": "tannerpolley/ePC-SAFT",
            "lock_manifest": "data/packet-locks.json",
            "packet_id": "mea-co2-h2o-nine-species-estimation",
            "packet_version": "1",
            "packet_path": "data/packets/mea-co2-h2o-nine-species-estimation/1",
            "packet_fingerprint": fingerprint,
            "materialization_mode": "byte-for-byte-validation-input",
            "destination": "analyses/2026-mea-parameter-estimation/inputs/data/mea-co2-h2o-nine-species-estimation/1",
            "files": [
                {"path": str(path), "sha256": sha256(PACKET / path)}
                for path in (
                    Path("packet.toml"),
                    Path("parameters/parameters.json"),
                    Path("provenance/source-origin.toml"),
                )
            ],
        },
    )
    return fingerprint


def main() -> None:
    parameters = promote_parameters()
    fingerprint = update_packet()
    failure = json.loads(FAILURE.read_text())
    fit_input = json.loads((CAMPAIGN / "inputs/fit-input.json").read_text())
    reaction = fit_input["observations"][0]["request"]["reaction_system"]
    active = {
        "schema_version": 1,
        "bundle_id": "mea-co2-h2o-nine-species-active",
        "bundle_version": 1,
        "lifecycle": "active_provisional",
        "scientific_status": "diagnostic_candidate",
        "parameter_document": "analyses/2026-mea-parameter-estimation/inputs/data/mea-co2-h2o-nine-species-estimation/1/parameters/parameters.json",
        "parameter_document_sha256": sha256(PACKET / "parameters/parameters.json"),
        "packet_fingerprint": fingerprint,
        "species_ids": reaction["species_ids"],
        "charges": reaction["charges"],
        "molar_masses_kg_per_mol": reaction["molar_masses_kg_per_mol"],
        "reaction_matrix_products_positive": reaction["reaction_matrix"],
        "reaction_standard_state": reaction["source_standard_state"],
        "domains": {
            "computational_temperature_k": [293.15, 393.15],
            "co2_source_fit_temperature_k": [216.0, 304.0],
            "cai_fit_temperature_k": [362.81, 431.85],
            "cai_fit_pressure_pa": [66660.0, 101330.0],
            "source_qualified": "record-specific; Provider result is 313.15 K and R4/R5 above 323 K are extrapolated",
            "validated_pressure": None,
            "validated_speciation": None,
        },
        "coordinates": coordinates(parameters),
        "model_alternatives": [
            {
                "identity": "original-born",
                "status": "candidate",
                "disposition": "active inherited baseline; reactive N=8 promotion failed",
            },
            {
                "identity": "ssm",
                "status": "rejected",
                "disposition": "no independent MEA electrolyte-sensitive evidence earned replacement",
            },
            {
                "identity": "ssm-plus-ds",
                "status": "rejected",
                "disposition": "added suppression complexity not independently identified",
            },
        ],
        "rejected_alternatives": [
            {
                "identity": f"reactive-candidate-{index}",
                "status": "rejected",
                "reason": "mixed N=8 certification/convergence gate failed",
            }
            for index in range(1, 5)
        ],
        "unsupported": [
            "calorimetry",
            "unique ionic-coordinate identification",
            "downstream transfer",
        ],
        "terminal_failure": {
            "runtime_seconds": 40133.06,
            "solver_status": "maximum_iterations_exceeded",
            "failed_gates": failure["failed_gates"],
            "source_commit": "f25d9c97bdf3f38ff806065f45d060e3cd4a0bfb",
        },
    }
    write_json(ACTIVE, active)
    print(
        json.dumps(
            {
                "active_bundle_sha256": sha256(ACTIVE),
                "packet_fingerprint": fingerprint,
                "coordinates": len(active["coordinates"]),
            }
        )
    )


if __name__ == "__main__":
    main()
