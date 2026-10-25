"""Write held-2014-strategy2/2 (`held`) or the stage-A candidate packets (`candidates`), their locks and the index."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
DATA = ROOT / "data"
sys.path.insert(0, str(DATA / "tools"))
from validate_packets import packet_fingerprint  # noqa: E402

CONVENTION = (
    "Eq.(21), p.2890 d_ion = 0.88 sigma_ion; closest approach inferred from the Figs. 1-2 reproduction "
    "(analyses/reproduction/2014-held), text Section 2.1 p.2886 says sigma"
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def value(item) -> str:
    return json.dumps(item, ensure_ascii=False) if not isinstance(item, (int, float)) else repr(item)


def table(name: str, entries: dict) -> str:
    return f"[[{name}]]\n" + "".join(f"{key} = {value(item)}\n" for key, item in entries.items())


def write(packet_id: str, version: str, manifest: dict, sources: list, domains: list, parameters: dict,
          inputs: list[str], transformation: str) -> None:
    root = DATA / "packets" / packet_id / version
    (root / "parameters").mkdir(parents=True, exist_ok=True)
    (root / "provenance").mkdir(exist_ok=True)
    document = root / "parameters/parameters.json"
    document.write_text(json.dumps(parameters, indent=2))
    source_ids = [s["source_id"] for s in sources]
    domain_ids = [d["domain_id"] for d in domains]
    origin = root / "provenance/source-origin.toml"
    origin.write_text(
        'provenance_kind = "derived-record-provenance-v1"\nschema_version = 1\n'
        f"transformations = {value([transformation])}\n\n"
        + table("records", {"retained_path": "parameters/parameters.json", "sha256": sha(document)})
        + "".join("\n" + table("sources", {"source_id": f"input-{i}", "role": "Retained fit input or numerical result",
                                           "retained_path": path, "sha256": sha(ROOT / path)})
                  for i, path in enumerate(inputs)))
    files = [
        {"path": "parameters/parameters.json", "role": "parameters", "media_type": "application/json",
         "sha256": sha(document), "source_refs": source_ids, "domain_refs": domain_ids,
         "units_basis": "Explicit units per coefficient; sigma, DH, Born and packing lengths in angstrom, energies in kelvin, "
                        "segment counts and k_ij dimensionless."},
        {"path": "provenance/source-origin.toml", "role": "provenance", "media_type": "application/toml",
         "sha256": sha(origin), "source_refs": source_ids, "domain_refs": domain_ids,
         "units_basis": "Non-numeric source identities and SHA-256 hashes."},
    ]
    (root / "packet.toml").write_text(
        "".join(f"{key} = {value(item)}\n" for key, item in
                {"schema_version": 2, "packet_id": packet_id, "packet_version": version, **manifest}.items())
        + "".join("\n" + table("sources", s) for s in sources)
        + "".join("\n" + table("domains", d) for d in domains)
        + "".join("\n" + table("files", f) for f in files))
    locks_path = DATA / "packet-locks.json"
    locks = json.loads(locks_path.read_text())
    entry = {"packet_id": packet_id, "packet_version": version, "packet_path": f"data/packets/{packet_id}/{version}",
             "packet_fingerprint": packet_fingerprint(root), "materialization": "byte-for-byte"}
    rows = locks["packets"]
    same = [i for i, p in enumerate(rows) if p["packet_id"] == packet_id]
    at = next((i for i in same if rows[i]["packet_version"] == version), None)
    if at is None:
        rows.insert(same[-1] + 1 if same else len(rows), entry)
    else:
        rows[at] = entry
    locks_path.write_text(json.dumps(locks, indent=2) + "\n")


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def toml_tables(path: Path, name: str) -> list[dict]:
    import tomllib

    return tomllib.loads(path.read_text())[name]


def held() -> None:
    old = DATA / "packets/held-2014-strategy2/1"
    parameters = json.loads((old / "parameters/parameters.json").read_text())
    parameters["document_version"] = 2
    for component in parameters["components"]:
        coefficients = {c["family"]: c for c in component.get("coefficients", [])}
        if "debye_huckel_diameter" in coefficients:
            dh = coefficients["debye_huckel_diameter"]
            dh["value"]["magnitude"] = coefficients["packing_diameter"]["value"]["magnitude"]
            dh["provenance"]["locator"] = CONVENTION
    manifest = {
        "title": "Held 2014 strategy 2 chloride parameters at 298.15 K, Eq. (21) closest approach",
        "description": (
            "Packet /1 with the Debye-Hueckel closest approach of all six ions set to the Eq. (21) diameter 0.88 sigma, "
            "the convention that reproduces the authors' Figs. 1-2 (1:1 salts; analyses/reproduction/2014-held); applying "
            "it to Mg2+ and Ca2+ is a model-level convention. Section 2.1 says sigma; /1 keeps that transcription. "
            "k_ij(Na+/water) 0.00045485 and k_ij(K+/water) 0.1997222 are the Table 2 footnote a/b laws evaluated at "
            "298.15 K only; FeOs's constant 0.0045 is a factor-10 evaluation of footnote a. Water sigma(T) uses "
            "exp(-0.01146 T) (Held 2008 Eq. 7), correcting the printed sign. Relative permittivity 78.09 is a declared "
            "Figiel 2025 supplement."),
        "license_or_use_basis": "Numeric source transcription and project-authored metadata for private scientific "
                                "verification; no full article content.",
        "tolerance_rationale": "Transcription with a settled convention, not a fitted or qualified set; the mixed "
                               "NaCl-MgCl2/NaCl-CaCl2 water properties fail the #556 limits (analyses/2026-chloride-"
                               "brine-refit). No numerical acceptance tolerance is claimed.",
    }
    write("held-2014-strategy2", "2", manifest, toml_tables(old / "packet.toml", "sources"),
          toml_tables(old / "packet.toml", "domains"), parameters,
          ["data/packets/held-2014-strategy2/1/parameters/parameters.json",
           "analyses/reproduction/2014-held/results/figure-line-comparison.csv"],
          "Packet /1 with debye_huckel_diameter = packing_diameter (0.88 sigma) for every ion; no fitted value.")


FIT = "issue-575-stage-a-joint-fit"
FIT_DOMAIN = "issue-575-chloride-brine-298-15-k"


def fitted_document(spec: dict, labels: list[str], physical: list[float], start: str, use_basis: str) -> dict:
    """The base packet's document with every fitted coordinate replaced once; the single substitution owner."""
    parameters = json.loads((DATA / "packets" / spec["base_packet"] / "parameters/parameters.json").read_text())
    parameters["document_id"], parameters["document_version"] = spec["candidate_packet"], 1
    fitted = dict(zip(labels, physical))
    provenance = {"source_id": FIT, "locator": f"stage-A joint fit, start {start} (analyses/2026-chloride-brine-refit)",
                  "domain_id": FIT_DOMAIN}
    stored = 0
    for owner, coefficients in [(c["component_id"], c.get("coefficients", [])) for c in parameters["components"]] + [
            ("+".join(sorted((q["component_id_a"], q["component_id_b"]))), q["coefficients"]) for q in parameters["pairs"]]:
        for c in coefficients:
            if f"{c['family']}/{owner}" in fitted:
                c["value"]["magnitude"] = fitted[f"{c['family']}/{owner}"]
                c["provenance"] = provenance
                stored += 1
    assert stored == len(fitted), "every fitted coordinate is stored exactly once"
    parameters["sources"].append({"source_id": FIT, "citation": "ePC-SAFT issue 575 stage-A native joint fit, 2026-10-16",
                                  "use_basis": use_basis})
    parameters["domains"].append({"domain_id": FIT_DOMAIN, "kind": "reported-conditions",
                                  "temperature_min": {"magnitude": 298.15, "unit": "kelvin"},
                                  "temperature_max": {"magnitude": 298.15, "unit": "kelvin"},
                                  "pressure_min": {"magnitude": 1, "unit": "bar"}, "pressure_max": {"magnitude": 1, "unit": "bar"}})
    return parameters


def candidates() -> None:
    declaration = json.loads((HERE / "inputs/fit-declaration.json").read_text())
    run = f"analyses/2026-chloride-brine-refit/results/{declaration['run']}"
    selection = json.loads((ROOT / run / "selection.json").read_text())
    owner = declaration["held_out"]["owner_qualification"]
    for form, chosen in selection["forms"].items():
        spec = declaration["forms"][form]
        qualified = owner["run"] == declaration["run"] and owner["form"] == form
        if qualified:
            (miss,) = read_csv(ROOT / run / "held-out" / chosen["label"] / "dinane-mounir-failures.csv")
            rms, limit = float(miss["RMS"]), float(miss["limit_RMS"])
            folder = ROOT / run / "held-out" / chosen["label"]
            training = [(r, float(r["m_NaCl"]) + float(r["m_KCl"]) + 3 * float(r["m_divalent"]))
                        for r in read_csv(ROOT / run / f"training-residuals-{form}.csv")]
            held = [float(r["I_mol_kg"]) for r in read_csv(ROOT / run / f"nominal-residuals-{form}.csv")
                    if r["realization_id"] == "r0000" and r["role"] == "mixture"
                    and r["source"] in ("robinson-bower-1966", "dinane-mounir-2003")]
            above = max(abs(float(r["residual"])) for r, i in training if r["source"] == "rard-miller-1987" and i > max(held))
            mgcl2 = [float(r["residual"]) for r in read_csv(folder / "binary-reference-checks.csv")
                     if r["declaration"] == "candidate" and r["role"] == "MgCl2 binary"]
            phi = next(v for v in read_csv(folder / "verdicts.csv") if v["declaration"] == "candidate" and v["quantity"] == "Phi")
            scope = (
                f"Held-out evidence covers only the two-salt NaCl-CaCl2 (Robinson-Bower) and NaCl-MgCl2 (Dinane-Mounir) "
                f"mixtures at I = {min(held):.2f}-{max(held):.2f} mol/kg, against training at I = "
                f"{min(i for _, i in training):.3f}-{max(i for _, i in training):.2f} mol/kg, so ternary Na-Mg-Ca brines are a "
                f"composition extrapolation; above {max(held):.2f} mol/kg only training residuals and the Rard-Miller Table III "
                f"MgCl2 reference check exist, and for MgCl2-rich brines they exceed the frozen Phi limits (Rard-Miller |dPhi| up to "
                f"{above:.3f} against the {float(phi['limit_maximum']):.2f} maximum; pure-MgCl2 Phi RMS "
                f"{(sum(x * x for x in mgcl2) / len(mgcl2)) ** 0.5:.3f} against the {float(phi['limit_RMS']):.2f} RMS limit). ")
            status = (
                f"Qualified by owner decision {owner['date']} for aqueous NaCl-MgCl2-CaCl2 chloride brines at 298.15 K and 1 bar over "
                f"the ionic-strength range of its training and held-out data (domain {FIT_DOMAIN}). Qualification limit: passes every "
                f"held-out bucket in {chosen['realizations_passing_every_bucket']:,} of 1,467 printed-rounding realizations under the "
                f"source-agreement criterion; fails only {miss['realization_id']} (every Dinane-Mounir value lowered by half a printed "
                f"unit), where dinane-mounir-2003 {miss['I_bucket']} RMS {rms:.6f} exceeds its limit {limit:.6f} by "
                f"{100 * (rms / limit - 1):.1f} % and the maximum {float(miss['maximum_absolute_error']):.4f} is within "
                f"{float(miss['limit_maximum']):.4f}. " + scope + "Refit: "
            )
        else:
            status = "Unqualified stage-A candidate: "
        summary = (
            status + "Na+, Mg2+ and Ca2+ refit jointly on "
            "Hamer-Wu NaCl, Goldberg-Nuttall MgCl2/CaCl2 (<= 3.3 mol/kg), Rard-Miller 1987 NaCl-MgCl2 and Saad-Padova-Marcus 1980 "
            f"NaCl-KCl-CaCl2 mixture osmotic coefficients at 298.15 K on {spec['base_packet']}; water, Cl-, K+, Li+ and every u_i unchanged. Held-out "
            "water properties under the frozen #556 limits, Dinane-Mounir judged by source agreement (owner amendment 2026-10-16: per "
            "realization, the larger of the frozen limit and the Rard-Miller/Dinane-Mounir disagreement) (RMS/max, robust over 1,467 "
            "realizations): "
            + "; ".join(f"{b['source']} {b['I_bucket']} {b['RMS']:.4f}/{b['maximum_absolute_error']:.4f} "
                        + ("pass" if b["robust_passed"] else f"fail ({b['passing_realizations']:,}/1,467 realizations pass; failing: "
                                                              f"{b['failing_realizations']})") for b in chosen["held_out"])
            + ". 298.15 K only; no K+/Li+ mixture, temperature or physical-validation claim.")
        sources = toml_tables(DATA / "packets" / spec["base_packet"] / "packet.toml", "sources") + [
            {"source_id": FIT, "citation": "ePC-SAFT issue 575 stage-A native joint fit, 2026-10-16", "use_basis": summary},
            {"source_id": "goldberg-nuttall-1978", "citation": "Goldberg and Nuttall, J. Phys. Chem. Ref. Data 1978, 7, 263.",
             "doi": "10.1063/1.555569", "use_basis": "Evaluated MgCl2 Table 17 and CaCl2 Table 20 osmotic coefficients <= 3.3 mol/kg"},
            {"source_id": "hamer-wu-1972", "citation": "Hamer and Wu, J. Phys. Chem. Ref. Data 1972, 1, 1047.",
             "doi": "10.1063/1.3253108", "use_basis": "Evaluated NaCl Table 16 osmotic coefficients"},
            {"source_id": "rard-miller-1987", "citation": "Rard and Miller 1987, NaCl-MgCl2 isopiestic mixtures; packet "
             "mixed-chloride-water-observations/1", "use_basis": "287 mixture osmotic coefficients, nominal #556 conversion"},
            {"source_id": "marcus-1991", "citation": "Marcus, Y. J. Chem. Soc., Faraday Trans. 1991, 87, 2995-2999.",
             "doi": "10.1039/FT9918702995", "use_basis": "Table 1 crystal diameters 2r as lower bounds of fitted diameters"},
            {"source_id": "saad-padova-marcus-1980", "citation": "Saad, Padova and Marcus, J. Solution Chem. 1980, 9, 259-268; packet "
             "saad-padova-marcus-1980-nacl-kcl-cacl2-observations/1", "doi": "10.1007/BF00644551",
             "use_basis": "56 balance-consistent NaCl-KCl-CaCl2 mixture osmotic coefficients (training; K+ at settled values)"}]
        domains = [{"domain_id": FIT_DOMAIN, "kind": "reported-conditions", "temperature_min": 298.15, "temperature_max": 298.15,
                    "temperature_unit": "kelvin", "pressure_min": 1, "pressure_max": 1, "pressure_unit": "bar",
                    "description": "Aqueous NaCl-MgCl2-CaCl2 at 298.15 K and 1 bar, ionic strength <= 9.9 mol/kg."}]
        manifest = {"title": f"{spec['title']} (issue #575 stage A)", "description": summary,
                    "license_or_use_basis": "Project-derived parameters for private scientific verification; no full paper content.",
                    "tolerance_rationale": "Held-out errors under the frozen #556 limits, with the owner's 2026-10-16 "
                                           "source-agreement limit for Dinane-Mounir, are stated in the description"
                                           + ("; qualified by the owner with the one failing realization stated there as its limit"
                                              if qualified else "")
                                           + "; numerical fitting, no physical validation."}
        write(spec["candidate_packet"], "1", manifest, list({s["source_id"]: s for s in sources}.values()), domains,
              fitted_document(spec, chosen["labels"], chosen["physical"], chosen["start"], summary),
              [f"data/packets/{spec['base_packet']}/parameters/parameters.json",
               "analyses/2026-chloride-brine-refit/inputs/fit-declaration.json",
               f"{run}/fit-{form}.json", f"{run}/selection.json"],
              f"Native Ceres multistart joint fit ({form} form); the lowest-training-cost converged optimum is stored.")


if __name__ == "__main__":
    {"held": held, "candidates": candidates}[sys.argv[1]]()
    index = subprocess.run([sys.executable, DATA / "tools/parameter_index.py"], capture_output=True, text=True, check=True)
    (DATA / "docs/parameter-index.md").write_text(index.stdout)
