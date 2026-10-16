#!/usr/bin/env python3
"""One-time recipe (#569): four pinned FeOs PC-SAFT parameter files -> packets feos-*/1 and their locks.

Usage: python3 data/tools/import_feos.py <directory holding the four FeOs JSON files>
Every source record becomes an admitted Engine record or a refusal with its cause
(sources/record-accounting.csv). Decimal source values convert before binary64 rounding. No equation
of state is evaluated. Rerun data/tools/parameter_index.py afterwards."""
import csv, io, json, math, sys
from decimal import Decimal
from pathlib import Path

from validate_packets import ROOT, packet_fingerprint, sha256

REPO, COMMIT, TREE = "feos-org/feos", "af84e3f9da2fac2d102d449cd07cbbb48237617b", "1e7fb4770d9055b96d010e0b3f3f2d0241269508"
HASH = {"gross2001.json": "f4b4018c7f02341b937c086cbb38626b72c1f3f3d9677a26b463112e16dfd404",
        "gross2002.json": "9680ae8998565feabf6b3808bc2f11df288f55db049677ccebe08a937332c491",
        "esper2023.json": "87fdcd8733245d90cd6ad9e9d7d0583bd3f7ffb3c4d67a25197190f289e2a679",
        "rehner2023_binary.json": "8b7d008b2bdf7ec49542eb443e36078fbe46e2737775126589c1718dfa13c335"}
META = json.loads((Path(__file__).with_name("import_feos.json")).read_text())  # citations, licences, domains: data


class Refused(Exception):
    pass


def value(magnitude, unit):
    return {"magnitude": float(magnitude), "unit": unit}


def positive(*numbers):
    if not all(isinstance(x, Decimal) and x.is_finite() and x > 0 for x in numbers):
        raise Refused("missing, nonpositive or nonfinite parameter")


def pure(record, locator, source, polar):
    """Component, sites and self edge of one pure record; FeOs na are donor (a), nb acceptor (b) sites."""
    if unknown := set(record) - {"identifier", "molarweight", "m", "sigma", "epsilon_k", "association_sites", "mu"}:
        raise Refused(f"unconsumed field {','.join(sorted(unknown))}")  # e.g. a quadrupole: no polar expansion here
    cid = record["identifier"]["name"]
    positive(*(record.get(k) for k in ("molarweight", "m", "sigma", "epsilon_k")))
    prov = lambda field: {"source_id": source, "locator": f"{locator}.{field}", "domain_id": source}
    coef = lambda family, field, magnitude, unit: {"identity": f"component/{cid}/{family}", "family": family,
                                                  "value": value(magnitude, unit), "provenance": prov(field)}
    coefficients = [coef("segment_count", "m", record["m"], "dimensionless"), coef("segment_diameter", "sigma", record["sigma"], "angstrom"),
                    coef("dispersion_energy_over_k", "epsilon_k", record["epsilon_k"], "kelvin")]
    if polar:  # Gross-Vrabec 2006 via FeOs: one dipole per molecule; mu absent is nonpolar (Esper 2023 Sec. III)
        mu = record.get("mu", Decimal(0))
        if not (mu.is_finite() and mu >= 0):
            raise Refused("invalid dipole moment")
        coefficients += [coef("dipole_moment", "mu", mu, "debye"), coef("dipole_count", "mu", int(mu > 0), "dimensionless")]
    sites, edges, sigma = [], [], float(record["sigma"])
    if len(record.get("association_sites", [])) > 1:
        raise Refused("more than one association site type")
    for site in record.get("association_sites", []):
        na, nb = site.get("na", Decimal(0)), site.get("nb", Decimal(0))
        if set(site) - {"na", "nb", "kappa_ab", "epsilon_k_ab"} or na % 1 or nb % 1 or na < 0 or nb < 0 or na + nb == 0:
            raise Refused("unsupported site record")
        self_associating = "kappa_ab" in site or "epsilon_k_ab" in site
        if self_associating != (na > 0 < nb):
            raise Refused("association parameters without a donor/acceptor pair")
        sites += [{"component_id": cid, "site_id": s, "site_role": role, "multiplicity": int(n), "provenance": prov(f"association_sites.n{s}")}
                  for s, role, n in (("a", "donor", na), ("b", "acceptor", nb)) if n]
        if self_associating:
            positive(site.get("kappa_ab"), site.get("epsilon_k_ab"))
            edges.append(edge(cid, cid, site["epsilon_k_ab"], site["kappa_ab"], prov("association_sites")))
    fixed = {"molar_mass": {"value": value(record["molarweight"] / 1000, "kg / mol"), "provenance": prov("molarweight")},
             "charge_number": {"value": value(0, "elementary-charge"), "provenance": prov("identifier")}}
    return {"component_id": cid, "identifier": record["identifier"], "fixed": fixed, "coefficients": coefficients}, sites, edges, sigma


def edge(a, b, energy, volume, provenance, source=None):
    names = f"association/{a}/a/{b}/b"
    return {"endpoint_a": {"component_id": a, "site_id": "a"}, "endpoint_b": {"component_id": b, "site_id": "b"},
            **({} if source else {"energy_over_k": {"identity": f"{names}/energy_over_k", "value": value(energy, "kelvin")},
                                  "volume": {"identity": f"{names}/volume", "value": value(volume, "dimensionless")}}),
            "source": source or {"kind": "explicit", "provenance": [provenance]}}


def rule(a, b, provenance):
    inputs = {k: [f"association/{c}/a/{c}/b/{k}" for c in (a, b)] for k in ("energy_over_k", "volume")}
    return edge(a, b, 0, 0, None, {"kind": "combining-rule", "rule_id": "wolbach-sandler", "inputs": inputs, "provenance": [provenance]})


def convert(directory):
    loaded, rows, packets = {}, [], {}
    for name, digest in HASH.items():
        if sha256(directory / name) != digest:
            raise SystemExit(f"{name}: SHA-256 differs from the pinned FeOs {COMMIT} file")
        loaded[name] = json.loads((directory / name).read_text(), parse_float=Decimal, parse_int=Decimal)
    for name, packet_id in (("gross2001.json", "feos-gross-sadowski-2001"), ("gross2002.json", "feos-gross-sadowski-2002"),
                            ("esper2023.json", "feos-esper-2023")):
        admitted = {}
        for index, record in enumerate(loaded[name]):
            try:
                admitted[record["identifier"]["cas"]] = pure(record, f"parameters/pcsaft/{name}[{index}]", packet_id, name == "esper2023.json")
                rows.append((packet_id, name, index, record["identifier"]["name"], "admitted", ""))
            except Refused as cause:
                rows.append((packet_id, name, index, record["identifier"].get("name"), "refused", str(cause)))
        packets[packet_id] = (admitted, [], [], [])
    esper, gross2002 = packets["feos-esper-2023"][0], packets["feos-gross-sadowski-2002"]
    associating = [cas for cas, item in gross2002[0].items() if item[2]]
    provenance = {"source_id": "feos-gross-sadowski-2002", "locator": "Gross-Sadowski 2002 Eqs. 2-3 (Wolbach-Sandler)",
                  "domain_id": "feos-gross-sadowski-2002"}
    for i, a in enumerate(associating):  # explicit source rules; k_ij stays absent, so a binary selection refuses on it
        for b in associating[i + 1:]:
            ida, idb = gross2002[0][a][0]["component_id"], gross2002[0][b][0]["component_id"]
            gross2002[2].extend((rule(ida, idb, provenance), rule(idb, ida, provenance)))
    used, pairs, edges, none = {}, [], [], []
    for index, record in enumerate(loaded["rehner2023_binary.json"]):
        locator, source = f"parameters/pcsaft/rehner2023_binary.json[{index}]", "feos-rehner-2023"
        ends = [esper.get(record.get(k, {}).get("cas")) for k in ("id1", "id2")]
        try:
            if set(record) - {"id1", "id2", "k_ij", "association_sites"}:
                raise Refused("unconsumed binary field")
            if None in ends:
                raise Refused("endpoint absent from or refused in esper2023.json")
            if ("k_ij" in record) == ("association_sites" in record):
                raise Refused("expected exactly one fitted binary parameter")  # Rehner 2023 Sec. 3: k_ij or eps^AiBj
            (a, sa, ea, sig_a), (b, sb, eb, sig_b) = ends
            ida, idb = a["component_id"], b["component_id"]
            prov = lambda note: {"source_id": source, "locator": f"{locator}{note}", "domain_id": source}
            k = record.get("k_ij", Decimal(0))
            if not k.is_finite():
                raise Refused("nonfinite k_ij")
            note = ".k_ij" if "k_ij" in record else ": k_ij not adjusted when eps^AiBj is fitted (Rehner 2023 Sec. 3)"
            pair_edges, declared = [], None
            donors = [(x, y) for x, y, sx, sy in ((ida, idb, sa, sb), (idb, ida, sb, sa))
                      if any(s["site_id"] == "a" for s in sx) and any(s["site_id"] == "b" for s in sy)]
            if "association_sites" in record:
                [site] = record["association_sites"]
                if set(site) != {"kappa_ab", "epsilon_k_ab"} or not donors:
                    raise Refused("unsupported induced-association record")
                if site["epsilon_k_ab"] <= 0:
                    raise Refused("nonpositive fitted cross-association energy: Eq. 17 strength below zero")
                positive(site["kappa_ab"], site["epsilon_k_ab"])
                # Eq. 17 volume basis (sigma_i sigma_j)^(3/2) -> the Engine's ((sigma_i + sigma_j)/2)^3 at constant sigma
                volume = float(site["kappa_ab"]) * (sig_a * sig_b) ** 1.5 / ((sig_a + sig_b) / 2) ** 3
                pair_edges = [edge(x, y, site["epsilon_k_ab"], volume, prov(".association_sites")) for x, y in donors]
            elif ea and eb:
                pair_edges = [rule(ida, idb, prov(": Eqs. 17-18 combining rule")), rule(idb, ida, prov(": Eqs. 17-18 combining rule"))]
            elif sa and sb:  # no donor/acceptor pair, or an acceptor-only kappa absent: Eq. 17 gives zero cross strength
                declared = {"component_id_a": ida, "component_id_b": idb, "provenance": prov(": Eq. 17 cross strength zero")}
            pairs.append({"component_id_a": ida, "component_id_b": idb, "coefficients": [
                {"identity": f"pair/{ida}/{idb}/k_ij", "family": "k_ij", "value": value(k, "dimensionless"), "provenance": prov(note)}]})
            edges += pair_edges
            none += [declared] if declared else []
            used.update({c["identifier"]["cas"]: e for c, e in ((a, ends[0]), (b, ends[1]))})
            rows.append((source, "rehner2023_binary.json", index, f"{record['id1']['cas']}/{record['id2']['cas']}", "admitted", ""))
        except Refused as cause:
            rows.append((source, "rehner2023_binary.json", index, f"{record['id1'].get('cas')}/{record['id2'].get('cas')}", "refused", str(cause)))
    packets["feos-rehner-2023"] = (used, pairs, edges, none)
    for packet_id, (admitted, pairs, extra, none) in packets.items():
        write(packet_id, admitted, pairs, extra, none, [r for r in rows if r[0] == packet_id])
    return rows


def write(packet_id, admitted, pairs, extra, none, rows):
    meta, root = META[packet_id], ROOT / "packets" / packet_id / "1"
    components = [item[0] for item in admitted.values()]
    sites = [s for item in admitted.values() for s in item[1]]
    families = [{"family_id": f"model/{k}", "kind": k, "choice": c} for k, c in (
        ("base", "pc-saft"), ("association", "general-site" if sites else "none"), ("electrolyte", "none"), ("permittivity", "none"))]
    families += [{"family_id": "model/polar", "kind": "polar", "choice": "dipole-dipole"}] if packet_id in ("feos-esper-2023", "feos-rehner-2023") else []
    document = {"schema": "epcsaft.parameters", "schema_version": 1, "document_id": packet_id, "document_version": 1, "purpose": "external-document",
                "sources": [{k: s[k] for k in ("source_id", "citation", "use_basis")} for s in meta["sources"]],
                "domains": [{"domain_id": packet_id, "kind": "source-parametrization"}], "components": components, "pairs": pairs,
                "model_families": families, "model_coefficients": [], "correlations": [],
                "topology": {"presets": [], "sites": sites, "edges": [e for item in admitted.values() for e in item[2]] + extra,
                             "no_cross_association": none}}
    files = {"parameters/parameters.json": json.dumps(document, indent=1, ensure_ascii=False) + "\n"}
    table = io.StringIO()
    csv.writer(table, lineterminator="\n").writerows([("source_file", "source_index", "identity", "disposition", "cause")] + [r[1:] for r in rows])
    files["sources/record-accounting.csv"] = table.getvalue()
    for relative, text in files.items():
        (root / relative).parent.mkdir(parents=True, exist_ok=True)
        (root / relative).write_text(text)
    origin = "\n".join(f'\n[[records]]\nretained_path = "{p}"\nsha256 = "{sha256(root / p)}"\norigin_repository = "{REPO}"\n'
                       f'origin_commit = "{COMMIT}"\norigin_tree = "{TREE}"\norigin_path = "parameters/pcsaft/{meta["file"]}"\n'
                       f'origin_sha256 = "{HASH[meta["file"]]}"\nrecipe = "data/tools/import_feos.py"' for p in files)
    (root / "provenance").mkdir(exist_ok=True)
    (root / "provenance/source-origin.toml").write_text(f'provenance_kind = "source-record-provenance-v1"\nschema_version = 1\n'
                                                         f'transformations = {json.dumps(meta["transformations"])}\n{origin}\n')
    q = json.dumps
    toml = [f"schema_version = 2\npacket_id = {q(packet_id)}\npacket_version = \"1\""] + [f"{k} = {q(meta[k])}" for k in (
        "title", "description", "license_or_use_basis", "attribution", "tolerance_rationale")]
    toml += ["\n[[sources]]\n" + "\n".join(f"{k} = {q(v)}" for k, v in s.items()) for s in meta["sources"]]
    toml += [f"\n[[domains]]\ndomain_id = {q(packet_id)}\ndescription = {q(meta['domain'])}"]
    ids = q([s["source_id"] for s in meta["sources"]])
    for path, role, media, units in (("parameters/parameters.json", "parameters", "application/json", meta["units_basis"]),
                                     ("sources/record-accounting.csv", "source", "text/csv", "Non-numeric source indices, identities and dispositions."),
                                     ("provenance/source-origin.toml", "provenance", "application/toml", "Non-numeric repository, commit, path and hash identifiers.")):
        toml.append(f"\n[[files]]\npath = {q(path)}\nrole = {q(role)}\nmedia_type = {q(media)}\nsha256 = {q(sha256(root / path))}\n"
                    f"source_refs = {ids}\ndomain_refs = [{q(packet_id)}]\nunits_basis = {q(units)}")
    (root / "packet.toml").write_text("\n".join(toml) + "\n")
    locks = json.loads((ROOT / "packet-locks.json").read_text())
    locks["packets"] = [p for p in locks["packets"] if p["packet_id"] != packet_id] + [
        {"packet_id": packet_id, "packet_version": "1", "packet_path": f"data/packets/{packet_id}/1",
         "packet_fingerprint": packet_fingerprint(root), "materialization": "byte-for-byte"}]
    (ROOT / "packet-locks.json").write_text(json.dumps(locks, indent=2) + "\n")


if __name__ == "__main__":
    counts = {}
    for packet_id, _, _, _, disposition, _ in convert(Path(sys.argv[1])):
        counts[packet_id, disposition] = counts.get((packet_id, disposition), 0) + 1
    print(json.dumps({f"{p}/{d}": n for (p, d), n in sorted(counts.items())}, indent=1))
