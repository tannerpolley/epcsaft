"""Read source declarations, run native Figure 2 and retain its values."""

import csv, hashlib, json, sys
from pathlib import Path
import epcsaft as ep

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
PACKET = ROOT / "data/packets/khudaida-2026-figure-2-electrolyte-lle/7"
SOURCE = LEAF / "inputs/source"
COMPONENTS = ("water", "ethanol", "isobutanol", "nacl")


def declarations():
    names = "water ethanol isobutanol sodium-cation chloride-anion".split()
    parameters = ep.Parameters.from_json(PACKET / "parameters/parameters.json", components=names)
    return ep.Mixture(parameters)


def main():
    import figure2

    model = declarations()
    rows = lambda p: list(csv.DictReader(p.open()))
    source = rows(LEAF / "figures/figure-02/source/source-points.csv")
    pairing = rows(SOURCE / "figure-02-feed-to-table-3-tie-lines.csv")
    feeds = [r for r in source if r["dataset"] == "feed"]
    assert len(feeds) == 8
    assert [r["experimental_tie_line_id"] for r in pairing] == list(map(str, range(8, 0, -1)))
    diagnostics = []
    fields = "success message solver_status amounts mole_fractions molar_densities coordinates residuals classified_residuals iterations evaluations requested_tolerance_met".split()
    for feed, mapping in zip(feeds, pairing):
        assert feed["point_id"] == mapping["feed_point_id"]
        states = figure2.run(
            model._model, [float(feed["x_" + k]) for k in COMPONENTS], float(feed["temperature_K"]), 100000
        )
        for state in states:
            solved = state["solved"]
            state["solved"] = {k: getattr(solved, k) for k in fields}
            state["solved"]["rows_kind_phase_index_eqid"] = [
                (r.kind.name, r.phase, r.index, r.eqid) for r in solved.rows
            ]
        diagnostics.append(dict(source=feed, pairing=mapping, states=states))
        print(feed["point_id"], states[-1]["admitted"], flush=True)
    output = LEAF / "results"
    (output / "figure-2-diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    library = ROOT / "build/engine-meson/shared/libepcsaft_native.so"
    driver = Path(figure2.__file__)
    files = list(PACKET.rglob("*")) + list(SOURCE.rglob("*"))
    files += [LEAF / "figures/figure-02/source/source-points.csv", LEAF / "scripts/figure2.cpp", library, driver]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.is_file()}
    identity = dict(engine=list(ep.provenance()), packet=str(PACKET.relative_to(ROOT)), sha256=hashes)
    (output / "figure-2-run-identity.json").write_text(json.dumps(identity, indent=2) + "\n")


if __name__ == "__main__":
    main()
