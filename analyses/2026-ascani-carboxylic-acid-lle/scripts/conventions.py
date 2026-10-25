"""Convention probes quoted in the notebook: pure-water density by water sigma form, and the aqueous-only sweep
without water-acetic acid cross-association. Native Engine calls only; writes results/conventions.json and
results/sweep-aq-no-cross.csv."""
import json

import epcsaft

import parameters
import run

water = {}
for variant in ("literal", "constant-sigma"):
    state = run.model(variant, run.S[:2]).state(run.T, P=run.P, x=[1.0, 0.0], phase="liquid")
    water[variant] = state.molar_density * run.M_W  # kg/m3
document = parameters.document("literal")
dropped = [e for e in document["topology"]["edges"] if e["source"]["kind"] != "explicit"]
document["topology"]["edges"] = [e for e in document["topology"]["edges"] if e["source"]["kind"] == "explicit"]
# The Engine admits a dropped cross pair only when it is declared (#569); this probe declares each one.
document["topology"]["no_cross_association"] = [
    {"component_id_a": a, "component_id_b": b, "provenance": e["source"]["provenance"][0]}
    for (a, b), e in {tuple(sorted((e["endpoint_a"]["component_id"], e["endpoint_b"]["component_id"]))): e
                      for e in dropped}.items()]
run.model = lambda variant, species=run.S: epcsaft.Mixture(epcsaft.Parameters.from_mapping(document, components=species))
rows = run.sweep("aq-no-cross", "no-cross", False, 161)
x_acid = [r["x_acetic-acid_aq"] for r in rows if r["success"]]
(run.RESULTS / "conventions.json").write_text(json.dumps(dict(
    water_density_kg_m3=water, no_cross_solves=len(rows), no_cross_success=len(x_acid),
    no_cross_x_acid_molecular_aq_min=min(x_acid), no_cross_x_acid_molecular_aq_max=max(x_acid)), indent=1) + "\n")
