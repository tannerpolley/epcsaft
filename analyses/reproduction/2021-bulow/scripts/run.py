"""Bülow, Ascani & Held 2021 (ePC-SAFT advanced II): Engine states at the paper's conditions, 298.15 K and 1 bar.

Python declares parameter documents, observations, molalities and compositions and reads native results:
transfer energies are native TransferGibbsEnergy predictions; solvent activities and sum_i nu_i ln(x_i phi_i)
come from the C++ driver bulow2021_states.cpp linked to the Engine library. The only Python
arithmetic on model values is the labelled closed-form Born checks of Eq. (24), the permittivity mapping
check and the unit-conversion candidates. Ion pairing (Bjerrum) and the solid-liquid solve are not Engine
capabilities (notebook "Capability gaps").
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import subprocess
import time

import epcsaft
from epcsaft import Mixture, Parameters, regression as r

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[2]
T, P, R = 298.15, 1.0e5, 8.314462618
# Tolerances, fixed before the first run (README "Tolerances"): digitization-limited comparisons
# against the paper's own model values.
TOL_TRANSFER_KJ = 0.5                                   # Fig. 8b bar, +-1.5 px outline reading = +-0.21 kJ/mol
TOL_GAMMA = {"fig4-licl": 0.005, "fig4-libr": 0.010, "figS1-licl-ethanol": 0.010}  # ~3 px per panel
PIXEL_HALF_WIDTH_FIG7 = 1.5                             # Fig. 7 bar-end reading, px
# Fig. 9 cross-solvent test (#418): every bar is evaluated at these bar-end offsets (px); the outer two bound
# the reading and the inner ones check that sum_i nu_i ln(x_i phi_i) rises monotonically across it.
READING_OFFSETS_PX = (-1.5, -0.75, 0.0, 0.75, 1.5)
TOL_BORN_IDENTITY = 1e-10                               # I1, absolute in ln phi
TOL_PERMITTIVITY = 1e-12                                # I2, relative

COMPONENTS = {r["component_id"]: r for r in csv.DictReader((HERE / "inputs/parameters/components.csv").open())}
PAIRS = list(csv.DictReader((HERE / "inputs/parameters/pairs.csv").open()))
SOLVENTS = ("water", "methanol", "ethanol")
# Inorganic-ion diameter conventions as (hard-sphere packing, Debye-Hueckel and Born) factors on sigma;
# dispersion always uses sigma. Part II Section 2.4 states the Born distance is sigma; Held 2014 Eq. 21
# (cited in Section 3.1) gives d = 0.88 sigma. "hard-sphere-0.88" is the earlier packing-only probe.
DIAMETERS = {"sigma": (1.0, 1.0), "hard-sphere-0.88": (0.88, 1.0), "unified-0.88": (0.88, 0.88)}


def value(magnitude, unit):
    return {"value": {"magnitude": float(magnitude), "unit": unit}}


def document(model, diameter="sigma", water_sigma="correlation", kij_set="prediction", dispersion="exclude-same-sign-ion-pairs"):
    """One parameter document; model is 'advanced' (Born + component permittivity mixing) or 'revised-2014'."""
    advanced = model == "advanced"
    packing_factor, ion_factor = DIAMETERS[diameter]
    components, correlations, sites, edges = [], [], [], []
    for cid, c in COMPONENTS.items():
        charge = int(c["charge"])
        sigma = float(c["sigma_angstrom"])
        coefficients = [dict(identity=f"component/{cid}/{f}", family=f, **value(v, u)) for f, v, u in (
            ("segment_count", c["segment_count"], "dimensionless"),
            ("dispersion_energy_over_k", c["u_over_k_K"], "kelvin"))]
        if advanced:  # inert with c_shell = 0 (plain Born term, Eq. 24)
            coefficients.append(dict(identity=f"component/{cid}/solvation_factor", family="solvation_factor", **value(1, "dimensionless")))
        if charge:
            inorganic = c["hard_sphere_diameter"] == "sigma"
            ion_d = ion_factor * sigma if inorganic else sigma
            coefficients += [dict(identity=f"component/{cid}/{f}", family=f, **value(v, "angstrom")) for f, v in (
                ("segment_diameter", sigma), ("debye_huckel_diameter", ion_d))]
            # inorganic ions d = packing_factor sigma; IL ions the PC-SAFT d(T) (Section 3.1), declared at 298.15 K
            # because the Engine takes a fixed packing diameter for charged species (exact for these isothermal states)
            d = packing_factor * sigma if inorganic else sigma * (1 - 0.12 * math.exp(-3 * float(c["u_over_k_K"]) / T))
            coefficients.append(dict(identity=f"component/{cid}/packing_diameter", family="packing_diameter", **value(d, "angstrom")))
            if advanced:
                coefficients += [dict(identity=f"component/{cid}/born_diameter", family="born_diameter", **value(ion_d, "angstrom")),
                                 dict(identity=f"component/{cid}/relative_permittivity", family="relative_permittivity",
                                      **value(c["relative_permittivity_298_15"], "dimensionless"))]
        elif cid == "water" and water_sigma == "correlation":
            correlations.append(dict(correlation_id="component/water/segment_diameter", component_id=cid,
                family="segment_diameter", form="constant-plus-sum-of-exponentials",
                constant=dict(identity="component/water/segment_diameter/constant", **value(sigma, "angstrom")),
                terms=[dict(amplitude=dict(identity=f"component/water/segment_diameter/term-{i}/amplitude", **value(a, "angstrom")),
                            exponent_coefficient=dict(identity=f"component/water/segment_diameter/term-{i}/exponent", **value(b, "1 / kelvin")))
                       for i, (a, b) in enumerate(((10.11, -0.01775), (-1.417, -0.01146)))]))
        else:
            coefficients.append(dict(identity=f"component/{cid}/segment_diameter", family="segment_diameter", **value(sigma, "angstrom")))
        if not charge:
            correlations.append(dict(correlation_id=f"component/{cid}/relative_permittivity", component_id=cid,
                family="relative_permittivity", form="constant",
                constant=dict(identity=f"component/{cid}/relative_permittivity/constant",
                              **value(c["relative_permittivity_298_15"], "dimensionless"))))
            sites += [dict(component_id=cid, site_id=s, site_role=role, multiplicity=1) for s, role in (("a", "donor"), ("b", "acceptor"))]
            edges.append(dict(endpoint_a=dict(component_id=cid, site_id="a"), endpoint_b=dict(component_id=cid, site_id="b"),
                              energy_over_k=dict(identity=f"association/{cid}/energy", **value(c["association_energy_K"], "kelvin")),
                              volume=dict(identity=f"association/{cid}/volume", **value(c["association_volume"], "dimensionless")),
                              source=dict(kind="explicit")))
        components.append(dict(component_id=cid, fixed=dict(molar_mass=value(c["molar_mass_g_mol"], "gram / mole"),
                                                            charge_number=value(charge, "elementary-charge")),
                               coefficients=coefficients))
    charges = {cid: int(c["charge"]) for cid, c in COMPONENTS.items()}
    same_sign = lambda p: charges[p["component_a"]] * charges[p["component_b"]] > 0
    pairs = [dict(component_id_a=p["component_a"], component_id_b=p["component_b"],
                  coefficients=[dict(identity=f"pair/{p['component_a']}/{p['component_b']}/k_ij", family="k_ij",
                                     **value(p["k_ij_298_15"] if p["set"] in ("all", kij_set) else 0.0, "dimensionless"))])
             for p in PAIRS if dispersion == "all-ion-pairs" or not same_sign(p)]  # prediction: ion-alcohol k_ij zero (Section 2.2)
    # solvent-solvent k_ij: inert, since every evaluated state contains one solvent (transfer mixtures hold both)
    pairs += [dict(component_id_a=a, component_id_b=b, coefficients=[dict(identity=f"pair/{a}/{b}/k_ij", family="k_ij", **value(0.0, "dimensionless"))])
              for a, b in (("water", "methanol"), ("water", "ethanol"), ("methanol", "ethanol"))]
    electrolyte = dict(kind="electrolyte", ionic_dispersion=dispersion,
                       **(dict(choice="born", c_shell=0.0, c_dielectric=0.0) if advanced else dict(choice="fully-dissociated-debye-huckel")))
    families = [dict(kind="base", choice="pc-saft"), dict(kind="association", choice="general-site"), electrolyte,
                dict(kind="permittivity", choice="component-permittivity-mixing" if advanced else "solvent-only")]
    model_coefficients = [dict(identity="model/ionic_region_relative_permittivity", family="ionic_region_relative_permittivity",
                               **value(1.0, "dimensionless"))] if advanced else []
    return dict(schema="epcsaft.parameters", schema_version=1, components=components, correlations=correlations,
                pairs=pairs, model_families=families, model_coefficients=model_coefficients,
                topology=dict(sites=sites, edges=edges, no_cross_association=[  # the solvent-solvent k_ij = 0 reason above
                    dict(component_id_a=a, component_id_b=b, provenance=dict(locator="solvent-solvent pairs: one solvent per state"))
                    for a, b in (("water", "methanol"), ("water", "ethanol"), ("methanol", "ethanol"))]))


def mixture(doc, ids):
    return Mixture(Parameters.from_mapping(doc, components=list(ids)))


def write(name, rows):
    with (HERE / "results" / name).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def transfer(doc, solvent, cation, anion, ion):
    """Native Eq. (26): Delta_tr G_ion water -> solvent (kJ/mol), infinite dilution, mole-fraction basis."""
    p = Parameters.from_mapping(doc, components=["water", solvent, cation, anion])
    item = r.electrolyte(p, "transfer_gibbs_energy", T, P, {solvent: 1.0}, ions=(cation, anion), stoichiometry=(1, 1),
                         molality=0.0, species=ion, reference={"water": 1.0}, observed=0.0, form="difference", scale=1.0)
    prediction = r.predict(Mixture(p), item)
    assert prediction.status == r.ObservationStatus.Available, prediction.message
    return prediction.value / 1000


def transfer_table(doc):
    """{(solvent, ion): kJ/mol} for the Fig. 8b ions, and the counter-ion dependence of Na+ (NaCl against NaI)."""
    found, invariance = {}, 0.0
    for solvent in ("methanol", "ethanol"):
        for cation, anion in (("sodium-cation", "chloride-anion"), ("potassium-cation", "bromide-anion"), ("sodium-cation", "iodide-anion")):
            for ion in (cation, anion):
                v = transfer(doc, solvent, cation, anion, ion)
                if (solvent, ion) in found:
                    invariance = max(invariance, abs(v - found[(solvent, ion)]))
                found[(solvent, ion)] = v
    return found, invariance


def transfer_gibbs():
    """Fig. 8b against the paper's model bars and experimental error bars, per diameter convention."""
    digitized = {(d["solvent"], d["ion"], d["series"]): d for d in csv.DictReader((HERE / "inputs/digitized/fig8b-transfer-gibbs.csv").open())}
    ions = {"Na+": "sodium-cation", "K+": "potassium-cation", "Cl-": "chloride-anion", "Br-": "bromide-anion", "I-": "iodide-anion"}
    # Independent check (Python, closed form): Born part of Delta_tr G from Eq. (24) at infinite dilution, J/mol per 1/m
    born_prefactor = 6.02214076e23 * 1.602176634e-19 ** 2 / (4 * math.pi * 8.8541878128e-12)
    eps = {s: float(COMPONENTS[s]["relative_permittivity_298_15"]) for s in SOLVENTS}
    rows, invariance = [], 0.0
    conventions = [("sigma", "correlation"), ("hard-sphere-0.88", "correlation"), ("sigma", "constant"),
                   ("hard-sphere-0.88", "constant"), ("unified-0.88", "correlation")]
    for diameter, water_sigma in conventions:
        for kij_set in ("prediction", "correlation"):
            engine, a = transfer_table(document("advanced", diameter, water_sigma, kij_set))
            no_born, b = transfer_table(document("revised-2014", diameter, water_sigma, kij_set))
            invariance = max(invariance, a, b)
            for solvent in ("methanol", "ethanol"):
                for label, ion in ions.items():
                    if (solvent, label, kij_set) not in digitized:
                        continue
                    g, n = engine[(solvent, ion)], no_born[(solvent, ion)]
                    born_analytic = born_prefactor / (DIAMETERS[diameter][1] * float(COMPONENTS[ion]["sigma_angstrom"]) * 1e-10) \
                        * (1 / eps[solvent] - 1 / eps["water"]) / 1000
                    paper = float(digitized[(solvent, label, kij_set)]["value_kJ_mol"])
                    exp = digitized[(solvent, label, "experiment")]
                    low, high = float(exp["lower_cap_kJ_mol"]), float(exp["upper_cap_kJ_mol"])
                    rows.append(dict(diameter=diameter, water_sigma=water_sigma, kij_set=kij_set, solvent=solvent, ion=label,
                                     engine_kJ_mol=g, paper_model_kJ_mol=paper, engine_minus_paper_kJ_mol=g - paper,
                                     within_tolerance=abs(g - paper) <= TOL_TRANSFER_KJ,
                                     engine_no_born_kJ_mol=n, born_analytic_kJ_mol=born_analytic,
                                     born_engine_minus_analytic_kJ_mol=g - n - born_analytic,
                                     experiment_kJ_mol=float(exp["value_kJ_mol"]), experiment_lower_kJ_mol=low,
                                     experiment_upper_kJ_mol=high, engine_within_experimental_bar=low <= g <= high,
                                     paper_within_experimental_bar=low <= paper <= high))
    write("fig8b-transfer-gibbs.csv", rows)
    return invariance


def unit_candidates():
    """Unit and standard-state conversions water -> alcohol (kJ/mol) that a Fig. 8b convention could add.

    The pure-solvent liquid densities are native; the logarithms are unit conversions, not model values."""
    rho = {s: mixture(document("advanced"), (s, "sodium-cation", "chloride-anion")).state(T, P=P, x=[1.0, 0.0, 0.0], phase="liquid").molar_density
           for s in SOLVENTS}
    M = {s: float(COMPONENTS[s]["molar_mass_g_mol"]) for s in SOLVENTS}
    rt = R * T / 1000
    rows = []
    for s in ("methanol", "ethanol"):
        for name, v in (("RT ln(M_water/M_alcohol)", rt * math.log(M["water"] / M[s])),
                        ("RT ln(mass density water/alcohol)", rt * math.log(rho["water"] * M["water"] / (rho[s] * M[s]))),
                        ("RT ln(molar density water/alcohol)", rt * math.log(rho["water"] / rho[s])),
                        ("RT ln 1000", rt * math.log(1000)), ("RT ln 1e5", rt * math.log(1e5))):
            rows.append(dict(solvent=s, candidate=name, value_kJ_mol=v, molar_density_water_mol_m3=rho["water"],
                             molar_density_alcohol_mol_m3=rho[s]))
    write("fig8b-unit-candidates.csv", rows)


def printed_slope_probe():
    """Table 3 as printed (slopes -7.981e-4, -4.012e-4) against the evaluated Held 2014 slopes (1e-3)."""
    printed = {"sodium-cation": -7.981e-4 * T + 2.38, "potassium-cation": -4.012e-4 * T + 1.396}
    doc = document("advanced")
    for pair in doc["pairs"]:
        ids = {pair["component_id_a"], pair["component_id_b"]}
        for cation, k in printed.items():
            if ids == {cation, "water"}:
                pair["coefficients"][0]["value"]["magnitude"] = k
    rows = []
    for cation, anion in (("sodium-cation", "chloride-anion"), ("potassium-cation", "bromide-anion")):
        for solvent in ("methanol", "ethanol"):
            evaluated = transfer(document("advanced"), solvent, cation, anion, cation)
            literal = transfer(doc, solvent, cation, anion, cation)
            rows.append(dict(ion=cation, solvent=solvent, k_ij_printed=printed[cation], transfer_evaluated_kJ_mol=evaluated,
                             transfer_printed_kJ_mol=literal, printed_minus_evaluated_kJ_mol=literal - evaluated))
    write("printed-slope-probe.csv", rows)


def molality(x_salt, basis, solvent):
    """Paper salt axis -> molality (mol/kg) of the 1:1 salt; the driver forms x_ion = mM/(1 + 2mM)."""
    M = float(COMPONENTS[solvent]["molar_mass_g_mol"]) / 1000
    if basis == "formula":                    # x_salt = n_salt / (n_salt + n_solvent)
        return x_salt / ((1 - x_salt) * M)
    return x_salt / ((1 - 2 * x_salt) * M)    # "species": x_salt = n_cation / (n_solvent + 2 n_salt)


def native_states(blocks):
    """Compile and run bulow2021_states.cpp on [(label, doc, ids, [("activity", molality) | ("sum", x)])]; its CSV rows."""
    job = ROOT / "build/native-drivers/bulow2021_states/jobs.txt"
    job.parent.mkdir(parents=True, exist_ok=True)
    q = json.dumps
    with job.open("w") as out:
        out.write(f"{T!r} {P!r} {len(blocks)}\n")
        for label, doc, ids, states in blocks:
            d = Parameters.from_mapping(doc, components=list(ids))._document
            out.write(f"{q(label)}\n{len(d.component_ids)} {len(d.site_component)} {len(d.entries)}\n")
            out.write(" ".join(map(q, d.component_ids)) + "\n")
            out.write("".join(f"{comp} {mult}\n" for comp, mult in zip(d.site_component, d.site_multiplicity)))
            out.write(" ".join(q(getattr(d, n)) for n in ("base_choice", "association_choice", "electrolyte_choice",
                                                          "ionic_dispersion_choice", "permittivity_choice", "polar_choice")) + "\n")
            for e in d.entries:
                out.write(f"{q(e.identity)} {q(e.field)} {len(e.index)} {' '.join(map(str, e.index))} {e.magnitude!r} {q(e.unit)}\n")
            out.write(f"{len(states)}\n" + "".join(f"{kind} {' '.join(map(repr, v)) if kind == 'sum' else repr(v)}\n"
                                                      for kind, v in states))
    driver = subprocess.check_output([str(ROOT / "tools/native-driver"), "regression_checks",
                                      str(Path(__file__).with_name("bulow2021_states.cpp"))], text=True).strip()
    return list(csv.DictReader(subprocess.check_output([driver, str(job)], text=True).splitlines()))


PANELS_4 = {"fig4-licl": ("methanol", "lithium-cation", "chloride-anion"),
            "fig4-libr": ("methanol", "lithium-cation", "bromide-anion"),
            "figS1-licl-ethanol": ("ethanol", "lithium-cation", "chloride-anion")}
SALTS_7 = {"NaCl": ("sodium-cation", "chloride-anion"), "NaBr": ("sodium-cation", "bromide-anion"),
           "KCl": ("potassium-cation", "chloride-anion"), "KBr": ("potassium-cation", "bromide-anion"),
           "KI": ("potassium-cation", "iodide-anion")}


def alcohol_x(x_salt):
    """Formula-basis salt mole fraction in one solvent -> (solvent, cation, anion) species mole fractions."""
    x_ion = x_salt / (1 + x_salt)
    return [1 - 2 * x_ion, x_ion, x_ion]


def finite_concentration():
    """Fig. 4 / S1 solvent activity coefficients and Fig. 7 sum_i nu_i ln(x_i phi_i), in one native driver run."""
    curves = list(csv.DictReader((HERE / "inputs/digitized/fig4-curves.csv").open()))
    markers = list(csv.DictReader((HERE / "inputs/digitized/fig4-markers.csv").open()))
    blocks, meta4, meta7 = [], [], []
    for diameter, kij_set in (("sigma", "prediction"), ("hard-sphere-0.88", "prediction"), ("unified-0.88", "prediction"), ("sigma", "correlation")):
        for model in ("revised-2014", "advanced"):
            doc = document(model, diameter, kij_set=kij_set)
            for panel, ids in PANELS_4.items():
                points = [("paper-model", float(c["x"]), float(c["y"])) for c in curves if c["panel"] == panel and c["series"] == model]
                points += [("experiment", float(m["x"]), float(m["y"])) for m in markers if m["panel"] == panel]
                states = []
                for basis in ("formula", "species"):
                    for kind, x_salt, y in points:
                        states.append(("activity", molality(x_salt, basis, ids[0])))
                        meta4.append(dict(diameter=diameter, kij_set=kij_set, model=model, panel=panel, basis=basis, kind=kind,
                                          x_salt=x_salt, digitized_gamma=y))
                blocks.append((f"fig4/{diameter}/{kij_set}/{model}/{panel}", doc, ids, states))
    bars = {(b["solvent"], b["salt"], b["series"]): b for b in csv.DictReader((HERE / "inputs/digitized/fig7-solubility.csv").open())}
    for diameter in ("sigma", "unified-0.88"):
        doc = document("advanced", diameter)
        for (solvent, salt, series), bar in bars.items():
            if series == "advanced-bjerrum" or not bar["x_salt"]:
                continue
            x, h = float(bar["x_salt"]), float(bar["height_px"])
            # a +-1.5 px bar-end reading moves x by the same fraction of the bar height
            xs = (x, x * (1 - PIXEL_HALF_WIDTH_FIG7 / h), x * (1 + PIXEL_HALF_WIDTH_FIG7 / h))
            blocks.append((f"fig7/{diameter}/{solvent}/{salt}/{series}", doc, (solvent, *SALTS_7[salt]),
                           [("sum", alcohol_x(v)) for v in xs]))
            meta7.append(dict(diameter=diameter, solvent=solvent, salt=salt, series=series, x_salt=x, height_px=h))
    out = native_states(blocks)
    activity = [row for row in out if row["kind"] == "activity"]
    sums = [float(row["sum_log_x_phi"]) for row in out if row["kind"] == "sum"]
    assert len(activity) == len(meta4) and len(sums) == 3 * len(meta7)
    rows4 = [dict(**m, molality_mol_kg=float(a["molality"]), osmotic_coefficient=float(a["osmotic_coefficient"]),
                  gamma_species_basis=float(a["gamma_solvent_species_basis"]), gamma_formula_basis=float(a["gamma_solvent_formula_basis"]),
                  log_solvent_activity=float(a["log_solvent_activity"]),
                  osmotic_minus_fugacity_ratio=float(a["log_solvent_activity"]) - float(a["log_solvent_activity_fugacity_ratio"]))
             for m, a in zip(meta4, activity)]
    write("fig4-solvent-activity.csv", rows4)
    rows7 = [dict(**m, sum_log_x_phi=sums[3 * i], sum_at_low_reading=sums[3 * i + 1], sum_at_high_reading=sums[3 * i + 2])
             for i, m in enumerate(meta7)]
    write("fig7-sum-log-x-phi.csv", rows7)
    pairs = []
    for diameter in ("sigma", "unified-0.88"):
        for salt in SALTS_7:
            for series in ("advanced", "experiment"):
                got = {g["solvent"]: g for g in rows7 if g["diameter"] == diameter and g["salt"] == salt and g["series"] == series}
                if len(got) < 2:
                    continue
                m, e = got["methanol"], got["ethanol"]
                # sum_i nu_i ln(x_i phi_i) rises with x here, so the reading interval maps to these endpoints
                low, high = m["sum_at_low_reading"] - e["sum_at_high_reading"], m["sum_at_high_reading"] - e["sum_at_low_reading"]
                pairs.append(dict(diameter=diameter, salt=salt, series=series, x_methanol=m["x_salt"], x_ethanol=e["x_salt"],
                                  methanol_minus_ethanol=m["sum_log_x_phi"] - e["sum_log_x_phi"],
                                  interval_low=low, interval_high=high, interval_contains_zero=low <= 0 <= high))
    write("fig7-cross-solvent.csv", pairs)
    assert all(g["sum_at_low_reading"] < g["sum_log_x_phi"] < g["sum_at_high_reading"] for g in rows7)
    return max(abs(row["osmotic_minus_fugacity_ratio"]) for row in rows4)


IL = ("c4mim-cation", "triflate-anion")
# R1 and the closed diagnostic set R2 (a)-(d) (#418 design section 4), plus the pair-rule check. One document
# serves both solvent sides. The paper's pair rule (owner choice A: inorganic like-charged pairs excluded,
# inorganic-IL like-charged and IL self pairs included) is not expressible: k_ij cannot address the Na+-Na+
# self pair, so R1 uses all-ion-pairs on both sides and "pair-rule-exclude-same-sign" brackets the effect.
VARIANTS = ("R1", "a-il-ions-out-of-born", "b-inorganic-il-like-charged-off", "c-il-ion-diameter-0.88",
            "d-species-basis", "pair-rule-exclude-same-sign")


def il_document(diameter, variant):
    doc = document("advanced", diameter, dispersion="exclude-same-sign-ion-pairs" if variant.startswith("pair-rule") else "all-ion-pairs")
    charge = {cid: int(c["charge"]) for cid, c in COMPONENTS.items()}
    for pair in doc["pairs"]:
        a, b = pair["component_id_a"], pair["component_id_b"]
        if variant.startswith("b-") and (a in IL) != (b in IL) and charge[a] * charge[b] > 0:
            pair["coefficients"][0]["value"]["magnitude"] = 1.0       # epsilon_ij = 0, the pair's dispersion off
    for comp in doc["components"]:
        for c in comp["coefficients"] if comp["component_id"] in IL else ():
            if variant.startswith("a-") and c["family"] == "born_diameter":
                c["value"]["magnitude"] *= 1e6                       # limiting case: IL ions out of the Born sum
            if variant.startswith("c-") and c["family"] in ("born_diameter", "debye_huckel_diameter"):
                c["value"]["magnitude"] *= 0.88
    return doc


def il_x(x_salt, basis):
    """Salt mole fraction in [C4mim][TfO] -> ([C4mim]+, [TfO]-, cation, anion) species mole fractions.
    formula: x_salt = n_salt/(n_salt + n_IL), so x_ion = x_salt/2; species: x_ion = x_salt."""
    x_ion = x_salt / 2 if basis == "formula" else x_salt
    return [0.5 - x_ion, 0.5 - x_ion, x_ion, x_ion]


def ionic_liquid():
    """Fig. 9 (#418): sum_i nu_i ln(x_i phi_i) in [C4mim][TfO] at the Fig. 9 bars minus the same sum in the
    alcohol at the Fig. 7 bars; zero when one K_SP serves both solvents (Eqs. 5, 28). Identities I1, I2."""
    bars = {}
    for name in ("fig7-solubility.csv", "fig9-solubility-ionic-liquid.csv"):
        for b in csv.DictReader((HERE / "inputs/digitized" / name).open()):
            if b["salt"] in ("NaCl", "NaBr", "KCl", "KBr") and b["series"] in ("advanced", "experiment") and b["x_salt"]:
                bars[(b["solvent"], b["salt"], b["series"])] = (float(b["x_salt"]), float(b["height_px"]))
    blocks, meta = [], []
    for diameter in ("sigma", "unified-0.88"):
        for variant in VARIANTS:
            doc = il_document(diameter, variant)
            basis = "species" if variant.startswith("d-") else "formula"
            for salt, ions in SALTS_7.items():
                for solvent in ("C4mim-TfO", "methanol", "ethanol"):
                    ids = (*IL, *ions) if solvent == "C4mim-TfO" else (solvent, *ions)
                    states = []
                    for series in ("advanced", "experiment"):
                        if (solvent, salt, series) not in bars:
                            continue
                        x, h = bars[(solvent, salt, series)]
                        for offset in READING_OFFSETS_PX:
                            v = x * (1 + offset / h)      # a bar-end offset moves x by the same fraction of the bar height
                            states.append(("sum", il_x(v, basis) if solvent == "C4mim-TfO" else alcohol_x(v)))
                            meta.append(dict(diameter=diameter, variant=variant, solvent=solvent, salt=salt, series=series,
                                             basis=basis if solvent == "C4mim-TfO" else "formula", offset_px=offset,
                                             x_salt=v, mole_fractions=states[-1][1], ids=ids))
                    if states:
                        blocks.append((f"fig9/{diameter}/{variant}/{solvent}/{salt}", doc, ids, states))
    out = native_states(blocks)
    assert len(out) == len(meta)
    rows = []
    for m, o in zip(meta, out):
        x = m.pop("mole_fractions")
        eps_mixed = sum(xi * float(COMPONENTS[c]["relative_permittivity_298_15"]) for xi, c in zip(x, m.pop("ids")))
        eps = float(o["bulk_relative_permittivity"])
        rows.append(dict(**m, log_phi_cation=float(o["log_phi_cation"]), log_phi_anion=float(o["log_phi_anion"]),
                         sum_log_x_phi=float(o["sum_log_x_phi"]), bulk_relative_permittivity=eps,
                         permittivity_mapping_relative_error=eps / eps_mixed - 1, molar_density_mol_m3=float(o["molar_density"]),
                         pressure_residual_relative=float(o["pressure_residual_relative"])))
    write("fig9-sum-log-x-phi.csv", rows)
    by = {(r["diameter"], r["variant"], r["solvent"], r["salt"], r["series"], r["offset_px"]): r for r in rows}
    # the interval endpoint rule needs the sum to rise across each bar's reading
    for key in {k[:5] for k in by}:
        sums = [by[(*key, o)]["sum_log_x_phi"] for o in READING_OFFSETS_PX]
        assert all(a < b for a, b in zip(sums, sums[1:])), key
    lo, mid, hi = READING_OFFSETS_PX[0], 0.0, READING_OFFSETS_PX[-1]
    pairs = []
    for (diameter, variant, solvent, salt, series, offset), il in by.items():
        if solvent != "C4mim-TfO" or offset != mid:
            continue
        for alcohol in ("methanol", "ethanol"):
            if (diameter, variant, alcohol, salt, series, mid) not in by:
                continue
            a = lambda o: by[(diameter, variant, alcohol, salt, series, o)]
            i = lambda o: by[(diameter, variant, solvent, salt, series, o)]
            d = il["sum_log_x_phi"] - a(mid)["sum_log_x_phi"]
            low = i(lo)["sum_log_x_phi"] - a(hi)["sum_log_x_phi"]
            high = i(hi)["sum_log_x_phi"] - a(lo)["sum_log_x_phi"]
            pairs.append(dict(diameter=diameter, variant=variant, series=series, salt=salt, alcohol=alcohol,
                              x_ionic_liquid=il["x_salt"], x_alcohol=a(mid)["x_salt"], il_minus_alcohol=d,
                              interval_low=low, interval_high=high, interval_contains_zero=low <= 0 <= high,
                              first_order_solubility_factor=math.exp(-d / 2)))
    write("fig9-cross-solvent.csv", pairs)
    # I1: Born is density independent, so removing the IL ions from the Born sum (a') changes each salt ion's ln phi by
    # -L0 (eps_i - eps)/eps^2 sum_{j in IL} x_j z_j^2 (1/a_j - 1/a'_j), L0 = e^2/(4 pi eps0 k_B T), at equal density.
    L0 = 1.602176634e-19 ** 2 / (4 * math.pi * 8.8541878128e-12 * 1.380649e-23 * T) * 1e10   # angstrom
    identity = []
    for (diameter, variant, solvent, salt, series, offset), r in by.items():
        if variant != "R1" or solvent != "C4mim-TfO":
            continue
        a = by[(diameter, "a-il-ions-out-of-born", solvent, salt, series, offset)]
        x_ion = r["x_salt"] / 2
        il_sum = sum((0.5 - x_ion) * (1 / float(COMPONENTS[j]["sigma_angstrom"]) - 1 / (1e6 * float(COMPONENTS[j]["sigma_angstrom"])))
                     for j in IL)
        eps = r["bulk_relative_permittivity"]
        closed = -L0 * (8.0 - eps) / eps ** 2 * il_sum            # both salt ions have eps_i = 8
        for ion in ("cation", "anion"):
            engine = r[f"log_phi_{ion}"] - a[f"log_phi_{ion}"]
            identity.append(dict(diameter=diameter, salt=salt, series=series, offset_px=offset, ion=ion,
                                 engine_log_phi_R1_minus_a=engine, closed_form=closed, difference=engine - closed,
                                 molar_density_relative_difference=r["molar_density_mol_m3"] / a["molar_density_mol_m3"] - 1))
    write("fig9-born-identity.csv", identity)
    born = max(abs(r["difference"]) for r in identity)
    density = max(abs(r["molar_density_relative_difference"]) for r in identity)
    permittivity = max(abs(r["permittivity_mapping_relative_error"]) for r in rows)
    closure = max(abs(r["pressure_residual_relative"]) for r in rows)
    assert born <= TOL_BORN_IDENTITY and density <= 1e-9 and permittivity <= TOL_PERMITTIVITY, (born, density, permittivity)
    return dict(fig9_born_identity_tolerance_abs=TOL_BORN_IDENTITY, fig9_permittivity_tolerance_relative=TOL_PERMITTIVITY,
                fig9_born_identity_max_abs=born, fig9_born_identity_density_max_relative=density,
                fig9_permittivity_mapping_max_relative=permittivity, fig9_pressure_residual_max_relative=closure)


def declarations():
    return [
        mixture(document(model, diameter, water, kij), (s, *ions))
        for model in ('advanced', 'revised-2014')
        for diameter in DIAMETERS
        for water in ('correlation', 'constant')
        for kij in ('prediction', 'correlation')
        for s in SOLVENTS
        for ions in (
            ('lithium-cation', 'chloride-anion'),
            ('lithium-cation', 'bromide-anion'),
            ('sodium-cation', 'chloride-anion'),
            ('sodium-cation', 'bromide-anion'),
            ('sodium-cation', 'iodide-anion'),
            ('potassium-cation', 'chloride-anion'),
            ('potassium-cation', 'bromide-anion'),
            ('potassium-cation', 'iodide-anion'),
        )
    ]


if __name__ == "__main__":
    start = time.perf_counter()
    invariance = transfer_gibbs()
    unit_candidates()
    printed_slope_probe()
    identity_gap = finite_concentration()
    il_checks = ionic_liquid()
    inputs = sorted(p for p in (HERE / "inputs").rglob("*") if p.is_file() and p.suffix in (".csv", ".json"))
    (HERE / "results/run-identity.json").write_text(json.dumps({
        "engine_commit": epcsaft.provenance()[0],
        "engine_tree_clean": epcsaft.provenance()[1],
        "script_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (Path(__file__), Path(__file__).with_name("bulow2021_states.cpp"))},
        "input_sha256": {p.relative_to(HERE).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        "counter_ion_invariance_max_abs_kJ_mol": invariance,
        "solvent_activity_osmotic_minus_fugacity_ratio_max_abs": identity_gap,
        "salt_mole_fraction_basis_fig7": "formula",
        **il_checks,
        "conditions": "298.15 K, 1e5 Pa, liquid branch",
        "wall_seconds": time.perf_counter() - start,
    }, indent=2) + "\n")
    print("done", time.perf_counter() - start, "s; counter-ion invariance", invariance, "kJ/mol; osmotic identity", identity_gap)
