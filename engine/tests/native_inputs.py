"""Explicit inputs of the retained native fixtures, restated through the public package (#59).

`ROWS` holds the C++ values stored by oracle/record_native_rows.py for these same inputs.
"""
import gzip
import hashlib
import itertools
import json
import math
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from epcsaft import (
    ActiveParameter,
    GAS_CONSTANT_J_PER_MOL_K as R,
    IdealCorrelation,
    IdealDIPPR107,
    IdealInterval,
    IdealNASA7,
    IdealNASA9,
    IdealPolynomial,
    IdealShomate,
    Mixture,
    ParameterFamily,
    Parameters,
    ThermochemistryRecord,
    _core,
    regression,
)
from epcsaft.equilibrium import (
    Amounts, Free, Molalities, NeutralReference, Phase, Pinned, Problem, Reaction, ReactionLogPolynomial, ReactionReference,
    ReferenceSourceBasis, Solid,
)
from tests.packets import parameters_json

ROWS = {
    line.split("\t")[1]: float(line.split("\t")[3])
    for line in (Path(__file__).parent / "oracle" / "native-rows.tsv").read_text().splitlines()[1:]
}
def same(evidence: Any, name: str, model: float | None, reference: float | None = None,
         atol: float = 1e-12, rtol: float = 1e-12, **inputs: Any) -> None:
    """Record and assert |model - reference| <= max(atol, rtol max(|model|, |reference|)).

    Without `reference`, the stored C++ row `name` is the reference at the 1e-10 equivalence
    tolerance of CODING_STANDARDS.md: a transport fault (unit, order, field) is far larger, and
    bit identity belongs to the corpora gate. A row below 1e-9 compares relatively only.
    """
    if reference is None:
        reference = ROWS[name]
        atol, rtol = (0.0 if abs(reference) < 1e-9 else 1e-10), 1e-10
    evidence.compare(name, model=model, reference=reference, **inputs, criterion={"atol": atol, "rtol": rtol})
    assert model is not None and math.isclose(model, reference, rel_tol=rtol, abs_tol=atol), name


def packet_metadata(path: Path, components: Sequence[str]) -> dict[str, Any]:
    """Evidence identity of data/packets/<id>/<version>/parameters/parameters.json: id, version and hash."""
    return {
        "packet_id": path.parts[-4],
        "packet_version": int(path.parts[-3]),
        "parameter_file": {
            "path": str(path.relative_to(path.parents[5])),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "size_bytes": path.stat().st_size,
        },
        "selected_components": list(components),
    }


def parameters(
    species: Sequence[tuple[float, float, float, float]],  # g/mol, segments, sigma / angstrom, epsilon / K
    pairs: Sequence[tuple[str, tuple[int, ...], float]] = (),  # k_ij and association site-pair entries
    sites: tuple[Sequence[int], Sequence[int]] = ((), ()),
    charges: Sequence[int] | None = None,
    extra: Sequence[tuple[str, tuple[int, ...], float, str]] = (),
    choices: tuple[str, str, str] = ("none", "", "none"),  # electrolyte, ionic dispersion, permittivity
) -> Parameters:
    document = _core.ParameterDocument()
    document.component_ids = [f"s{i}" for i in range(len(species))]
    document.site_component, document.site_multiplicity = list(sites[0]), list(sites[1])
    document.base_choice = "pc-saft"
    document.association_choice = "general-site" if sites[0] else "none"
    document.electrolyte_choice, document.ionic_dispersion_choice, document.permittivity_choice = choices
    rows: list[tuple[str, tuple[int, ...], float, str]] = list(extra)
    for i, (mass, segments, sigma, epsilon) in enumerate(species):
        rows += [("molar_mass", (i,), mass, "gram / mole"), ("segment_count", (i,), segments, "dimensionless"),
                 ("sigma_constant", (i,), sigma, "angstrom"), ("dispersion_energy_over_k", (i,), epsilon, "kelvin")]
        rows.append(("charge_number", (i,), float(charges[i]) if charges else 0.0, "elementary-charge"))
    rows += [(field, index, value, "kelvin" if field == "association_energy_over_k" else "dimensionless")
             for field, index, value in pairs]
    entries = []
    for field, index, magnitude, unit in rows:
        item = _core.ParameterEntry()
        item.identity, item.field, item.index = f"{field}/{index}", field, list(index)
        item.magnitude, item.unit = magnitude, unit
        entries.append(item)
    document.entries = entries
    site_ids = tuple((document.component_ids[c], f"site{k}") for k, c in enumerate(sites[0]))
    return Parameters(_core.build_parameters(document), tuple(document.component_ids), site_ids)


def mea_cp() -> Any:
    """A fresh copy of the frozen manufactured MEA packet, tools/fixtures/mea-cp.json.gz."""
    return json.loads(gzip.decompress((Path(__file__).resolve().parents[2] / "tools/fixtures/mea-cp.json.gz").read_bytes()))


def decode(v: Any) -> Any:
    """Rebuild the _core objects of a frozen fixture's {"type", "fields"} and {"enum", "name"} entries."""
    if isinstance(v, list):
        return [decode(x) for x in v]
    if not isinstance(v, dict):
        return v
    if "enum" in v:
        return getattr(getattr(_core, v["enum"]), v["name"])
    cls = getattr(_core, v["type"])
    fields = {k: decode(x) for k, x in v["fields"].items()}
    try:
        obj = cls()
    except TypeError:
        return cls(**fields)
    for k, x in fields.items():
        setattr(obj, k, x)
    return obj


ETHANE = (30.07, 1.6069, 3.5206, 191.42)


def budget() -> Any:
    """The former implicit fit budget, now explicit: 100 iterations and 600 s."""
    return _core.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0)


def aqueous_electrolyte(*more: tuple[str, tuple[int, ...], float, str], c_shell: float = 0.0) -> Parameters:
    """properties_support::aqueous_electrolyte(): water, Na+ and Cl- with Debye-Huckel and Born terms."""
    extra: list[tuple[str, tuple[int, ...], float, str]] = [*more,
        ("relative_permittivity", (0,), 78.09, "dimensionless"), ("ionic_region_relative_permittivity", (), 8.0, "dimensionless"),
        ("born_c_shell", (), c_shell, "dimensionless"), ("born_c_dielectric", (), 0.0, "dimensionless"),
        ("ion_suppression_coefficient", (), 7.01, "dimensionless"),
        ("sigma_term_amplitude", (0, 0), 10.11, "angstrom"), ("sigma_term_exponent", (0, 0), -0.01775, "1 / kelvin"),
        ("sigma_term_amplitude", (1, 0), -1.417, "angstrom"), ("sigma_term_exponent", (1, 0), -0.01146, "1 / kelvin")]
    for i, (solvation, diameter, born) in enumerate(((1.5, None, None), (1.0, 2.484416, 3.445), (1.0, 2.42528, 4.1))):
        extra.append(("solvation_factor", (i,), solvation, "dimensionless"))
        if diameter is not None and born is not None:
            extra += [("packing_diameter", (i,), diameter, "angstrom"), ("debye_huckel_diameter", (i,), diameter, "angstrom"),
                      ("born_diameter", (i,), born, "angstrom")]
    pairs = [("k_ij", (0, 1), -0.3), ("k_ij", (0, 2), 0.0), ("k_ij", (1, 2), 0.8),
             ("association_energy_over_k", (0, 1), 2425.7), ("association_volume", (0, 1), 0.04509)]
    species = ((18.0153, 1.2047, 2.7927, 353.95), (22.98977, 1.0, 2.8232, 230.0), (35.4527, 1.0, 2.756, 170.0))
    return parameters(species, pairs, ((0, 0), (1, 1)), (0, 1, -1), extra,
                      ("born", "exclude-same-sign-ion-pairs", "ion-fraction-suppression"))


def associating_binary(shift: float = 0.0) -> Parameters:
    """properties_support::associating_binary(); `shift` / K moves the reciprocal site pair (0, 1) energy."""
    edges = (((0, 1), 2383.4744 + shift, 0.118488), ((0, 3), 2093.8222, 0.1540467874409226),
             ((1, 2), 2093.8222, 0.1540467874409226), ((2, 3), 1804.17, 0.2039))
    pairs: list[tuple[str, tuple[int, int], float]] = [("k_ij", (0, 1), -0.052)]
    for index, energy, volume in edges:
        pairs += [("association_energy_over_k", index, energy), ("association_volume", index, volume)]
    return parameters(((61.08, 4.5354, 2.6019, 204.0438), (18.01528, 2.1945, 2.229, 141.66)), pairs,
                      ((0, 0, 1, 1), (2, 1, 2, 2)))


def p3_forms() -> dict[str, Any]:
    """properties_support::p3_forms(): cp = 30 J/(mol K), h(300 K) = 1000 J/mol, s(300 K) = 50 J/(mol K)."""
    nasa = (1000.0 / R - 30.0 * 300.0 / R, 50.0 / R - 30.0 * math.log(300.0) / R)
    return {
        "polynomial": IdealPolynomial([30.0, 0.0], 300.0, 1000.0, 50.0),
        "shomate": IdealShomate([30.0, 0.0, 0.0, 0.0, 0.0, -9.0, 50.0 - 30.0 * math.log(0.3), 0.0], 1000.0),
        "nasa7": IdealNASA7([30.0 / R, 0.0, 0.0, 0.0, 0.0, *nasa]),
        "nasa9": IdealNASA9([0.0, 0.0, 30.0 / R, 0.0, 0.0, 0.0, 0.0, *nasa]),
        "dippr107": IdealDIPPR107([0.0, 30.0, 0.0, 0.0, 0.0], 300.0, 1000.0, 50.0),
    }


def record(*forms: Any, lower: float = 250.0, upper: float = 350.0) -> ThermochemistryRecord:
    """One closed interval per component at 1 bar, in component order."""
    return ThermochemistryRecord(1.0e5, [[IdealInterval(lower, upper, True, True, IdealCorrelation(form, 1.0e5))]
                                         for form in forms])


def surface_model(**inputs: Any) -> Mixture:
    """The C++ example's model: constant cp with anchors at 300 K and 1 bar on [250, 450] K."""
    forms = (IdealPolynomial([80.0, 0.0], 300.0, 1000.0, 50.0), IdealPolynomial([35.0, 0.0], 300.0, 2000.0, 70.0))
    return Mixture(associating_binary(), thermochemistry=record(*forms, upper=450.0), **inputs)


def bubble(T: float) -> Problem:
    """equilibrium_support::boundary(): L free, V pinned 0, P free from 1e5 Pa, equimolar amounts."""
    return Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
                   T=T, P=Free(1.0e5), feed=Amounts({"s0": 0.5, "s1": 0.5}))


def mixed_difference(
    evaluate: Callable[[float, float, tuple[float, ...]], float],
    temperature: float,
    coordinate: float,
    composition: tuple[float, ...],
    directions: Sequence[tuple[float, float, tuple[float, ...]]],
    step: float,
) -> float:
    total = 0.0
    for signs in itertools.product((-1.0, 1.0), repeat=len(directions)):
        varied_temperature = temperature + step * sum(
            sign * direction[0] for sign, direction in zip(signs, directions)
        )
        varied_coordinate = coordinate + step * sum(
            sign * direction[1] for sign, direction in zip(signs, directions)
        )
        varied_composition = tuple(
            value + step * sum(
                sign * direction[2][index]
                for sign, direction in zip(signs, directions)
            )
            for index, value in enumerate(composition)
        )
        total += math.prod(signs) * evaluate(
            varied_temperature, varied_coordinate, varied_composition
        )
    return total / (2.0 * step) ** len(directions)


# Reference-calorics fixture: identical ethane 2A <=> B; only A carries an ideal-gas record,
# cp = 52 + 0.12 (T - 300) J/(mol K) anchored at h = s = 0 at 300 K.
DIMER_LAW = ReactionLogPolynomial(math.log(0.3), 200.0, 3.0, -0.004, 300.0, 250.0, 350.0, "analytic-source")
CP_A = (52.0, 0.12)
A_RECORD = [IdealInterval(200.0, 400.0, True, True, IdealCorrelation(IdealPolynomial(CP_A, 300.0, 0.0, 0.0), 1.0e5))]


def records(*species: list[IdealInterval]) -> ThermochemistryRecord:
    return ThermochemistryRecord(1.0e5, list(species))


def twin_ethane(*species: list[IdealInterval]) -> Mixture:
    return Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]), thermochemistry=records(*species))


def dimer_problem(temperature: float, reference_pressure: float | None,
                  feed: dict[str, float] | None = None, pressure: float = 3.0e6) -> Problem:
    return Problem(
        T=temperature, P=pressure, feed=Amounts(feed or {"s0": 2.0, "s1": 0.2}), phases=[Phase("L", kind="liquid")],
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"s0": -2.0, "s1": 1.0}, correlation=DIMER_LAW, reference=ReactionReference(
            ReferenceSourceBasis.RawMoleFractionInfiniteDilution, reference_pressure))],
    )


def ionic_problem(temperature: float, pressure: float, reference_pressure: float | None) -> Problem:
    """Manufactured formula reaction H2O ⇌ Na+ + Cl- (charge neutral, not atom balanced) on
    the native equilibrium_checks aqueous case."""
    law = ReactionLogPolynomial(math.log(1e-8), 0.0, 0.0, 0.0, 298.15, 250.0, 350.0, "analytic-source")
    return Problem(
        T=temperature, P=pressure, feed=Amounts({"s0": 1.0, "s1": 1e-6, "s2": 1e-6}),
        phases=[Phase("L", kind="liquid")],
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0, "s2": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"s0": -1.0, "s1": 1.0, "s2": 1.0}, correlation=law, reference=ReactionReference(
            ReferenceSourceBasis.RawMoleFractionInfiniteDilution, reference_pressure))],
    )


WATER_ISOBUTANOL = ("water", "isobutanol")


def water_isobutanol(theta: Sequence[float] | None = None, active: bool = True) -> Mixture:
    """The engine/README.md water/isobutanol record (Khudaida 2026 packet 7). `theta` edits k_ij(isobutanol, water)
    and isobutanol eps/k, which are the active parameters unless `active` is False."""
    mapping = json.loads(parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7).read_text())
    if theta is not None:
        pair = next(p for p in mapping["pairs"] if p["coefficients"][0]["identity"] == "pair/isobutanol/water/k_ij")
        pair["coefficients"][0]["value"]["magnitude"] = theta[0]
        iso = next(c for c in mapping["components"] if c["component_id"] == "isobutanol")
        next(k for k in iso["coefficients"] if k["identity"].endswith("dispersion_energy_over_k"))["value"]["magnitude"] = theta[1]
    p = Parameters.from_mapping(mapping, components=WATER_ISOBUTANOL)
    return Mixture(p, active_parameters=[regression.active_parameter(p, "k_ij", ("isobutanol", "water")),
                   regression.active_parameter(p, "dispersion_energy_over_k", "isobutanol")] if active else None)


def water_isobutanol_problem(feed: Sequence[float] = (0.7, 0.3), T: float = 298.15) -> Problem:
    """Two liquids at 1e5 Pa from the guesses (0.46, 0.54) and (0.99, 0.01); `feed` in mol of water, isobutanol."""
    return Problem(phases=[Phase("L1", kind="liquid", support=WATER_ISOBUTANOL, composition_guess=(0.46, 0.54)),
                           Phase("L2", kind="liquid", support=WATER_ISOBUTANOL, composition_guess=(0.99, 0.01))],
                   T=T, P=1.0e5, feed=Amounts(dict(zip(WATER_ISOBUTANOL, feed, strict=True))))


# Packet 3/4 nine-species MEA/CO2/water: species order and reactions R1-R5.
MEA = ("carbon-dioxide", "monoethanolamine", "water", "protonated-monoethanolamine", "carbamate-anion",
        "bicarbonate-anion", "carbonate-anion", "hydronium-cation", "hydroxide-anion")
MEA_REACTIONS = (  # R1-R5: water, CO2 and bicarbonate dissociation, carbamate hydrolysis, MEAH+ dissociation
    {"water": -2.0, "hydronium-cation": 1.0, "hydroxide-anion": 1.0},
    {"carbon-dioxide": -1.0, "water": -2.0, "bicarbonate-anion": 1.0, "hydronium-cation": 1.0},
    {"water": -1.0, "bicarbonate-anion": -1.0, "carbonate-anion": 1.0, "hydronium-cation": 1.0},
    {"monoethanolamine": 1.0, "water": -1.0, "carbamate-anion": -1.0, "bicarbonate-anion": 1.0},
    {"protonated-monoethanolamine": -1.0, "water": -1.0, "monoethanolamine": 1.0, "hydronium-cation": 1.0})


def mea_mapping(path: Path) -> dict[str, Any]:
    """An MEA document with its historical CO2/MEA omission declared explicitly (#569): no cross edge, not a
    claim that the physical interaction is absent."""
    mapping: dict[str, Any] = json.loads(path.read_text())
    mapping["topology"].setdefault("no_cross_association", []).append({
        "component_id_a": "carbon-dioxide", "component_id_b": "monoethanolamine",
        "provenance": {"source_id": mapping["document_id"], "locator": f"{path.relative_to(Path(__file__).resolve().parents[2])}: no CO2/MEA edge"}})
    return mapping


def nominally_unloaded_mea(thermochemistry: bool = True) -> tuple[Any, list[Reaction], NeutralReference, list[float], dict[str, Any]]:
    """The Engine-owned reference declaration engine/data/reference/mea-nominally-unloaded: model (with ideal
    calorics unless thermochemistry is False), reactions, the application's neutral reference, the shipped liquid guess
    and the reaction system."""
    root = Path(__file__).resolve().parents[1] / "data" / "reference"
    request = json.loads((root / "mea-nominally-unloaded" / "request.json").read_text())
    system = request["reaction_system"]
    mapping = mea_mapping(root / "mea-nominally-unloaded" / "parameters.json")
    for family in mapping["model_families"]:
        if family["kind"] == "electrolyte":
            family.update(c_shell=1.0, c_dielectric=1.0)
    calorics = json.loads((root / "mea-zero-inventory" / "ideal-calorics.json").read_text())
    p0, components = calorics["reference_pressure_pa"], calorics["components"]
    intervals: list[list[IdealInterval]] = [[] for _ in MEA]
    for name in ("carbon-dioxide", "water"):
        lower, upper = components[name]["range_k"]
        intervals[MEA.index(name)] = [IdealInterval(lower, upper, True, True, IdealCorrelation(
            IdealShomate(tuple(components[name]["coefficients"]), components[name]["formation_enthalpy_j_per_mol"]), p0))]
    mea = components["monoethanolamine"]  # Cp(T) re-centred on the h = s = 0 anchor
    centre = mea["enthalpy_entropy_anchor_k"]
    shifted = [sum(a * math.comb(j, k) * centre ** (j - k) for j, a in enumerate(mea["coefficients"]) if j >= k)
               for k in range(4)]
    lower, upper = mea["range_k"]
    intervals[1] = [IdealInterval(lower, upper, True, True, IdealCorrelation(IdealPolynomial(shifted, centre, 0.0, 0.0), p0))]
    model = Mixture(Parameters.from_mapping(mapping), thermochemistry=ThermochemistryRecord(p0, intervals) if thermochemistry else None)
    reactions = []
    for record in system["engine_reactions"]:
        law, reference = record["engine_correlation"], record["engine_reference"]
        reactions.append(Reaction(dict(zip(MEA, record["stoichiometry"])), name=record["reaction_id"],
            correlation=ReactionLogPolynomial(*(law[key] for key in (
                "a", "b", "c", "d", "reference_temperature", "temperature_min", "temperature_max", "standard_state_id"))),
            reference=ReactionReference(getattr(ReferenceSourceBasis, reference["source_basis"]),
                                        reference_pressure_pa=reference["reference_pressure_pa"])))
    # The application's neutral reference: neutral solutes 1, source ion pairs summed, scaled by 1/8.
    ratios = dict(zip(MEA[:2] + MEA[3:], (0.125, 0.125, 0.625, 0.25, 0.125, 0.125, 0.125, 0.125)))
    guess = request["continuation"]["state"]["phases"][0]["mole_fractions"]
    reference = NeutralReference(solvent="water", solute_molality_ratios=ratios, standard_molality_mol_per_kg=1.0)
    return model, reactions, reference, guess, system


def nominally_unloaded_mea_liquid() -> tuple[Any, Problem]:
    """The reference declaration as one liquid at 353.15 K and 101325 Pa; tools/time-key-solves times it too."""
    model, reactions, reference, guess, system = nominally_unloaded_mea()
    return model, Problem(phases=[Phase("L", kind="liquid", support=MEA, composition_guess=guess)], T=353.15,
                          P=101325.0, feed=Amounts(dict(zip(MEA, system["feed_amounts_mol"]))), reactions=reactions,
                          neutral_reference=reference, topology_declared=True)


# Halite in water, Held 2014 strategy 2: Clarke and Glew 1985 (doi:10.1063/1.555730) Table 34 dissolution constants
# at 25 C (Gibbs energy, enthalpy in J/mol; heat capacity in J/(mol K)); halite molar volume 27.1 cm3/mol from
# PHREEQC pitzer.dat (phreeqpython 1.6.2).
HALITE_T, HALITE_P = 298.15, 101325.0
HALITE_SPECIES = ("water", "sodium-cation", "chloride-anion")
_DG, _DH, _DCP = -9038.7, 3858.4, -134.03
HALITE_LN_K = -_DG / (R * HALITE_T)
_B = -(_DH - HALITE_T * _DCP) / R
HALITE_LAW = ReactionLogPolynomial(HALITE_LN_K - _B / HALITE_T, _B, _DCP / R, 0.0, HALITE_T, 288.15, 308.15,
                                   "molality-infinite-dilution")
HALITE_V_SOLID = 27.1e-6
HALITE_PARAMETERS = (
    ("component/sodium-cation/dispersion_energy_over_k", ActiveParameter(ParameterFamily.DispersionEnergyOverK, (1,)), 230.0),
    ("component/chloride-anion/segment_diameter", ActiveParameter(ParameterFamily.SigmaConstant, (2,)), 2.756),
    ("pair/chloride-anion/water/k_ij", ActiveParameter(ParameterFamily.Kij, (2, 6)), -0.25))
WATER_MOLAR_MASS = 18.01528e-3  # kg/mol


def held_parameters(shift: tuple[str, float] = ("", 0.0)) -> Parameters:
    """The Held 2014 set with one document coefficient shifted by an absolute amount."""
    document = json.loads(parameters_json("held-2014-strategy2", 1).read_text())
    for record in (*document["components"], *document["pairs"]):
        for item in record["coefficients"]:
            item["value"]["magnitude"] += shift[1] if item["identity"] == shift[0] else 0.0
    return Parameters.from_mapping(document, components=list(HALITE_SPECIES))


def held(shift: tuple[str, float] = ("", 0.0), active: tuple[ActiveParameter, ...] = ()) -> Mixture:
    return Mixture(held_parameters(shift), active_parameters=list(active))


def halite_in_water(molality: float, amount: Free | Pinned = Free(), temperature: float = HALITE_T,
                   pressure: float = HALITE_P, correlation: ReactionLogPolynomial = HALITE_LAW) -> Problem:
    """One liquid beside halite, Na+ and Cl- at `molality` in 1 kg water; the solid reference pressure is HALITE_P."""
    halite = Solid("halite", formula={"sodium-cation": 1.0, "chloride-anion": 1.0}, correlation=correlation,
                   reference=ReactionReference(ReferenceSourceBasis.CommonMolalityInfiniteDilution, HALITE_P),
                   molar_volume_m3_per_mol=HALITE_V_SOLID, amount=amount)
    return Problem(T=temperature, P=pressure, phases=[Phase("L", kind="liquid")], solids=[halite],
                   feed=Molalities({"sodium-cation": molality, "chloride-anion": molality}, solvent=("water",),
                                   solvent_mass_kg=1.0),
                   neutral_reference=NeutralReference(solvent="water", standard_molality_mol_per_kg=1.0,
                                                      solute_molality_ratios={"sodium-cation": 1.0, "chloride-anion": 1.0}))
