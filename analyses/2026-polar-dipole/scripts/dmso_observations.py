"""Measured DMSO saturation rows from the NIST ThermoML Archive (#28, #145).

Archive: ThermoML.v2020-09-30.tgz, https://doi.org/10.18434/mds2-2422
(SHA-256 231161b5e443dc1ae0e5da8429d86a88474cb722016e5b790817bb31c58d7ec2). Extract it, then
    python analyses/2026-polar-dipole/scripts/dmso_observations.py /path/to/extracted/archive
Vapor pressures: the liquid-gas vapor-pressure tables of three papers; replicate measurements at one
temperature are averaged. Liquid density: Paknejad et al. 2019 at the lowest measured pressure of each
isotherm (100 kPa to 333 K, 500 kPa above); the compression from p_sat to that pressure is below the
1 kg/m3 reporting resolution. Uncertainties are the papers' expanded uncertainties as recorded in ThermoML.
"""
import csv
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "inputs/dmso-saturation.csv"
INCHIKEY = "IAZDPXIOMUYVGZ-UHFFFAOYSA-N"
M = 78.13  # g/mol, GV2006 Table 2 (the packet value)
SOURCES = {  # DOI: (property, citation)
    "10.1016/j.fluid.2011.02.002": ("Vapor or sublimation pressure, kPa", "Fulem, Ruzicka and Ruzicka, Fluid Phase Equilib. 303 (2011) 205, static"),
    "10.1016/j.jct.2018.04.014": ("Vapor or sublimation pressure, kPa", "de Pablo et al., J. Chem. Thermodyn. 2018, static"),
    "10.1021/je400531a": ("Vapor or sublimation pressure, kPa", "Zhang et al., J. Chem. Eng. Data 2013, ebulliometric"),
    "10.1021/acs.jced.9b00396": ("Mass density, kg/m3", "Paknejad, Mohammadkhani and Zarei, J. Chem. Eng. Data 64 (2019) 5174, vibrating tube"),
}


def points(root: ET.Element, prop: str) -> list[tuple[float, float, float, float | None]]:
    """(T/K, value, expanded uncertainty, pressure/kPa or None) for pure-DMSO liquid(-gas) rows of one property."""
    dmso = next(c.find(".//{*}nOrgNum").text for c in root.findall("{*}Compound") if c.findtext(".//{*}sStandardInChIKey") == INCHIKEY)
    out = []
    for data in root.findall("{*}PureOrMixtureData"):
        if [c.findtext(".//{*}nOrgNum") for c in data.findall("{*}Component")] != [dmso]:
            continue
        if data.findtext(".//{*}ePropName") != prop or "Crystal" in [p.text for p in data.findall(".//{*}ePhase")]:
            continue
        kinds = {v.findtext(".//{*}nVarNumber"): v.findtext(".//{*}eTemperature") or v.findtext(".//{*}ePressure") for v in data.findall("{*}Variable")}
        fixed = next((float(c.findtext("{*}nConstraintValue")) for c in data.findall("{*}Constraint") if c.findtext(".//{*}ePressure")), None)
        for row in data.findall("{*}NumValues"):
            values = {kinds[v.findtext("{*}nVarNumber")][:1]: float(v.findtext("{*}nVarValue")) for v in row.findall("{*}VariableValue")}
            value = row.find("{*}PropertyValue")
            u = value.findtext(".//{*}nCombExpandUncertValue")
            out.append((values["T"], float(value.findtext("{*}nPropValue")), float(u), values.get("P", fixed)))
    return out


def main() -> None:
    archive = Path(sys.argv[1])
    pressure, density = [], []
    for doi, (prop, citation) in SOURCES.items():
        found = points(ET.parse(archive / f"{doi}.xml").getroot(), prop)
        if prop.startswith("Vapor"):
            groups = defaultdict(list)
            for T, p, u, _ in found:
                groups[T].append((p, u))
            pressure += [(T, 1e3 * sum(p for p, _ in g) / len(g), max(u / p for p, u in g), None, doi, citation) for T, g in sorted(groups.items())]
        else:
            lowest = {}
            for T, rho, u, P in found:
                if T not in lowest or P < lowest[T][2]:
                    lowest[T] = (rho, u, P)
            density += [(T, None, None, (rho / M * 1e3, u / rho, P), doi, citation) for T, (rho, u, P) in sorted(lowest.items())]
    rows = []
    for series in (sorted(pressure), sorted(density)):
        for k, (T, p, u_p, rho, doi, citation) in enumerate(series):
            rows.append({"component_id": "dimethyl-sulfoxide", "T_K": T, "p_sat_Pa": "" if p is None else round(p, 4),
                         "rho_liq_mol_m3": "" if rho is None else repr(rho[0]), "role": "calibration" if k % 2 == 0 else "validation",
                         "reference_eos": f"{citation}; ThermoML doi:{doi}" + ("" if rho is None else f"; P = {rho[2]:g} kPa"),
                         "u_p_sat_rel": "" if u_p is None else round(u_p, 4), "u_rho_liq_rel": "" if rho is None else round(rho[1], 6)})
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(pressure)} vapor-pressure and {len(density)} density rows -> {OUT}")


if __name__ == "__main__":
    main()
