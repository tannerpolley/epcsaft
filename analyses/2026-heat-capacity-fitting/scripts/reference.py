"""Verify chemicals 1.5.2 against IAPWS R6-95(2018), pp. 11,15, then retain the nine targets."""
import csv
import itertools
import json
from decimal import Decimal
from importlib.metadata import version
from pathlib import Path

from chemicals import iapws as water

SET = Path(__file__).resolve().parents[1]
ROOT = SET.parents[1]
MASS = 0.018015268  # IAPWS-95 conversion; the release keeps its specific R = 461.51805.
R = 461.51805
assert version("chemicals") == "1.5.2"
rows = []


def compare(locator, T, rho, quantity, model, reference, tolerance):
    error = abs(model - reference)
    rows.append(dict(locator=locator, T_K=T, rho_kg_m3=rho, quantity=quantity,
                     model=model, reference=reference, tolerance=tolerance, absolute_error=error,
                     passed=error <= tolerance))


def half_digit(text):
    return float(Decimal(5).scaleb(Decimal(text).as_tuple().exponent - 1))


# Table 3 cp/R expression, using only Table 6's printed derivatives.
def release_cp(derivatives):
    ideal_tt, residual_tt, residual_d, residual_dd, residual_dt = derivatives
    tau, delta = 647.096 / 500, 838.025 / 322
    return R * (-tau**2 * (ideal_tt + residual_tt)
                + (1 + delta * residual_d - delta * tau * residual_dt)**2
                / (1 + 2 * delta * residual_d + delta**2 * residual_dd))


printed = [-1.93249185, -2.23440737, -0.364366650, 0.856063701, -1.12176915]
half_digits = [5e-9, 5e-9, 5e-10, 5e-10, 5e-9]
oracle = release_cp(printed)
# Monotone dependence at this state: extrema over printed-digit interval corners bound cp rounding.
corners = [release_cp([v + sign * h for v, h, sign in zip(printed, half_digits, signs)])
           for signs in itertools.product((-1, 1), repeat=5)]
cp_tolerance = max(abs(value - oracle) for value in corners)
P = water.iapws95_P(500, 838.025)
compare("Table 6 / Table 3", 500, 838.025, "cp_J_kg_K",
        water.iapws95_properties(500, P)[5], oracle, cp_tolerance)

# Table 7, all eleven single-phase rows. Strings retain the published precision.
table7 = [
    (300, 996.5560, "0.0992418352", "4.13018112", "1501.51914", "0.393062643"),
    (300, 1005.308, "20.0022515", "4.06798347", "1534.92501", "0.387405401"),
    (300, 1188.202, "700.004704", "3.46135580", "2443.57992", "0.132609616"),
    (500, 0.4350000, "0.0999679423", "1.50817541", "548.314253", "7.94488271"),
    (500, 4.532000, "0.999938125", "1.66991025", "535.739001", "6.82502725"),
    (500, 838.0250, "10.0003858", "3.22106219", "1271.28441", "2.56690919"),
    (500, 1084.564, "700.000405", "3.07437693", "2412.00877", "2.03237509"),
    (647, 358.0000, "22.0384756", "6.18315728", "252.145078", "4.32092307"),
    (900, 0.2410000, "0.100062559", "1.75890657", "724.027147", "9.16653194"),
    (900, 52.61500, "20.0000690", "1.93510526", "698.445674", "6.59070225"),
    (900, 870.7690, "700.000006", "2.66422350", "2019.33608", "4.17223802"),
]
for T, rho, pressure, cv, speed, entropy in table7:
    exact_pressure = water.iapws95_P(T, rho)
    values = water.iapws95_properties(T, exact_pressure)
    for quantity, model, text, multiplier in (
        ("pressure_Pa", exact_pressure, pressure, 1e6),
        ("cv_J_kg_K", values[4], cv, 1e3),
        ("sound_m_s", values[6], speed, 1),
        ("entropy_J_kg_K", values[2], entropy, 1e3),
    ):
        compare("Table 7", T, rho, quantity, model, float(text) * multiplier, half_digit(text) * multiplier)

for T, pressure, rho in [(275, "0.000698451167", "999.887406"),
                          (450, "0.932203564", "890.341250")]:
    compare("Table 8", T, float(rho), "pressure_Pa", water.iapws95_Psat(T),
            float(pressure) * 1e6, half_digit(pressure) * 1e6)
    compare("Table 8", T, float(rho), "rho_kg_m3", water.iapws95_rhol_sat(T), float(rho), half_digit(rho))

nist = ROOT / "analyses/2026-shared-derivative-qualification/inputs/nist-srd69-solvent-isobars.tsv"
with nist.open() as stream:
    for row in csv.DictReader(stream, delimiter="\t"):
        if row["solvent"] != "water":
            continue
        T = float(row["T_K"])
        rho = float(row["mass_density_kg_m3"])
        compare("NIST SRD 69 retained isobar", T, rho, "rho_kg_m3",
                water.iapws95_rho(T, 1e5), rho, half_digit(row["mass_density_kg_m3"]))

with (SET / "results/reference-verification.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
assert all(row["passed"] for row in rows), [row for row in rows if not row["passed"]]

(SET / "results/reference-summary.json").write_text(json.dumps(dict(
    chemicals_version=version("chemicals"), checks=len(rows), cp_oracle_J_kg_K=oracle,
    cp_oracle_tolerance_J_kg_K=cp_tolerance, cp_oracle_relative_tolerance=cp_tolerance/oracle), indent=2) + "\n")
targets = []
for T in range(280, 441, 20):
    properties = water.iapws95_properties(T, 1e6)
    targets.append(dict(T_K=T, P_Pa=1e6, rho_mol_m3=properties[0]/MASS,
                        cp_J_mol_K=properties[5]*MASS, P_sat_Pa=water.iapws95_Psat(T)))
with (SET / "inputs/iapws95-water.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(targets[0]))
    writer.writeheader()
    writer.writerows(targets)
print(json.loads((SET / "results/reference-summary.json").read_text()))
