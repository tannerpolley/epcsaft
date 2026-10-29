"""Join ThermoML groups 13/14 by gas-water fraction; preserve compiler uncertainty."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    source = ROOT / "inputs/bernatova-2006.json"
    data = json.loads(source.read_text())
    groups = {g["nPureOrMixtureDataNumber"]: g for g in data["PureOrMixtureData"]}
    pressure, composition = groups[13], groups[14]
    for group in (pressure, composition):
        assert group["Variable"][0]["VariableID"]["RegNum"]["nOrgNum"] == 1
        assert group["Variable"][0]["VarPhaseID"]["eVarPhase"] == "Gas"
        assert group["Constraint"][0]["nConstraintValue"] == 323.15
        assert group["Property"][0]["CombinedUncertainty"]["nCombUncertLevOfConfid"] == 95
    by_y = {r["VariableValue"][0]["nVarValue"]: r["PropertyValue"][0] for r in pressure["NumValues"]}
    assert len(by_y) == len(pressure["NumValues"]), "Nonunique pressure join key"
    rows = []
    for index, row in enumerate(composition["NumValues"], 1):
        y = row["VariableValue"][0]["nVarValue"]
        x, p = row["PropertyValue"][0], by_y[y]
        selected = index in (9, 10, 11)
        rows.append(dict(row=index, T_K=323.15, x_water=x["nPropValue"], y_water_reported=y,
                         P_Pa=1000 * p["nPropValue"],
                         U95_P_Pa=1000 * p["CombinedUncertainty"]["nCombExpandUncertValue"],
                         U95_x=x["CombinedUncertainty"]["nCombExpandUncertValue"], selected=selected,
                         exclusion="" if selected else "Acid-rich state outside adopted comparison"))
    assert len(rows) == 11 and sum(r["selected"] for r in rows) == 3
    target = ROOT / "inputs/measurements.csv"
    with target.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lineage = dict(doi="10.1016/j.fluid.2006.06.005",
                   source_url="https://trc.nist.gov/ThermoML/10.1016/j.fluid.2006.06.005.json",
                   source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                   groups=[13, 14], join_key="Gas-water mole fraction, compound 1",
                   pressure_conversion="kPa to Pa, multiply pressure and its uncertainty by 1000",
                   uncertainty_authority="Data-file compiler; combined expanded uncertainty, 95 percent",
                   uncertainty_method="Propagation of evaluated standard uncertainties",
                   temperature_uncertainty=None, covariance=None,
                   calibration_label="Fixed-parameter comparison; historical calibration overlap unresolved",
                   output_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    (ROOT / "inputs/measurement-lineage.json").write_text(json.dumps(lineage, indent=2) + "\n")
    print("Joined 11 records; retained selected and excluded rows.")


if __name__ == "__main__":
    main()
