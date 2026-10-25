"""Read retained reference diagnostics into notebook tables without running a model."""
from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parent
components = pd.read_csv(root.parents[1] / "results/ionic/numerical-components.csv").set_index("method")
terminal = pd.read_csv(root / "truth-terminal.csv")
ladder = pd.read_csv(root / "truth-ladder.csv")
wide = pd.read_csv(root / "truth-wide.csv")
quantities = [
    ("Mean reference-value contribution, input and output errors", components.loc["eiv", "reference_value"], "scaled parameter coordinate"),
    ("Preceding reference molality", ladder.iloc[2]["m"], "mol/kg"),
    ("Returned terminal molality at the 2026-10-14 freeze", ladder.iloc[8]["m"], "mol/kg"),
    ("Cation reference change", terminal.iloc[0]["old_difference"], "dimensionless ln(phi)"),
    ("Anion reference change", terminal.iloc[1]["old_difference"], "dimensionless ln(phi)"),
    ("Maximum absolute gamma precision discrepancy", wide["gamma_difference"].abs().max(), "dimensionless gamma"),
]
pd.DataFrame(quantities, columns=["Quantity", "Value", "Basis"]).to_html(
    root / "quantities.html", index=False, float_format="{:.6g}".format
)
allowances = pd.read_csv(root / "allowance-predictions.csv").set_index("option")
selected = allowances.loc[["fine_1.21e-12", "new_terminal_1e-16"], ["epsilon_conditional"]]
selected.index = ["Closer preceding rung, 1.21e-12 mol/kg", "Lower terminal, 1e-16 mol/kg"]
selected.columns = ["Conditional allowance (scaled parameter coordinate)"]
selected.to_html(root / "conditional-allowances.html", float_format="{:.6g}".format)
