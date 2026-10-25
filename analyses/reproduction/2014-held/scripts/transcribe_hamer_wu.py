"""Extract Hamer-Wu 1972 osmotic coefficients for the Held 2014 Table 5 salts.

Reads the Zotero Markdown transcription (attachment VYE43TLN of 28NPJ5H6) given
as the only argument and writes inputs/hamer-wu-1972-osmotic.csv. Rows marked
'b' (supersaturated, calculated) are dropped; 'a' (saturated) rows are kept.
"""
import csv
import hashlib
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
# Hamer-Wu table number -> salt label used in Held 2014 Table 5.
TABLES = {4: "HCl", 5: "HBr", 6: "HI", 7: "HClO4", 10: "LiBr", 11: "LiI", 12: "LiOH", 13: "LiClO4",
          14: "LiNO3", 15: "NaF", 17: "NaBr", 18: "NaI", 19: "NaOH", 21: "NaClO4", 23: "NaNO3",
          24: "NaH2PO4", 27: "KF", 28: "KCl", 29: "KBr", 30: "KI", 31: "KOH", 34: "KNO3", 35: "KH2PO4",
          44: "CsF", 45: "CsCl", 46: "CsBr", 47: "CsI", 48: "CsOH", 49: "CsNO3", 55: "NH4Cl",
          56: "NH4ClO4", 57: "NH4NO3", 61: "NaFo", 62: "LiAc", 63: "NaAc", 64: "KAc", 66: "CsAc"}


def number(cell):
    flag = "b" if "{b}" in cell or "mathrm{b}" in cell else "a" if "{a}" in cell or "text {a}" in cell else ""
    cell = re.sub(r"\{ *\}\^\{[^}]*\{?[ab]\}?\}|\\\(|\\\)|\\text|\\mathrm", "", cell)
    digits = re.sub(r"[^0-9.]", "", cell)
    return (float(digits), flag) if digits.strip(".") else (None, "")


def main(source):
    text = Path(source).read_text()
    parts = re.split(r"\\caption\{(?:Table|TABLE) (\d+)\. [^\n]*\}", text)
    rows = []
    for number_text, body in zip(parts[1::2], parts[2::2]):
        salt = TABLES.get(int(number_text))
        if salt is None:
            continue
        tabular = body.split("\\end{tabular}")[0]
        for line in [r for r in tabular.split("\\hline") if "&" in r and "=" not in r][1:]:
            cells = [c.strip().rstrip("\\").strip() for c in line.split("&")]
            for i in range(0, len(cells) - 2, 3):
                (m, fm), (phi, fp), (gamma, _) = number(cells[i]), number(cells[i + 1]), number(cells[i + 2])
                if m is None or phi is None or "b" in (fm, fp):
                    continue
                rows.append(dict(salt=salt, hamer_wu_table=int(number_text), molality_mol_kg=m,
                                 osmotic_coefficient=phi, gamma_pm_m=gamma, saturated=fm == "a"))
    rows.sort(key=lambda r: (r["hamer_wu_table"], r["molality_mol_kg"]))
    with (HERE / "inputs/hamer-wu-1972-osmotic.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(len(rows), "rows; source sha256", hashlib.sha256(Path(source).read_bytes()).hexdigest())


if __name__ == "__main__":
    main(sys.argv[1])
