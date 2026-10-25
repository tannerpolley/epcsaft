"""Source consistency (design section 7): Dinane-Mounir water activity against the value Rard-Miller's isopiestic Phi implies
at the same NaCl-MgCl2 composition. Source arithmetic only (no model): Phi is interpolated linearly in ionic strength along
each Rard-Miller series, then linearly in the NaCl ionic-strength fraction, and a_w = exp(-M_w q Phi). Rows come from the
frozen #556 compare.realize; evaluate.py applies the same conversion in every frozen realization (owner amendment,
2026-10-16)."""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit  # noqa: E402

RARD_MILLER, DINANE_MOUNIR = "rard-miller-1987", "dinane-mounir-2003"


def along(points: list[tuple[float, float]], ionic: float) -> float:
    i = max(1, min(len(points) - 1, next((k for k, (x, _) in enumerate(points) if x >= ionic), len(points) - 1)))
    (x0, y0), (x1, y1) = points[i - 1], points[i]
    return y0 + (ionic - x0) * (y1 - y0) / (x1 - x0)


def implied(rard_miller: list[dict], dinane_mounir: list[dict], buckets: dict[str, str], mw: float) -> list[dict]:
    """Per Dinane-Mounir row of one realization: measured a_w minus the a_w Rard-Miller's Phi implies; buckets maps the
    row to its nominal ionic-strength bucket, as #556 does."""
    series = defaultdict(list)
    for r in rard_miller:
        if r["role"] == "mixture" and r["included"]:
            series[float(r["fraction_label"])].append((r["I_mol_kg"], r["observed"]))
    fractions = sorted(series)
    for y in fractions:
        series[y].sort()
    out = []
    for r in dinane_mounir:
        m_na, m_mg = r["m_NaCl"], r["m_divalent"]
        ionic, q, y = m_na + 3 * m_mg, 2 * m_na + 3 * m_mg, m_na / (m_na + 3 * m_mg)
        hi = next(f for f in fractions if f >= y)
        lo = fractions[fractions.index(hi) - 1]
        phi = along(series[lo], ionic) + (y - lo) * (along(series[hi], ionic) - along(series[lo], ionic)) / (hi - lo)
        aw = math.exp(-mw * q * phi)
        out.append({"row_id": r["row_id"], "I_bucket": buckets[r["row_id"]], "y_NaCl": y, "I_mol_kg": ionic,
                    "aw_dinane_mounir": r["observed"], "phi_rard_miller_interpolated": phi,
                    "aw_implied_by_rard_miller": aw, "difference": r["observed"] - aw})
    return out


def nominal_buckets(nominal_dinane_mounir: list[dict]) -> dict[str, str]:
    return {r["row_id"]: "I<=3" if r["I_mol_kg"] <= 3 else "I>3" for r in nominal_dinane_mounir}


def main() -> None:
    compare = fit.comparison()
    nominal = compare.read(compare.INPUTS / "frozen-realizations.csv")[0]
    rows = {source: compare.realize(source, nominal) for source in (RARD_MILLER, DINANE_MOUNIR)}
    out = implied(rows[RARD_MILLER], rows[DINANE_MOUNIR], nominal_buckets(rows[DINANE_MOUNIR]), compare.MW_SOURCE)
    with (fit.HERE / "results/source-consistency.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(out[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)
    for bucket in ("I<=3", "I>3"):
        d = [o["difference"] for o in out if o["I_bucket"] == bucket]
        print(bucket, len(d), f"bias {sum(d) / len(d):+.5f} RMS {math.sqrt(sum(x * x for x in d) / len(d)):.5f} max {max(map(abs, d)):.5f}")


if __name__ == "__main__":
    main()
