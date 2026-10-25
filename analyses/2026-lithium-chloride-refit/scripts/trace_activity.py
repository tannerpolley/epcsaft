"""Trace log10 gamma+-(LiCl) in chloride brines, read from native states (Lithium's #213 route)."""

from __future__ import annotations

import math
from datetime import UTC

from epcsaft import Mixture, Parameters

T, P = 298.15, 1.0e5
COMPONENTS = ("water", "lithium-cation", "sodium-cation", "chloride-anion")
LN10 = math.log(10)


def _x(p: Parameters, m_li: float, m_background: float, charge: int) -> list[float]:
    amounts = {"water": 1 / p.molar_masses[0], "lithium-cation": m_li,
               p.component_ids[2]: m_background, "chloride-anion": m_li + charge*m_background}
    total = math.fsum(amounts.values())
    return [amounts[c] / total for c in p.component_ids]


class Activities:
    def __init__(self, p: Parameters, background: str = "NaCl"):
        ion, self.charge = {"NaCl": ("sodium-cation", 1), "MgCl2": ("magnesium-cation", 2),
                            "CaCl2": ("calcium-cation", 2)}[background]
        assert p.component_ids == ("water", "lithium-cation", ion, "chloride-anion")
        self.p, self.mixture = p, Mixture(p)
        water = [1.0, 0.0, 0.0, 0.0]
        self.refs = {salt: self.mixture.reference(T, P, solvent_fractions=water, ions=(cation, "chloride-anion"), stoichiometry=(1, z))
                     for salt, cation, z in (("LiCl", "lithium-cation", 1), (background, ion, self.charge))}
        # Offset of 1/2[ln(x rho) + mu_res/RT] for Li+ + Cl-, fixed so the binary matches the native molality-scale activity.
        self.offset = self._bracket(0.1, 0.0) - 0.5 * math.log(0.1 * 0.1) - self.refs["LiCl"].activity(0.1).mean_ionic_log_molality_scale

    def _bracket(self, m_li: float, m_background: float) -> float:
        x = _x(self.p, m_li, m_background, self.charge)
        s = self.mixture.state(T, P=P, x=x, phase="liquid")
        mu = s.residual_chemical_potential_over_rt
        return 0.5 * sum(math.log(x[j] * s.molar_density) + mu[j] for j in (1, 3))

    def end_member(self, salt: str, m: float) -> float:
        return self.refs[salt].activity(m).mean_ionic_log_molality_scale / LN10

    def trace(self, m_background: float, m_li: float = 1e-6) -> float:
        """log10 gamma+-(LiCl) on the molality scale at the declared background salt molality."""
        ln = self._bracket(m_li, m_background) - 0.5 * math.log(m_li * (m_li + self.charge*m_background)) - self.offset
        return ln / LN10

    def binary_check(self, m: float) -> float:
        """Route minus native binary activity (decade); must vanish at every m if the offset convention holds."""
        return (self._bracket(m, 0.0) - math.log(m) - self.offset) / LN10 - self.end_member("LiCl", m)

    def row(self, m: float) -> dict:
        trace, li, na = self.trace(m), self.end_member("LiCl", m), self.end_member("NaCl", m)
        return {"m": m, "log_gamma_trace": trace, "log_gamma_LiCl": li, "log_gamma_NaCl": na,
                "theta_eff": LN10 * (trace - 0.5 * (li + na)) / m}


def divalent():
    """Post-fit, report-only supplement to design §5; never rerun the fit or its scores."""
    import json
    import subprocess
    import sys
    from datetime import datetime

    import epcsaft
    import fit

    fit.verify()
    record = json.loads((fit.RESULTS / "fit.json").read_text())
    selected = next(s for s in record["starts"] if s["start"] == record["selected_start"])
    assert selected["physical"] == json.loads((fit.RESULTS / "evaluation.json").read_text())["selected"]["physical"]
    song = fit.read(fit.RESULTS / "trace-divalent-song.csv")
    rows = []
    for salt, ion in (("MgCl2", "magnesium-cation"), ("CaCl2", "calcium-cation")):
        components = ("water", "lithium-cation", ion, "chloride-anion")
        baseline = Activities(fit.parameters(components=components), salt)
        stage_a = Activities(fit.parameters(selected["physical"], components=components), salt)
        for ionic_strength in (1.0, 3.0, 6.0, 8.3):
            m = ionic_strength/3
            errors = [a.binary_check(m) for a in (baseline, stage_a)]
            assert max(map(abs, errors)) <= 1e-10
            rows.append({"background": salt, "I_background_mol_kg": ionic_strength,
                         "m_background_mol_kg": m, "m_LiCl_mol_kg": 1e-6,
                         "baseline_log_gamma_trace": baseline.trace(m),
                         "stage_A_log_gamma_trace": stage_a.trace(m),
                         "Song_log_gamma_trace": next(float(r["log_gamma_trace"]) for r in song
                                                       if float(r["I_background_mol_kg"]) == ionic_strength) if salt == "MgCl2" else None,
                         "baseline_binary_error_decade": errors[0], "stage_A_binary_error_decade": errors[1]})
    assert all(math.isfinite(v) for row in rows for v in row.values() if isinstance(v, float))
    fit.write(fit.RESULTS / "trace-divalent.csv", rows)
    (fit.RESULTS / "trace-divalent.json").write_text(json.dumps({
        "evaluated_at_utc": datetime.now(UTC).isoformat(), "selected_start": selected["start"],
        "engine_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "engine_module": epcsaft._core.__file__, "python": sys.executable,
        "temperature_K": T, "pressure_Pa": P, "base_packet": fit.D["base_packet"],
        "meaning": "Post-fit report-only supplement, design §5 and readiness r1 C6; no limits or scores",
        "ionic_strength_basis": "Background salt only; total ionic strength also includes trace LiCl",
        "inputs_sha256": {str(path.relative_to(fit.ROOT)): fit.sha(path) for path in
                          (fit.DECL, fit.RESULTS / "fit.json", fit.RESULTS / "evaluation.json", fit.RESULTS / "trace-divalent-song.csv",
                           fit.ROOT / f"data/packets/{fit.D['base_packet']}/parameters/parameters.json")}
    }, indent=2) + "\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "--divalent":
        divalent()
    else:
        p = Parameters.from_json(sys.argv[1], components=COMPONENTS)
        a = Activities(p)
        print("binary route check (decade):", [f"{a.binary_check(m):.2e}" for m in (0.01, 1, 3, 6)])
        for m in (1, 3, 6):
            print(a.row(m))
