"""Fixed numerical acceptance, shared by timed-return checks and the report."""

import math
from copy import deepcopy

PRESSURE_RESIDUAL_LIMIT = 1e-12
DENSITY_RELATIVE_LIMIT = 1e-10


def density_values(result, fixed=False):
    values = dict(
        rho_mol_m3=result.molar_density,
        pressure_Pa=result.pressure,
        slope_Pa_m3_mol=result.pressure_density_derivative,
    )
    if not fixed:
        values.update(density_branch=result.density_branch.name, stable_root_count=result.stable_root_count)
    return values


def density_mutations(value):
    for field, replacement in (
        ("rho_mol_m3", 2),
        ("slope_Pa_m3_mol", -1),
        ("density_branch", "Liquid"),
        ("stable_root_count", 99),
        ("pressure_Pa", 2),
    ):
        changed = deepcopy(value)
        changed[field] = changed[field] * replacement if isinstance(changed[field], float) else replacement
        yield field, changed


def difference(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), "result fields differ"
        return max(
            (
                max((abs(x - y) for x, y in zip(a[k], b[k], strict=True)), default=0.0)
                if k == "tpd_over_rt"
                else difference(a[k], b[k])
                for k in a
            ),
            default=0.0,
        )
    if isinstance(a, list):
        assert len(a) == len(b), "result lengths differ"
        if all(isinstance(v, (int, float)) for v in a + b):
            assert all(math.isfinite(v) for v in a + b), "nonfinite result"
            scale = max(map(abs, a + b), default=0.0)
            return max(
                (
                    abs(x - y) / scale if scale else 0.0
                    for x, y in zip(a, b, strict=True)
                ),
                default=0.0,
            )
        return max((difference(x, y) for x, y in zip(a, b, strict=True)), default=0.0)
    if isinstance(a, (int, float)) and not isinstance(a, bool):
        assert math.isfinite(a) and math.isfinite(b)
        return abs(a - b) / max(abs(a), abs(b)) if a or b else 0.0
    assert a == b, f"{a!r} != {b!r}"
    return 0.0


def density_check(state, value, gas_constant, *, fixed=False):
    assert all(
        math.isfinite(value[k])
        for k in ("rho_mol_m3", "pressure_Pa", "slope_Pa_m3_mol")
    )
    delta = abs(value["rho_mol_m3"] / state["oracle_rho_mol_m3"] - 1)
    assert delta <= DENSITY_RELATIVE_LIMIT, f"density difference {delta}"
    assert value["slope_Pa_m3_mol"] > 0, "nonpositive pressure slope"
    if not fixed:
        expected = (
            "Unspecified"
            if state["oracle_stable_root_count"] == 1
            else state["branch"].capitalize()
        )
        assert value["density_branch"] == expected, (
            f"density branch {value['density_branch']} != {expected}"
        )
        assert value["stable_root_count"] == state["oracle_stable_root_count"], (
            "stable root count"
        )
    residual = abs(value["pressure_Pa"] - state["P_Pa"]) / (
        value["rho_mol_m3"] * gas_constant * state["T_K"]
    )
    assert residual <= PRESSURE_RESIDUAL_LIMIT, (
        f"reduced pressure residual {residual}"
    )
    return delta


def application_check(value, reference):
    if isinstance(value, dict) and "compositions" in value:
        # Stability minima have no trial-order identity after #647's Michelsen search.
        value, reference = deepcopy(value), deepcopy(reference)
        for result in (value, reference):
            pairs = sorted(zip(result["compositions"], result["tpd_over_rt"], strict=True))
            result["compositions"], result["tpd_over_rt"] = map(list, zip(*pairs, strict=True))
    delta = difference(value, reference)
    assert delta <= 1e-10, f"application difference {delta}"
    return delta
