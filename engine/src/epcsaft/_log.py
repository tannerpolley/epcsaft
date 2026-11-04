"""The ``epcsaft`` logger: quiet unless configured, fed only from records results already carry.

INFO writes one line per public solve or fit. DEBUG adds one line per Ipopt iteration of
an equilibrium solve, collected natively only while DEBUG is enabled (no trace capture),
and one line per fit optimizer iteration from ``iteration_costs``/``iteration_seconds``.
"""

import logging
import time
from collections.abc import Callable
from contextlib import nullcontext
from typing import Any, TypeVar

from . import _core

R = TypeVar("R")
log = logging.getLogger("epcsaft")
log.addHandler(logging.NullHandler())

_IPOPT = ("inf_pr", "inf_du", "mu", "step_norm", "regularization", "alpha_pr", "alpha_du")
_IPOPT_LINE = "%s Ipopt run %d iteration %d (%s): objective=%.17g, %s, line_search_trials=%d"
# Native result fields, read when present: equilibrium, multistart, continuation and fit results.
_SUMMARY = (
    "success status usable selected budget_exhausted solver_status iterations evaluations residual_evaluations "
    "jacobian_evaluations accepted_steps final_cost cause_stop"
).split()


def _summary(result: Any) -> str:
    fields = [f"{name}={getattr(result, name)}" for name in _SUMMARY if hasattr(result, name)]
    if hasattr(result, "convergence"):
        fields += [f"stop={result.convergence.stop.name}", f"cause={result.convergence.cause or '-'}"]
    elif hasattr(result, "stop"):
        fields.append(f"stop={result.stop}")
    if hasattr(result, "causes") and result.selected is None:  # multistart: why no start was selected
        fields.append(f"causes={result.causes}")
    message = getattr(result, "message", "")
    return ", ".join(fields) + (f": {message.split('; ')[0]}" if message else "")


def logged(solve: Callable[..., R], *args: Any) -> R:
    """Run one public native solve and log it; returns the native result unchanged."""
    if not log.isEnabledFor(logging.INFO):
        return solve(*args)
    name = solve.__name__
    debug = log.isEnabledFor(logging.DEBUG)
    begun = time.perf_counter()
    with _core.IpoptIterations() if debug else nullcontext() as collected:
        result = solve(*args)
    seconds = time.perf_counter() - begun
    for index, item in enumerate(result if isinstance(result, list) else [result]):
        label = f"{name} start {index}" if isinstance(result, list) else name
        log.info("%s: %s; %.6g s", label, _summary(item), seconds)
        for k, cost in enumerate(getattr(item, "iteration_costs", ())):
            log.debug("%s iteration %d: cost=%.17g, seconds=%.6g", label, k, cost, item.iteration_seconds[k])
    for record in collected.records if collected else ():
        i = record.iteration
        values = ", ".join(f"{field}={getattr(i, field):.6e}" for field in _IPOPT)
        mode = "restoration" if i.mode else "regular"
        log.debug(_IPOPT_LINE, name, record.run, i.iteration, mode, i.objective, values, i.line_search_trials)
    return result
