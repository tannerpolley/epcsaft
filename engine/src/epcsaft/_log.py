"""The ``epcsaft`` logger: quiet unless configured, fed only from records results already carry.

INFO writes one line per public solve or fit. DEBUG adds one line per Ipopt iteration of
an equilibrium solve, collected natively only while DEBUG is enabled (no trace capture),
and one line per fit optimizer iteration from ``iteration_costs``/``iteration_seconds``.
"""

import logging
import time
from collections.abc import Callable
from contextlib import nullcontext
from typing import Any

from . import _core

log = logging.getLogger("epcsaft")
log.addHandler(logging.NullHandler())

_SUMMARY = ("success", "status", "usable", "selected", "solver_status", "iterations", "evaluations", "final_cost")


def _summary(result: Any) -> str:
    fields = [f"{name}={getattr(result, name)}" for name in _SUMMARY if hasattr(result, name)]
    if hasattr(result, "convergence"):
        fields += [f"stop={result.convergence.stop.name}", f"cause={result.convergence.cause or '-'}"]
    elif hasattr(result, "stop"):
        fields.append(f"stop={result.stop}")
    message = getattr(result, "message", "")
    return ", ".join(fields) + (f": {message.split('; ')[0]}" if message else "")


def logged(name: str, call: Callable[[], Any]) -> Any:
    """Run one public solve and log it; returns the native result unchanged."""
    if not log.isEnabledFor(logging.INFO):
        return call()
    debug = log.isEnabledFor(logging.DEBUG)
    begun = time.perf_counter()
    with _core.IpoptIterations() if debug else nullcontext() as collected:
        result = call()
    seconds = time.perf_counter() - begun
    for index, item in enumerate(result if isinstance(result, list) else [result]):
        label = f"{name} start {index}" if isinstance(result, list) else name
        log.info("%s: %s; %.6g s", label, _summary(item), seconds)
        for k, (cost, wall) in enumerate(zip(getattr(item, "iteration_costs", ()), getattr(item, "iteration_seconds", ()))):
            log.debug("%s iteration %d: cost=%.17g, seconds=%.6g", label, k, cost, wall)
    for record in collected.records if debug else ():
        i = record.iteration
        log.debug("%s Ipopt run %d iteration %d (%s): objective=%.17g, inf_pr=%.6e, inf_du=%.6e, mu=%.6e, step_norm=%.6e, "
                  "regularization=%.6e, alpha_pr=%.6e, alpha_du=%.6e, line_search_trials=%d", name, record.run, i.iteration,
                  "restoration" if i.mode else "regular", i.objective, i.inf_pr, i.inf_du, i.mu, i.step_norm,
                  i.regularization, i.alpha_pr, i.alpha_du, i.line_search_trials)
    return result
