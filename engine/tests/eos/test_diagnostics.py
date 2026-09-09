import hashlib
import json
from collections.abc import Callable
from functools import partial
from pathlib import Path
from typing import Any

import pytest

from epcsaft import Convergence, Mixture, Parameters, Stop, trace
from tests.packets import parameters_json

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]


def _parameter_metadata(path: Path, packet_id: str, version: int, components: tuple[str, ...]) -> dict[str, object]:
    supplied_mapping = json.loads(path.read_text())
    return {
        "packet_id": packet_id,
        "packet_version": version,
        "parameter_file": {
            "path": str(path.relative_to(REPOSITORY_ROOT)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "size_bytes": path.stat().st_size,
        },
        "effective_parameters": {
            "supplied_mapping": supplied_mapping,
            "selected_components": list(components),
        },
    }


def _stopping_criterion(*, threshold: float, verdict: bool, name: str) -> dict[str, object]:
    return {
        "type": "custom",
        "helper": "engine.tests.eos.test_diagnostics._stopping_criterion",
        "rule": "residual(model) <= stopping criterion(reference)",
        "verdict": verdict,
        "threshold": threshold,
        "scale": "absolute residual",
        "rationale": "Use the existing density-root solver stopping criterion.",
        "criterion_name": name,
    }


@pytest.mark.evidence(
    test_id="python.diagnostics.trace-convergence-round-trip",
    title="Trace and convergence round trip",
    kind="Numerical",
    category="Bindings/Python",
    requirement_id="B1",
    execution_policy="Automatic fast",
    lifecycle="Implemented",
    evidence_expectation="Numerical",
    expected_comparisons=(
        "trace_identity/residual",
        "trace_identity/iterations",
        "trace_identity/residual_criterion",
    ),
    basis="Relationship",
    method="public trace convergence round trip",
    claim="Captured and uncaptured public convergence diagnostics retain the same residual and iteration values.",
)
def test_trace_and_convergence_round_trip(evidence: Any) -> None:
    components = ("methane", "ethane")
    temperature = 300.0
    pressure_target = 1.0e5
    composition = (0.5, 0.5)
    phase = "liquid"
    parameter_path = parameters_json("gross-2001-methane-ethane", 2)
    parameters = Parameters.from_json(parameter_path, components=components)
    mixture = Mixture(parameters)
    with trace() as captured:
        state = mixture.state(
            temperature, P=pressure_target, x=composition, phase=phase
        )
    convergence = state.density_convergence
    assert convergence is not None
    assert isinstance(convergence, Convergence)
    root = next(span for span in captured.spans if span.parent is None)
    assert root.convergence is not None
    inputs = {
        "components": list(components),
        "T_K": temperature,
        "pressure_target_Pa": pressure_target,
        "composition": list(composition),
        "phase": phase,
        "_parameters": _parameter_metadata(
            parameter_path,
            "gross-2001-methane-ethane",
            2,
            components,
        ),
    }
    exact = {
        "type": "scalar",
        "rule": "symmetric_close",
        "atol": 0.0,
        "rtol": 0.0,
        "scale": "max(abs(model), abs(reference))",
        "combining_rule": "max(atol, rtol*scale)",
    }
    evidence.compare(
        "trace_identity/residual",
        quantity="density root residual",
        model=float(root.convergence.residual),
        reference=float(convergence.residual),
        unit="1",
        inputs=inputs,
        criterion=exact,
        basis="Relationship",
        method="trace convergence round trip",
        claim="The root span exposes the same residual as the returned state convergence record.",
        relationship_latex=r"r_{trace} = r_{state}",
        quantity_latex="r",
        inputs_latex=r"T,P,\mathbf{x},\mathrm{phase}",
    )
    evidence.compare(
        "trace_identity/iterations",
        quantity="density root iterations",
        model=float(root.convergence.iterations),
        reference=float(convergence.iterations),
        unit="count",
        inputs=inputs,
        criterion=exact,
        basis="Relationship",
        method="trace convergence round trip",
        claim="The root span exposes the same iteration count as the returned state convergence record.",
        relationship_latex=r"k_{trace} = k_{state}",
        quantity_latex="k",
        inputs_latex=r"T,P,\mathbf{x},\mathrm{phase}",
    )
    if convergence.stop == Stop.converged:
        stop_criterion = convergence.criterion
        if stop_criterion is not None:
            evidence.compare(
                "trace_identity/residual_criterion",
                quantity="density root residual against stopping criterion",
                model=float(convergence.residual),
                reference=float(stop_criterion.value),
                unit="1",
                inputs={**inputs, "criterion_name": stop_criterion.name},
                criterion=_stopping_criterion(
                    threshold=float(stop_criterion.value),
                    verdict=convergence.residual <= stop_criterion.value,
                    name=stop_criterion.name,
                ),
                basis="Relationship",
                method="density root stopping criterion",
                claim="A converged density root satisfies its recorded stopping criterion.",
                relationship_latex=r"r \leq r_{criterion}",
                quantity_latex="r",
                inputs_latex=r"T,P,\mathbf{x},\mathrm{phase},r_{criterion}",
            )
        else:
            evidence.observe(
                "trace_identity/residual_criterion",
                quantity="density root residual",
                model=float(convergence.residual),
                unit="1",
                inputs={**inputs, "stop": str(convergence.stop), "criterion": "missing"},
                basis="Relationship",
                method="density root convergence diagnostic",
                claim="A converged result retains its residual when the solver criterion is unavailable.",
                relationship_latex=r"r_{observed}",
                quantity_latex="r",
                inputs_latex=r"T,P,\mathbf{x},\mathrm{phase}",
            )
    else:
        evidence.observe(
            "trace_identity/residual_criterion",
            quantity="density root residual",
            model=float(convergence.residual),
            unit="1",
            inputs={**inputs, "stop": str(convergence.stop), "criterion": "not available"},
            basis="Relationship",
            method="density root convergence diagnostic",
            claim="The non-converged residual is retained without inventing a stopping criterion.",
            relationship_latex=r"r_{observed}",
            quantity_latex="r",
            inputs_latex=r"T,P,\mathbf{x},\mathrm{phase}",
        )
    assert root.convergence.residual == convergence.residual
    assert root.convergence.iterations == convergence.iterations
    assert root.convergence.stop == convergence.stop
    if convergence.stop == Stop.converged:
        assert stop_criterion is not None
        assert convergence.residual <= stop_criterion.value
    assert "DENSITY_ROOT_RESIDUAL" in repr(convergence)


@pytest.mark.evidence(
    test_id="python.surface.reference-requires-stoichiometry",
    title="Reference requires stoichiometry",
    kind="Package",
    category="Surface",
    requirement_id="S1",
    execution_policy="Automatic fast",
    lifecycle="Implemented",
    evidence_expectation="Outcome only",
)
def test_reference_requires_stoichiometry() -> None:
    parameters = Parameters.from_json(
        parameters_json("figiel-2025-reference-electrolytes", 3),
        components=("water", "lithium-cation", "chloride-anion"),
    )
    with pytest.raises(TypeError):
        reference: Callable[..., object] = Mixture(parameters).reference
        call = partial(
            reference,
            298.15,
            1.0e5,
            solvent_fractions=(1.0, 0.0, 0.0),
            ions=("lithium-cation", "chloride-anion"),
        )
        call()
