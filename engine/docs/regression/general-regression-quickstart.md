# Parameter Fit quickstart

Regression has one source-neutral public API:

```python
import json
from pathlib import Path

from epcsaft.regression import fit

request = Path(
    "validation/campaigns/2026-pure-saturation-regression/"
    "propane-baseline.fit.json"
)
result = fit(request)

usable = result.best_usable_start
if usable is None:
    raise RuntimeError("no usable fit start")

print(result.fitted_values())
print(usable.rank, usable.condition_number, usable.active_bounds)

# Explicit outputs are optional; fit never mutates the input parameter file.
result.to_json("propane-baseline.fit-result.json")
Path("propane-candidate.parameters.json").write_text(
    json.dumps(result.to_parameters().to_mapping(), indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
```

The `best_usable_start` property means optimizer-usable start; it is not
a scientific acceptance or parameter-adoption status.

`fit` also accepts a strict `ParameterFit` or matching in-memory mapping.
Each declared optimizer coordinate uses
`physical = origin + scale * affine`; callers provide the exact ordered set of
physical starts with units implied by their parameter coordinates and stable
start identities. Regression performs the affine conversion once internally.
No starts are generated implicitly.

Optimizer stopping controls are optional values in the existing solver mapping:

```python
request["solver"]["controls"] = {
    "maximum_iterations": 25,
    "function_tolerance": 1.0e-6,
    "gradient_tolerance": 1.0e-6,
    "parameter_tolerance": 1.0e-6,
}
```

Those values are the defaults. Ipopt stops when a convergence tolerance is
reached or the iteration ceiling is exhausted. These controls do not alter EOS,
phase-equilibrium, reaction, density, or certification tolerances.

`FitResult` retains every start, residuals, exact Jacobians, rank, conditioning,
active bounds, parameter-space optimizer basin/cluster identity, and artifact
fingerprints. This optimizer cluster is not a thermodynamic physical-branch
identity.

`to_parameters()` is an explicit authority-neutral candidate export. A
downstream scientific campaign, not Regression, decides whether to adopt it.

The observation families are supplied by EOS and Equilibrium owners. Configured
generic-equilibrium problems retain their declared phases, outputs, local
continuation anchors, exact total parameter Jacobians, and typed non-evaluable
trials. Regression only maps matching data to those outputs and minimizes the
declared residuals; it does not select a phase topology or a case-specific
route.
