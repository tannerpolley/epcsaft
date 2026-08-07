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

accepted = result.accepted_start
if accepted is None:
    raise RuntimeError("no usable fit start")

print(result.fitted_values())
print(accepted.rank, accepted.condition_number, accepted.active_bounds)

# Explicit outputs are optional; fit never mutates the input parameter file.
result.to_json("propane-baseline.fit-result.json")
Path("propane-fitted.parameters.json").write_text(
    json.dumps(result.to_parameters().to_mapping(), indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)
```

`fit` also accepts a strict `ParameterFit` or matching in-memory mapping.
`FitResult` retains every start, residuals, exact Jacobians, rank, conditioning,
active bounds, basin identity, and artifact fingerprints. The sole public
exports are `ParameterFit`, `FitResult`, and `fit`.

The initial cutover admits pure-saturation and fixed-state-pressure
observations backed by exact EOS selected-parameter derivatives. It does not
admit mixture phase-equilibrium, charged/electrolyte, direct permittivity or
Born, or reaction observations.
