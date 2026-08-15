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
`FitResult` retains every start, residuals, exact Jacobians, rank, conditioning,
active bounds, parameter-space optimizer basin/cluster identity, and artifact
fingerprints. This optimizer cluster is not a thermodynamic physical-branch
identity.

`to_parameters()` is an explicit authority-neutral candidate export. A
downstream scientific campaign, not Regression, decides whether to adopt it.

The admitted families are EOS pure-saturation/fixed-state observations and the
Equilibrium-owned fixed-T,P homogeneous reactive and fixed-topology reactive
bubble/VLE families. Both reactive families use immutable local continuation
anchors, exact total parameter Jacobians, and typed non-evaluable trials. The
bubble family treats reported partial pressure as `y_i P` while phase
equilibrium uses the Provider fugacity. Neither family makes a global-equilibrium
or parameter-promotion claim. Nonreactive mixture phase-equilibrium, direct
permittivity, and standalone Born observations are not admitted.
