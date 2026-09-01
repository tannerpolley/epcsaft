# Native Regression quickstart

Regression accepts compiled EOS or Equilibrium owners and ordered observed
rows.  Owners retain thermodynamic values, exact physical Jacobians, and
structured failure diagnostics; Ceres receives only native owner handles.

```python
import epcsaft
from epcsaft import eos, regression

def fit_pressure(model, parameters, active_identity, observation_pa, bounds):
    active = epcsaft.ActiveParameterSet(parameters, (active_identity,))
    owner = eos.compile_regression_evaluator(
        {
            "model": model,
            "family": "fixed-state",
            "temperature_k": 310.0,
            "molar_density_mol_m3": 1050.0,
            "mole_fractions": (0.4, 0.6),
            "row_identities": ("pressure_pa",),
        },
        active,
    )
    return regression.least_squares(
        (owner,),
        observations=(observation_pa,),
        initial_parameters=active.values,
        bounds=bounds,
        affine_origin=(0.0,),
        affine_scale=(1.0,),
        residual_scales=(1000.0,),
    )
```

`result.parameters`, `result.predictions`, `result.residuals`, and
`result.physical_jacobian` are physical outputs. `result.optimizer_jacobian`
is independently retained after the exact affine and residual-scale chain
rule. `result.counters`, `result.failures`, `result.covariance`, and native
solver options preserve the complete diagnostic result. `curve_fit` is a
small alias for the same owner-neutral core.

Equilibrium owners are compiled with
`equilibrium.compile_regression_evaluator` only for anchored local D-034 and
D-035 declarations. Unsupported topology, missing anchors, and missing local
starts refuse before Ceres is constructed.
