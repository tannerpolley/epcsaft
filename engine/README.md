# Unified ePC-SAFT Engine

One Python distribution provides EOS, equilibrium, and regression:

```python
import epcsaft
from epcsaft import equilibrium, regression
```

Build and install the current source non-editably through the root environment:

```bash
bash .codex/environments/setup.sh
(cd engine && ../.venv/bin/python -m pytest -q)
```

Setup keeps the checkout-local Ninja build warm, uses `ccache` when available,
and reuses a verified machine-local wheel keyed by source, ABI, platform, and
toolchain. `tools/check-affected.sh` runs the affected scientific domain;
`tools/check-all.sh` adds the retained slow numerical tier.

The package owns one CMake graph and no installed cross-module SDK. Data
packet materializations are explicit and hash-bound. Installed-wheel campaigns
live in `validation/`. See `AGENTS.md`, `ARCHITECTURE.yaml`, and
`../governance/GOVERNANCE.md`.

Regression exports only `ParameterFit`, `FitResult`, and `fit`:

```python
from epcsaft.regression import fit

result = fit("parameter-fit.json")
parameters = result.to_parameters()
```

The fit request owns the selected parameter identities, data hashes, units,
residuals, bounds, scales, ties, and multistart policy. Engine owns the exact
derivatives and one private native solver; Validation owns literature
acceptance.

The canonical capability cutover and its EOS, Regression, and Equilibrium owner
plans begin at
[`docs/plans/2026-08-07-general-scientific-capability-hard-cutover.md`](docs/plans/2026-08-07-general-scientific-capability-hard-cutover.md).
