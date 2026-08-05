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

The package owns one CMake graph and no installed cross-module SDK. Data packet
materializations are explicit and hash-bound. Installed-wheel campaigns live
in `validation/`. See `AGENTS.md`, `ARCHITECTURE.yaml`, and
`../governance/GOVERNANCE.md`.
