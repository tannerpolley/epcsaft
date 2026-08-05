# ePC-SAFT

Private scientific monorepo for the unified ePC-SAFT engine and its governed
data, numerical evidence, validation campaigns, and research oracles.

## Start here

```bash
bash .codex/environments/setup.sh
bash tools/status.sh
bash tools/check-affected.sh
```

The setup creates one root `.venv`, builds `engine/` as a non-editable
`epcsaft` wheel, and installs the tools needed by every directory. It never
installs separate EOS, Equilibrium, or Regression distributions.

Run all retained scientific evidence and readiness checks with:

```bash
bash tools/check-all.sh
```

Run a canonical end-to-end workflow with:

```bash
bash tools/run-slice.sh khudaida
bash tools/run-slice.sh figiel
```

See `AGENTS.md` for agent routing and `governance/GOVERNANCE.md` for the
scientific and architectural doctrine.
