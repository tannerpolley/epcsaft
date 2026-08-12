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

Run the parallel fast installed-wheel gate with:

```bash
bash tools/check-all.sh
```

The wheel is built and installed before tests begin. Slow Engine stress
evidence, literature or figure campaigns, full native clang-tidy, and
sanitizers are deliberately separate; run only the exact owner command needed
for the change, such as a canonical slice below.

Run a canonical end-to-end workflow with:

```bash
bash tools/run-slice.sh khudaida
bash tools/run-slice.sh figiel
```

See `AGENTS.md` for agent routing and `governance/GOVERNANCE.md` for the
scientific and architectural doctrine.
