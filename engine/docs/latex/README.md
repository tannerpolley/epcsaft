# EOS equation mirror

`equations.tex` is the exhaustive human-readable mirror of mathematical
relations executed by the current EOS. Runtime C++ and Python remain the
equation owners, and `docs/science/full-eos.md` remains the normative
scientific contract.

The legacy catalogue is retained only as provenance and structural guidance:

- commit `7920e56dbb0fc1fb418bae44e67ef50f334b75cd`
- Git blob `2e970476258558de05391d9d77f31d1736d78de8`
- legacy file SHA-256
  `a0f2f3aff12042a8a50c5254c2af23e04b60af3e886b97b2756e6ecad800bacc`

Check synchronization from the repository root:

```bash
python scripts/docs/sync_equation_registry.py --check
```

Optionally render a local PDF:

```bash
cd docs/latex
latexmk -pdf -outdir=out equations.tex
```

The `out/` directory and generated PDFs are not tracked. CppAD remains the sole
nonlinear production derivative authority; this mirror does not introduce a
second derivative implementation.
