# PHREEQC calcium-complex laws and IAPWS carbonate sources

Current native results reproduce the two representable PHREEQC calcium laws.
The [notebook](notebook.qmd) computes the comparison, rounding diagnostic and
limits from [retained results](results/calcium-comparison.csv). This is numerical
verification against the PHREEQC source model, not physical validation or a
full carbonate-speciation reproduction. Independent delivered review is pending.

Run from the checkout root after the ordinary per-checkout setup in README
Checks:

```bash
nice -n 10 .venv/bin/python analyses/reproduction/2018-iapws/scripts/reproduce.py
```

The command uses this checkout's native library, builds the pinned upstream
PHREEQC reference when absent, unsets HSL, and retains the twelve independent comparisons and execution identity.
All source-law evaluation and logarithm conversion run in C++ through the Engine;
Python invokes the programs and compares their output. PHREEQC receives the
original WATEQ4F database. Its actual pressure at the hottest points differs
from the requested pressure; no reaction-volume term changes these source laws.
The upstream v3.9.0 commit still prints a PHREEQC_3.8.9 banner, recorded in the raw
console output. Source/database/executable hashes are retained separately.

The [source directory](inputs/cases/carbonate-complexation/source/) owns IAPWS
TGD9-18 Table 2, WATEQ4F coefficient transcriptions and PHREEQC Example 1 values.
IAPWS constants are source values only. The original case files retain their
historical identities; they are not fresh Engine calculations.

Sources: [IAPWS TGD9-18](https://iapws.org/public/documents/KuswX/AIL.pdf),
[PHREEQC manual](https://doi.org/10.3133/tm6A43), SOLUTION_SPECIES p.200, and the
[pinned WATEQ4F database](https://github.com/phreeqc-dev/phreeqc3/blob/b0b3be767158ccc3322d2c816625cf470045e67e/database/wateq4f.dat),
lines identified by the retained coefficient CSV.

Full water/carbonate temperature laws wait on [#790](https://github.com/tannerpolley/ePC-SAFT/issues/790);
ideal molarity chemistry waits on [#791](https://github.com/tannerpolley/ePC-SAFT/issues/791).
Complex-species EOS parameters, a common proton/hydronium pH convention and the
molarity-to-nonideal-activity conversion remain unavailable for their distinct
reasons in `results/unavailable.csv`. The carbonate-only EOS comparison is
excluded from this delivery by the lead's Ready review on #789.
