# Reaction references and neutral association

Numerical verification for #922, governed by the main checkout's
`build/issue-919/design-e3e4.md` r2 and `readiness-e3e4-r2.md`.
The lead's 2026-10-29 amendment permits 360 E3 gross executable lines within
620 total, with 100 test lines separate. Tolerances are unchanged.

The retained E3 comparisons pass. The first E4 MDEA comparison misses its
20% held-out improvement requirement in one pressure split. E4 is incomplete;
no Engine/default/packet change or physical improvement is adopted.
Read [the notebook](notebook.qmd) for numbers, interpretation and remaining gates.

In this checkout's prepared environment:

```bash
.venv/bin/python analyses/2026-reaction-reference-association/scripts/declare.py
nice -n 10 .venv/bin/meson compile -C build/engine-meson/shared -j 6 reaction_reference_e3
nice -n 10 build/engine-meson/shared/reaction_reference_e3 analyses/2026-reaction-reference-association/inputs/mdea.txt
nice -n 10 build/engine-meson/shared/reaction_reference_e3 analyses/2026-reaction-reference-association/inputs/mea.txt mea
nice -n 10 .venv/bin/python analyses/2026-reaction-reference-association/scripts/e4.py
.venv/bin/python analyses/2026-reaction-reference-association/scripts/present.py
```

The native CSV streams are retained as `results/e3.csv` and `results/mea.csv`
before presentation. `e4.py` runs only the MDEA neutral comparison; its result
explicitly records that the full E4 study is unqualified. Native Engine owners
perform all property calculations, equilibrium solves and fitting iterations.

E5 uses the frozen E4 coefficients without fitting. Both 3B reactive pressure
comparisons fail the four-row improvement and other-row preservation limits;
see the notebook and `results/e5-summary.json`. Its incremental allowance is
360 gross executable lines against `48b23a091`, with no Engine or test additions.
Generate frozen inputs with `declare.e5_inputs()` (the copied baseline rows must
exist under `build/issue-922/e5/`), run `reaction_reference_e3` on each
`inputs/e5-*-neutral.txt`, and on `inputs/e5-*-reactive.txt` with argument `e5`.
Retain the streams as `results/e5-neutral-*.csv` and `results/e5-reactive-*.csv`,
then call `present.e5_present()`. These functions declare inputs or read retained
values; every EOS and equilibrium calculation runs in the native executable.

The E5 evidence copy of `notebook.html` uses external assets in `notebook_files/`
and `results/`; keep these beside the HTML. Resource embedding failed (see the
evidence `render.log`), so this HTML is not self-contained.
