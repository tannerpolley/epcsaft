# Species-reference numerical verification

The native comparisons follow `neutral_reaction_reference_transfer`,
`reference_single_insertion_taylor` and ALG-REFERENCE. `reference_qualification.cpp`
declares the CO2 parameters used by the Uyan generator, the identical-ethane
reference cases and the existing record-anchored ethane caloric case. It calls
the Engine's native reference and pure-state calculations. Centered differences
only verify the exact derivatives. The accepted column in `qualification.tsv`
marks the refined derivative values and terminal gas limits; preceding values
retain the observed refinement. The oversized pure species is a manufactured
density-domain refusal, not a physical fluid prediction. The caloric case declares the
pure ethane liquid: at 260 K and 3 or 5 MPa the native isotherm has only its
liquid-side root, so a declared vapor is refused.

`reference_snapshot.cpp` serializes each existing reference value, lower/higher
column, ladder, change, criterion and availability as hexadecimal doubles.
The literal field order is in `dump`. The retained base output is
`results/data/native/species-reference/default-native-values.csv`; the candidate output has identical bytes.
The macro enables the full explicit default-map comparison without changing
the baseline serialization. `wheel_snapshot.py` reads the retained named
base-wheel values and compares omitted, empty and full default maps, including
first temperature, pressure and EOS-coordinate actions.

Run from the checkout root with its own non-editable wheel:

```bash
bash tools/setup.sh
bash analyses/reproduction/2015-uyan/scripts/species-reference/qualify.sh
.venv/bin/python analyses/reproduction/2015-uyan/scripts/qualify_species_reference.py
.venv/bin/python -m pytest -q engine/tests/equilibrium/test_neutral_reference.py
```

The historical diagnostic reuses the same native pressure starts and physical
input guesses as the source-ledger generator, with only the explicit CO2
reference declaration added. Source constants stay raw. Earlier tables and
their unknown wheel identity remain separate. The retained key-solve summary and
work table identify both wheel hashes and the merged base, and retain work
counts, state differences and paired timings. The original base wheel also remains under `build/issue-274/base` in
the delivery worktree; it is not committed.

Presentation reads retained values and never evaluates the model:

```bash
.venv/bin/python analyses/reproduction/2015-uyan/scripts/present_species_reference.py
quarto render analyses/reproduction/2015-uyan/notebook.qmd --to html --no-execute
```

The evidence is numerical verification under the declared inputs. Finite
reference changes do not certify infinite-dilution bias. Delivered review
and manuscript promotion have not occurred.
