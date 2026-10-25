# The like-ion dispersion convention in ePC-SAFT

Audit of every paper in the project's Zotero library on the ePC-SAFT lineage, to settle
which engine is right about the term the crosswalk isolated. **The literature is
unanimous and explicit: dispersion between two ions of the same sign is zero.**

Sources are the library's own markdown transcriptions under `~/Zotero/storage`; line
numbers are into those files.

## The rule, stated by the papers

| paper | what it says about ion–ion dispersion |
|---|---|
| **Cameretti, Sadowski & Mollerup 2005** — the founding ePC-SAFT | Applies PC-SAFT's "hard-core, dispersion, and association interactions … already taken into account by the original PC-SAFT", adding only the charging term (:45). States no ion–ion rule itself; the model is the one Held 2008 then describes. |
| **Held, Cameretti & Sadowski 2008** | The rule, in one sentence (:109): *"Eq. (6) is applied for interactions between water and ions only. **Dispersion between two ions is neglected in this work.**"* All ion–ion dispersion off, cation–anion included. No `k_ij` at all, because `u_ion` was already fitted against salt/water. |
| **Held et al. 2014, ePC-SAFT Revised** | The one deliberate change (:12): *"a new modeling approach is suggested that accounts explicitly also for **dispersion interactions between anions and cations**. This yields a much more precise description of electrolyte solutions at higher concentrations."* Cation–anion comes back, with a fitted `k_ij`. Like-ion is not part of what returns. |
| **Bülow, Ascani & Held 2020 / 2021, ePC-SAFT advanced I & II** | The change is the concentration-dependent dielectric constant `ε_r(x̄)`; the combining rules stay "Berthelot-Lorentz … Eqs. (12) and (13)" (:151). The ion–ion rule is not revisited. |
| **Figiel, Yu & Held 2025** — the packet's own source | Decisive (:50): *"the dispersion energy between **cation-cation pairs and anion-anion pairs is set to zero** in ePC-SAFT (`u_ii = 0`)"*. |
| **Schick et al. 2023; Khudaida et al. 2026** | Both use the same ePC-SAFT with a `k_ij` on the cross terms and neither restates the rule, which is consistent with it being settled. |

So the convention has one change in twenty years, and it is not the one at issue:

- **2005–2013**: no ion–ion dispersion of any kind.
- **2014 onward**: cation–anion dispersion, with a fitted `k_ij`. **Like-ion stays zero.**

## The packet proves it independently

Figiel's eq 3 is the Chen–Kreglewski temperature-dependent diameter
`d_i(T) = σ_i (1 − 0.12 exp(−3 u_ii / k_B T))`. Setting `u_ii = 0` collapses it to eq 4,
`d_ion = 0.88 σ_ion` — a *derived consequence* of the rule, not an independent choice.

Every ion in `figiel-2025-reference-electrolytes/2` satisfies it exactly:

| ion | σ | packing diameter | 0.88 σ |
|---|---|---|---|
| lithium-cation | 2.8449 | 2.503512 | 2.503512 |
| sodium-cation | 2.8232 | 2.484416 | 2.484416 |
| potassium-cation | 3.3417 | 2.940696 | 2.940696 |
| chloride-anion | 2.7560 | 2.425280 | 2.425280 |
| bromide-anion | 3.0707 | 2.702216 | 2.702216 |

The packet the greenfield engine reads is built on `u_ii = 0`.

## Verdict

**The superseded engine is right and the greenfield engine is wrong.** The greenfield
engine applies the full combining rule to every pair including cation–cation and
anion–anion; ePC-SAFT sets those to zero. That is precisely the `x_ion²` term the crosswalk
isolated, and it accounts for the whole disagreement:

- the two engines respond identically to a cation–anion `k_ij` perturbation (to ten
  figures), because both include that term — correct under Held 2014;
- they differ only in the like-ion entries, because only one of them zeroes them.

The engine is also **internally inconsistent**: it honours `u_ii = 0` in the *diameter*
(`physics/helmholtz.hpp`, `EqID: ionic_packing_diameter`, which takes
`packing_diameter_angstrom` for a charged component instead of evaluating eq 3) while
ignoring `u_ii = 0` in the dispersion energy. One consequence of the rule is implemented and
the rule itself is not.

## What has to change

In the dispersion double sum, the entry for a pair `(i, j)` where `charge_i ≠ 0` and
`charge_j ≠ 0` and `sign(charge_i) == sign(charge_j)` — which includes every diagonal ion
entry — takes `u_ij = 0`. Cation–anion keeps `u_ij = (1 − k_ij) √(u_i u_j)`.

This is a **Stage 2 defect in `native/physics/helmholtz.hpp`**, surfaced by Stage 4. It is
not a Stage 4 defect: the parameter builder maps every field correctly, which the crosswalk
established separately.

Stage 2's gates did not catch it and could not have. They are closed forms, exact
symmetries, limits and derivative agreement — all of which hold for any `u_ij` matrix,
including a wrong one. The rule is a modelling convention, and the only thing that sees a
wrong convention is a comparison against the model that defined it. That is an argument for
the oracle, not against the property tests.

## Addendum: the rule is already registered, and the packets disagree

Two things found after the audit above, which change *where* the fix belongs.

**The equation is already in the registry**, `same_sign_ionic_dispersion`, with both branches
and a citation for each:

```text
ε_ij = √(ε_i ε_j)(1 − k_ij),  z_i z_j ≤ 0, or all-ion-pair dispersion
     = 0,                      z_i z_j > 0 and same-sign exclusion
```

Source: `\cite{Figiel2025}` for the exclusion, `\cite{Ji2012}` for all-ion-pair dispersion.
Status "Implementation convention", change note *"The caller selects all ion pairs or
same-sign exclusion as part of the electrolyte model family."* And
`physics/helmholtz.hpp:240` already carries the `// EqID: same_sign_ionic_dispersion`
annotation on the line that computes `epsilon_ij`.

**So the rule was designed as a declared model-family choice, and only one branch was ever
written.** The annotation points at code that implements the other branch unconditionally.
This is not a kernel that forgot a rule; it is a selector that was never built.

**The packets encode it inconsistently, as data:**

| packet | same-sign pairs declared | value |
|---|---|---|
| `held-cameretti-sadowski-2008/2` | 198 of 435 | `k_ij = 1.0` — the exclusion, written as data |
| `figiel-2025-validation-parameters/2` | 6 of 21 | `k_ij = 0.0` — full dispersion |
| `mea-co2-h2o-nine-species-estimation/1` and `/2` | 7 of 21 | `k_ij = 0.0` |
| `figiel-2025-reference-electrolytes/2` | 0 of 15 | absent, so the matrix default 0.0 applies |
| `figiel-2025-aqueous-iodides/2`, `figiel-yu-held-2025/2`, `khudaida-2026-…/5` | 0 | absent |

The Held 2008 packet does exactly what a k_ij override would do, by hand, 198 times. Every
other electrolyte packet either omits the pairs or declares `0.0`, and under the current
kernel both mean **full like-ion dispersion** — which contradicts the paper each of them
comes from. The MEA bundle's seven `0.0` entries are in the same class.

**Consequence for the fix.** Forcing `k_ij = 1` for same-sign pairs is the right mechanism —
the loop at `helmholtz.hpp:228-260` runs over `i` and `j` including `i == j` and applies
`(1 − k_ij)` with no diagonal special case, so it zeroes the term exactly. But applying it
unconditionally would hardcode one branch of an equation the registry says is a caller
choice, and would silently overwrite the `0.0` that two packets declare on purpose.

Three constraints on where it goes:

1. **It must reach the native fixtures.** `tests/native/model_support.hpp` builds `Parameters`
   structs directly and never passes through `parameter_builder.cpp`, so a builder-only fix
   leaves `aqueous_electrolyte()` and `nine_species()` on the wrong convention.
2. **The temperature slope must go with it.** `k_ij` is
   `k_ij + slope·(1/T − 1/T_ref)` when a reference temperature is set
   (`helmholtz.hpp:237-241`); setting the constant to 1 without zeroing the slope leaves
   `ε_ij ≠ 0` at every temperature but `T_ref`.
3. **No silent default.** `engine/AGENTS.md` forbids hidden physical defaults. An electrolyte
   parameter set whose same-sign entries are neither declared nor resolved by a declared
   choice should be rejected by name, not quietly given one convention or the other.

## Addendum 2: measured against Figiel's own published curve

The paper's Figure 5 model line for LiBr is now digitized —
`analyses/reproduction/2025-figiel/figures/figure-05/source/libr-model-digitized.csv`, 27
points over 2.160–5.847 mol kg⁻¹, per-point uncertainty 0.0078–0.0122 in γ (0.45–0.8 %
relative), from the lossless embedded raster rather than a screen capture. It is a
digitization and is not an oracle row.

**The engine is 49.19 % AARD against it**, 27 of 27 points, monotone from −32.6 % at
2.16 mol kg⁻¹ to −61.2 % at 5.85. The digitization uncertainty is 0.6 %, so the gap is about
**82× the digitization error** — nothing here is an artifact of tracing a line off a page.

| m | Figiel's model | this engine | deviation |
|---|---|---|---|
| 2.160 | 1.161 | 0.7820 | −32.6 % |
| 3.168 | 1.482 | 0.8437 | −43.1 % |
| 4.340 | 1.935 | 0.9246 | −52.2 % |
| 5.847 | 2.693 | 1.0436 | −61.2 % |

The density solve raised no error at any of the 27 points, so this is not the frozen-slope
defect of `EOS.md` — it is the missing convention.

**Worth recording separately:** Figiel's published model is itself **6–24 % below** the
Hamer–Wu literature over 4–5.5 mol kg⁻¹ (−6.2 % at 4.0, −23.6 % at 5.5). Two errors in the
same direction is why comparing against the experimental points alone understated the gap —
below 2 mol kg⁻¹ they partly cancel. A reproduction target has to be the model, not the data
the model is fitted to.

## Addendum 3: the fix cannot live in the packet

Setting `k_ij = 1` for same-sign pairs is the right mechanism, but it **cannot be expressed
as packet data**. The builder rejects a self-pair:

```text
ParameterError: invalid component pair for pair/bromide-anion/bromide-anion/k_ij
```

which is correct — a "pair" between a component and itself is not a pair. But the like-ion
terms that matter *are* the diagonal: for a two-ion salt such as LiBr the only same-sign
combinations are Li–Li and Br–Br, and neither is addressable from a packet. The 198
`k_ij = 1.0` records in `held-cameretti-sadowski-2008/2` are all *distinct* same-sign pairs
and therefore leave every diagonal entry at 0 — that packet does not implement the exclusion
either, only two thirds of it.

So the resolution has to be code that writes the `k_ij` matrix, not a packet that declares it:
the declared model-family choice is read, and the code sets every entry with
`z_i z_j > 0` — diagonal included — to 1, with the temperature slope zeroed alongside.

## Addendum 4: borrowed, not new — the equation was carried and its decision was not

`physics/helmholtz.hpp:228-260` is **byte-identical** to
`9eccbfd:engine/src/epcsaft/_native/helmholtz.hpp:303-337`, down to the
`// EqID: same_sign_ionic_dispersion` annotation on the line computing `epsilon_ij`. The
kernel was carried verbatim and **never contained the rule**.

The rule lived one layer up, in `9eccbfd:engine/src/epcsaft/_resolve.py:735-753`:

```python
same_sign_ion_pairs = frozenset(
    frozenset((component_ids[left], component_ids[right]))
    for left, right in combinations(range(len(component_ids)), 2)
    if charges[left] * charges[right] > 0
)
k_ij = _pair_matrix(component_ids, pair_records, "k_ij",
                    required=True, not_applicable=same_sign_ion_pairs)
if ionic_dispersion == "exclude-same-sign-ion-pairs":
    # EqID: same_sign_ionic_dispersion
    size = len(component_ids)
    k_ij = tuple(
        1.0 if charges[left] * charges[right] > 0 else k_ij[left * size + right]
        for left in range(size) for right in range(size)
    )
```

A **declared model choice** — `ionic_dispersion`, default `"exclude-same-sign-ion-pairs"`,
validated against an enum at `_dictionary_io.py:750-753` — resolved into the full `k_ij`
matrix **including the diagonal** (`for left in range(size) for right in range(size)`, not
`combinations`), with same-sign pairs marked not-applicable so the loader does not demand a
record it cannot express.

Stage 0 deleted that layer. Stage 4's native parameter builder reimplemented the mapping and
not the resolution.

**So the defect is borrowed and the loss is a decision, not an equation.** This is
`TRANSFER.md`'s carry rule running in reverse: *"carry an equation, write a decision"*
assumes the two are separable, and here they were not. What looked like a self-contained
kernel was correct only in combination with a choice made in a layer that was deleted. The
`// EqID:` annotation came along and its selector did not.

The tell was in the tree the whole time. The registry entry for `same_sign_ionic_dispersion`
is a **case split** — two branches, a citation for each, and a change note saying the caller
selects between them — while the annotated line implements neither branch, only the
unconditional formula. A carried kernel whose registered equation has a case split needs its
selector carried with it, and nothing checks that.

## Addendum 5: the exclusion reproduces the published curve to 0.27 %

The superseded engine exposes the selector as an option, so both branches can be run against
the digitized Figure 5 LiBr curve on identical parameters. **Verified by the coordinator**,
27 points, 2.160–5.847 mol kg⁻¹:

| `ionic_dispersion` | AARD vs Figiel's published curve | max |
|---|---|---|
| `exclude-same-sign-ion-pairs` | **0.27 %** | 0.56 % |
| `all-ion-pairs` | **49.19 %** | 61.25 % |

| m | Figiel | exclusion | all pairs |
|---|---|---|---|
| 2.160 | 1.161 | 1.158 (−0.2 %) | 0.782 (−32.6 %) |
| 3.973 | 1.774 | 1.777 (+0.2 %) | 0.898 (−49.4 %) |
| 5.572 | 2.548 | 2.536 (−0.5 %) | 1.021 (−59.9 %) |

Three things at once, from one experiment:

1. **The diagnosis is exact.** The greenfield engine's 49.19 % is reproduced to four decimals
   by the superseded engine forced onto `all-ion-pairs`. It is running the wrong branch, and
   nothing else is wrong with it on this path.
2. **The fix is confirmed before it is written.** The exclusion lands inside the
   digitization uncertainty (0.45–0.8 % per point) across the whole range.
3. **The digitized curve is validated as a target.** A tracing off a page that agrees with an
   independent implementation to 0.27 % is measuring the published model, not the tracer.

