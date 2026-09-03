# The dielectric suppression coefficient, and what the activity gate measures

Status: findings, 2026-09-02. Written because Stage 4's activity gate passed while the
engine's mean ionic activity coefficients drifted low at high molality, and the first
account of that drift — recorded in `STAGE-LOG.md` and corrected below — was wrong in two
places.

Evidence is **measured** (run on this tree today), **derived** (follows from the equations),
or **guessed**.

## 1. The equation, and where the number comes from

Every ePC-SAFT electrolyte term reads a bulk relative permittivity that falls as salt is
added. The engine spells that with one equation (`EqID: ion_suppressed_bulk_permittivity`,
`physics/electrolyte.hpp`):

```text
ε = ε_solvent-mixture / (1 + Σ_i α_i x_i)
```

`α_i` is `ion_suppression_coefficient[i]`, and the sum runs over ions in ion mole fraction.
This is Zuber et al. 2014 Eq. 9 (Fluid Phase Equilibria 376, 116–123). **One fitted `α` per
ion per solvent** — Table 2 gives eighteen for water, Table 3 eight for methanol, Table 4
four for ethanol, and Table 9 a general cation and anion value per solvent for ions the
paper does not fit. Now in the tree as
`data/packets/zuber-2014-electrolyte-permittivity/1`.

Two published values are negative (Ca²⁺ −0.53, SO₄²⁻ −0.66 in water) and are retained as
fitted: the paper explains that the dielectric constant of water *rises* with MgSO₄ to a
point before falling, which Eq. 9 can only represent with a negative parameter.

**The engine has always supported per-ion values.** `parameter_builder.cpp` accepts the
`ion-specific-suppression` choice and requires a coefficient for every charged component
(`:374`); the `ion-fraction-suppression` choice is the one-number form, read from
`model_coefficients` and broadcast to every ion (`:366-370`). What was missing was the data,
not the capability.

## 2. Why Figiel uses 7.01 for every ion, and why that is not a free parameter

`figiel-2025-reference-electrolytes/2` declares `ion-fraction-suppression` with a single
`7.01`. That is not an arbitrary placeholder to be swapped. **measured**, from
`analyses/reproduction/2025-figiel/scripts/check_figiel_2025_source_packet.py:275-293` and
the parameter-provenance ledger it checks: the Born diameters of Li⁺, Na⁺, K⁺, Cl⁻ and Br⁻
are classified `fitted-in-Figiel-2025` and are staged **after** the dielectric suppression
stage, whose single reported value is `7.01`.

So the ion diameters in that packet were regressed **with `α = 7.01` held fixed**. The two
are one calibration. Changing `α` alone leaves diameters that were fitted to absorb its
error.

> **The rule.** `born_diameter` and the dielectric suppression coefficient are **paired**.
> A packet's Born diameters are only valid with the suppression coefficient they were
> regressed against, and this packet's are only valid with `7.01`. Never evaluate
> `figiel-2025-reference-electrolytes` with any other suppression value, and never
> substitute a suppression coefficient into a packet without re-fitting the diameters in
> the same step. Any consumer that changes one must change both, and a fit that varies one
> while holding the other must say so.

## 3. The experiment

Per-ion Zuber values substituted into the Figiel packet, nothing else changed, mean ionic
activity coefficient against the digitized literature of Figiel Figure 5, `m ≤ 2 mol/kg`,
21 points per salt (**measured**, `analyses/reproduction/2025-figiel/figures/figure-05/`):

| salt | AARD, α = 7.01 | AARD, Zuber per-ion |
|---|---|---|
| LiCl | 6.23 % | 7.02 % |
| LiBr | 6.83 % | 6.97 % |
| KCl | 5.29 % | 9.11 % |
| KBr | 7.54 % | 10.62 % |
| NaCl | 3.34 % | 6.89 % |
| NaBr | 2.18 % | 4.58 % |

**Physically correct per-ion coefficients make the fit worse on every salt.** That is the
expected outcome of §2 and is evidence for it, not against per-ion support: the diameters
are absorbing the difference. Using Zuber's values properly means re-regressing the Born
diameters against them, which is Stage 10's kind of work, on Stage 10's machinery.

## 4. The mechanism, and the direction of the error

A mean ionic activity coefficient falls steeply from 1 at low molality — Debye–Hückel — then
passes through a minimum and rises. The rise is what the Born term does when the permittivity
itself depends on composition. The Born contribution to a species chemical potential goes as
`−(z_i² e²/8πε₀ d_i)(1 − 1/ε)`; as salt is added, `ε` falls, `(1 − 1/ε)` falls, the Born
energy becomes less negative, and `γ±` is pushed back up. **Stronger suppression means a
stronger upturn.** `derived`, and confirmed by the experiment: lowering K⁺'s coefficient from
7.01 to Zuber's 2.43 weakens the upturn and moves KCl *further* below the data, not toward it.

Against the digitized literature, every deviation is one-signed and grows with molality
(**measured**, α = 7.01):

| salt | 0.001 m | 0.1 m | 0.5 m | 1.0 m | 2.0 m |
|---|---|---|---|---|---|
| KCl | −0.0 % | −0.5 % | −3.0 % | −7.8 % | −19.3 % |
| KBr | −0.0 % | −1.0 % | −5.2 % | −11.4 % | −24.9 % |
| LiBr | −0.0 % | −1.5 % | −4.6 % | −9.7 % | −23.6 % |
| NaCl | −0.0 % | −0.9 % | −1.9 % | −4.5 % | −12.4 % |

Exact at infinite dilution, and the upturn is too weak. Whether the cause is the Born
diameters, the suppression, or something else is a question for the planner; §5 says what has
already been ruled out.

## 5. Ruled out: the mapping

Every field this model reads is mapped and assigned by Stage 4's builder, checked line by
line: `ion_suppression_coefficient` per component (`:227`) and from the model-level record
(`:366-370`), `ionic_region_relative_permittivity` (`:292`), `born_c_shell` and
`born_c_dielectric` (`:297-298`), each rejected by name when the declared choice requires it
and it is absent (`:359-363`). The physics layer reads all four (`physics/born.hpp:42-76`).
The drift is not a transcription loss between packet and array.

## 6. Correcting the first account

`STAGE-LOG.md` first recorded this as the engine being 1.8–2.9× worse than the superseded one
and missing a qualitative feature. Both halves were wrong, for one reason: that comparison
used the six oracle rows against
`analyses/reproduction/2025-figiel/results/data/figiel-2025-miac.csv`, and **that file's rows
carry a different `parameter_fingerprint` per salt** — it is a per-salt fitted comparison, not
one model. Against the one consistent superseded run in
`figures/figure-05/results/data/figure-05-series.csv`, over 21 points per salt:

| salt | superseded | this engine |
|---|---|---|
| LiCl | 6.09 % | 6.23 % |
| LiBr | 7.10 % | 6.83 % |
| KCl | 6.10 % | 5.29 % |
| KBr | 5.26 % | 7.54 % |
| NaCl | 6.58 % | **3.34 %** |
| NaBr | 8.35 % | **2.18 %** |

Comparable or better on four of six. And at low molality the superseded engine returns
**exactly 1.0000** for NaCl and NaBr below 0.05 mol/kg — no electrostatics at all, up to 14.7 %
wrong — where this engine matches the literature to 0.02–0.34 %. The "missing upturn" claim
was also wrong: LiCl and LiBr are the only salts with an upturn inside 2 mol/kg and this
engine reproduces both, too weakly. KCl and KBr fall monotonically in the *data* as well;
their minimum is near 2.5–3 mol/kg, outside the oracle's range.

## 7. What the gate measures, and what it should

The gate asks whether the engine reproduces **the measured activity coefficients** within the
parameter source's own adequacy criterion. That is a fact about the world and a legitimate
oracle claim. It is not the question the reproduction is trying to answer, which is whether
the engine reproduces **the model Figiel published**.

Nothing in the tree currently answers the second. The digitized files under
`figures/figure-05/source/` are literature points; `figure-05-series.csv`'s model rows are the
superseded engine's own output, and `inputs/docs/legacy-derived-gates.json` keeps only scalar
per-figure bounds derived from it (figure 5: RMSE ≤ 0.151, max ≤ 1.276) with the model series
deliberately not retained. **Missing fact, named:** Figiel's own reported deviations on these
salts, per salt, from the paper. With them, the reproduction claim becomes "our AARD is within
a stated margin of the published AARD" — a comparison against the paper rather than against a
superseded implementation.

Two things follow for the briefs, neither of which is a change to Stage 4:

1. An average over a molality range cannot see a systematic drift that is exact at one end.
   A criterion of the form "AARD ≤ x" should be paired with a **shape** condition — the sign
   of the deviation must not be constant across the range, or the highest-molality row must
   be within its own bound.
2. Per-ion suppression is data plus a re-regression, not a substitution. The packet is in the
   tree; the fit that would make it usable belongs to Stage 10.
