# Electrolyte Reference States And MIAC Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the bounded Figiel reference-state interface, one provider Hamer--Wu anchor, and an installed-wheel six-salt MIAC validation campaign.

**Architecture:** `src/epcsaft/reference.py` is one deep module over `EPCSAFT.evaluate`; it owns dilution, reference-state equations, and immutable results without copying the EOS. The provider owns one physical anchor and source parameter records; `ePC-SAFT-validation` owns the full Hamer--Wu campaign, metrics, receipt, and plot.

**Tech Stack:** Python 3.11+, pybind11/C++20 provider kernel, Pint quantities, pytest, standard-library validation runner, Matplotlib for the retained campaign plot.

## Global Constraints

- Work directly on `main`; do not push without explicit user approval.
- Preserve `ParameterBundle -> ParameterSet -> EPCSAFT.evaluate -> EosResult` as the ordinary state lifecycle.
- Reference calculations must consume public EOS results and must not copy Helmholtz equations.
- Use Figiel 2025 Eqs. 17 and 19--22 with an explicitly converged electroneutral dilution path.
- The physical campaign is six aqueous salts at 298.15 K and 1 bar against Hamer--Wu tables, not digitized Figure 5 points.
- The provider retains exactly one NaCl physical tolerance test.
- Validation imports only the installed public provider surface.

---

### Task 1: Consolidate the Figiel source bundle

**Files:**
- Replace: `src/epcsaft/data/bundles/figiel-2025-water-sodium-chloride/1/**`
- Create: `src/epcsaft/data/bundles/figiel-2025-reference-electrolytes/1/**`
- Modify: `tests/test_parameters.py`

**Interfaces:**
- Produces: `ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)` with selectable water/methanol + one cation + one anion subsystems.

- [ ] Write a failing catalog test that selects each of LiCl, NaCl, KCl, LiBr, NaBr, and KBr in water and one NaBr methanol subsystem.
- [ ] Run `uv run pytest tests/test_parameters.py -q` and confirm failure because the consolidated bundle does not exist.
- [ ] Add only Figiel Tables 1--5 records needed by those selections, remove the displaced NaCl-only bundle, and update exact sparse-record assertions.
- [ ] Run `uv run pytest tests/test_parameters.py -q` and confirm all parameter tests pass.
- [ ] Commit provider source-record consolidation.

### Task 2: Add the reference-state module

**Files:**
- Create: `src/epcsaft/reference.py`
- Modify: `src/epcsaft/eos.py`
- Modify: `src/epcsaft/__init__.py`
- Create: `tests/test_reference_states.py`
- Modify: `docs/science/full-provider-eos.md`

**Interfaces:**
- Produces: `EPCSAFT.reference_state(*, temperature, pressure, solvent_mole_fractions, phase)`.
- Produces: immutable `ElectrolyteReference`, `IonicActivity`, and `TransferFreeEnergy` results.
- Consumes: `EPCSAFT.evaluate(...).log_fugacity_coefficient` only.

- [ ] Write a failing test for a converged water/NaCl reference and assert exact electroneutrality, finite limiting log fugacities, and convergence diagnostics.
- [ ] Run the focused test and confirm `EPCSAFT.reference_state` is absent.
- [ ] Implement deterministic formula-unit dilution and the immutable reference result with one private evaluation helper.
- [ ] Run the focused test and confirm it passes.
- [ ] Write failing tests for Figiel Eqs. 17 and 20--22, including the required mole-fraction-to-molality conversion.
- [ ] Implement `solvation_free_energy` and `at_molality` by composing existing EOS fugacity output.
- [ ] Write failing tests for water-to-methanol transfer identity, reverse-sign symmetry, and mismatched-reference rejection.
- [ ] Implement `transfer_to` using Figiel Eq. 19.
- [ ] Add failing rejection tests for multi-salt, lone-ion, charged solvent reference, missing fugacity, and invalid quantities, then implement only those checks.
- [ ] Run `uv run pytest tests/test_reference_states.py tests/test_eos.py tests/test_properties.py -q` and confirm all pass.
- [ ] Update the scientific registry and public documentation with the admitted reference equations and limits.
- [ ] Commit the provider reference-state module.

### Task 3: Add the one provider physical anchor

**Files:**
- Modify: `tests/test_reference_states.py`

**Interfaces:**
- Consumes: water/NaCl `ElectrolyteReference.at_molality`.
- Produces: one Hamer--Wu Table 16 tolerance assertion at 1 mol/kg.

- [ ] Add the NaCl test asserting `abs(gamma_pm_m - 0.657) <= 0.08` at 298.15 K and 1 bar.
- [ ] Run the focused test and confirm it fails before any numerical adjustment.
- [ ] Correct only source tracing, composition conversion, or reference-state implementation defects revealed by the test; do not fit a parameter or loosen the approved tolerance.
- [ ] Run the full provider test suite and confirm all tests pass.
- [ ] Commit the physical anchor.

### Task 4: Build the six-salt installed-wheel campaign

**Files in `../ePC-SAFT-validation`:**
- Create: `data/hamer-wu-1972-aqueous-alkali-halides.csv`
- Create: `data/hamer-wu-1972-aqueous-alkali-halides.yaml`
- Create: `campaigns/figiel_2025_miac.py`
- Modify: `README.md`
- Modify: `CONTEXT.md`

**Interfaces:**
- Consumes: installed public `epcsaft` wheel and the catalog/reference interface from Tasks 1--3.
- Produces: one JSON receipt and one prediction CSV supplied by command-line output paths.

- [ ] Copy only the six Hamer--Wu table columns `salt,molality_mol_kg,gamma_pm_m` and record DOI, table numbers, transcription, units, and source hash.
- [ ] Write the campaign first so it fails against the currently installed provider because the new public reference interface is absent.
- [ ] Build a provider wheel, install it in an isolated campaign environment, and rerun until the campaign evaluates every source row through public interfaces.
- [ ] Implement the approved per-salt and pooled metrics and acceptance gates without importing lab or provider source paths.
- [ ] Emit artifact hash, parameter fingerprints, every prediction, metric, tolerance, and limitation in the receipt.
- [ ] Update validation repository context from Slice 1-only scope to the two owned campaigns.
- [ ] Commit the campaign source and source-data record.

### Task 5: Retain and verify campaign evidence

**Files in `../ePC-SAFT-validation`:**
- Create: `results/figiel-2025-miac.json`
- Create: `results/figiel-2025-miac.csv`
- Create: `results/figiel-2025-miac.png`

**Interfaces:**
- Consumes: the immutable provider wheel and Task 4 campaign.
- Produces: retained physical acceptance evidence and plot.

- [ ] Run the campaign against the exact built wheel and retain the JSON and CSV outputs.
- [ ] Generate one plot with all six source datasets and corresponding clean-provider predictions.
- [ ] Visually inspect the plot and verify source/model labels, axes, units, and the initial decrease/later increase behavior.
- [ ] Rebuild provider wheel and sdist, test the installed artifacts, and rerun the campaign from an isolated directory with no lab checkout on `PYTHONPATH`.
- [ ] Run the cleanup hook in both repositories and inspect both worktrees.
- [ ] Commit retained validation evidence; do not push.

