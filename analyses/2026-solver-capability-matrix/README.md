# Solver capability matrix: the ePC-SAFT Engine and comparable tools

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. This is a survey
of pinned software sources and retained numerical evidence, rather than one current-main
numerical producer. Its dated native derivative inventory requires the exact release
checkout.

**In the pinned sources surveyed, no other tool was identified with this built-in combination: the ePC-SAFT equation of state, chemical-reaction and phase equilibrium of charged species solved together in one calculation, and exact derivatives of the solved state.** The Engine's derivative orders and restrictions are stated separately in row 17: first-order feed responses on the reactive extraction case, and orders one to three on the cited neutral, ionic single-phase and reacting cases. Among the SAFT-family tools checked (Clapeyron.jl, FeOs, teqp, thermopack, pcsaft, the CoolProp PC-SAFT backend, gSAFT and Multiflash), the Engine is the only one shown to have chemical-reaction equilibrium; the gSAFT and Multiflash cells are unknown, not absent. Reactive electrolyte liquid–liquid equilibrium as such is not unique to the Engine. Aspen Plus and OLI document it with activity-coefficient models. Reaktoro has a documented example of aqueous speciation beside a separate liquid. The Held and Sadowski group has published simultaneous ePC-SAFT phase and dissociation equilibria that include liquid–liquid systems, and Tsanas et al. a reactive eNRTL flash with ions in both liquids (see the related-work note). On small neutral equilibrium problems, the specialized solvers in FeOs, Clapeyron, teqp and thermopack run 1.7 to 130 times faster than the Engine's single Gibbs-minimization driver in the retained measurements. FeOs is slower than the Engine on one water/pentanol LLE case. On electrolyte LLE with constant permittivity, Clapeyron's electrolyte flash finds the same states as the Engine and runs about 6 to 140 times faster.

The matrix has 20 problem classes and 15 tools. One row per (class, tool) cell is in [results/capability-matrix.csv](results/capability-matrix.csv), giving the method, evidence type, exact locator and notes. Sources and access dates are in [results/sources.md](results/sources.md). The eight published-case reproductions of #387 are folded into the Engine rows; they are summarized, with the ion-diameter convention and the capability gaps they found, under [published-case reproductions](#published-case-reproductions-387).

| Problem class | Engine | Clap. | FeOs | teqp | tpack | pcsaft | CoolP. | Reakt. | PHREEQC | Cantera | IDAES | Aspen† | OLI† | gSAFT† | Mflash† |
|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| pure saturation and density roots | G | S* | S* | S* | S* | S | S | ? | ? | S | P | S | ? | S | S |
| neutral multicomponent PT flash | G | G | S* | P* | S* | ? | P | P | ? | ? | S | S | ? | S | S |
| bubble/dew points | G | G* | S* | S | S* | S | S | ? | ? | ? | P | S | ? | S | S |
| phase envelopes and critical points | G | S | S* | S* | S | ? | ? | ? | ? | ? | ? | S | ? | S | S |
| stability analysis (tangent-plane distance) | G | S | S | P | S | ? | ? | ? | ? | ? | ? | S | ? | ? | ? |
| neutral LLE and VLLE | G | S* | S* | S | P | ? | ? | ? | ? | ? | P | S | ? | ? | S |
| associating mixtures (water/alcohol) | G | S* | S* | S | S | S | S | ? | ? | ? | ? | S | ? | S | S |
| electrolyte activity and osmotic coefficients | S | S* | U* | ? | X | P | U | S | S | P | S | S | S | ? | ? |
| Born solvation and concentration-dependent permittivity | S | S* | S* | ? | X | P | X | P | ? | ? | ? | P | ? | ? | ? |
| electrolyte VLE with chemical speciation (amine-CO2-water) | G | X | ? | ? | X | ? | ? | U | U | ? | P | S | S | ? | ? |
| electrolyte LLE (ions partitioning between two liquids) | G | S* | ? | ? | X | ? | ? | ? | ? | ? | ? | S | S | ? | ? |
| liquid-phase chemical-reaction equilibrium (speciation only) | G | ? | ? | ? | ? | ? | ? | G | S | G | S | S | S | ? | ? |
| simultaneous reactive and phase equilibrium (reactive flash) | G | X | ? | ? | ? | ? | ? | G | S | G | U | G | S | ? | ? |
| reactive electrolyte LLE | G‡ | X | ? | ? | X | ? | ? | P | ? | ? | ? | S | S | ? | ? |
| trace and zero-inventory species | G | P | ? | ? | ? | ? | ? | ? | ? | P | ? | ? | ? | ? | ? |
| caloric properties and PH/PS flash | G | G | S* | G | S | P | P | G | ? | G | S | S | ? | ? | ? |
| exact derivatives of equilibrium results (state and parameters) | G | P* | P* | P | P | ? | ? | S | ? | ? | P | ? | ? | ? | ? |
| equation-oriented process coupling (IDAES/Pyomo) | P | ? | ? | ? | ? | ? | ? | ? | ? | ? | G | P | P | P | ? |
| parameter regression: pure and binary neutral | S | P | P | P | ? | U | ? | ? | ? | ? | P | S | ? | ? | ? |
| parameter regression: electrolyte and reactive | S | P | ? | ? | ? | ? | ? | ? | ? | ? | P | S | ? | ? | ? |

**Legend.**
- **G** native-general: one general formulation solves this class and is shared with other specifications, for example a Gibbs minimization or a generalized flash used for several specification pairs. For the Engine, G marks equilibrium specifications solved by its one driver. Its property evaluations, density roots and fitting routines are separate and are marked S.
- **S** native-specialized: a dedicated routine for that class.
- **P** partial: a material scope limit (see notes).
- **U** user-assembled: the user builds the calculation from lower-level outputs.
- **X** not-supported: the source or documentation excludes it.
- **?** not-found: the evidence checked does not show it. This is unknown, not proof of absence.
- **\*** a retained measurement compares this tool with the Engine on the same case. A tool can offer a general path and still have been timed on a specialized one; the CSV notes name the timed path (for example Clapeyron's `bubble_pressure` for C2).
- **‡** scope restriction in the CSV notes: the Engine witnesses for row 14 have ions confined to the aqueous liquid, or only trace ions in the organic liquid.
- **†** commercial; documentation only.

Tool columns: Clap. = Clapeyron.jl; tpack = thermopack; pcsaft = the Zmeri `pcsaft` Python package; CoolP. = the CoolProp PC-SAFT backend; Reakt. = Reaktoro; IDAES = the IDAES generic property framework; gSAFT = the gSAFT property package in gPROMS (row 18 refers to the gPROMS host instead); Mflash = KBC Multiflash.

Abbreviations: EOS = equation of state; LLE, VLE, VLLE = liquid–liquid, vapor–liquid and vapor–liquid–liquid equilibrium; PT, PH, PS flash = phase split at specified pressure with temperature, enthalpy or entropy; ePC-SAFT = electrolyte perturbed-chain statistical associating fluid theory; MSE = OLI's mixed-solvent electrolyte model; HKF = Helgeson–Kirkham–Flowers aqueous model; eNRTL = electrolyte non-random two-liquid model; IDAES = Institute for the Design of Advanced Energy Systems process modeling framework; Pyomo = Python optimization modeling language; ASL = AMPL Solver Library interface; Ipopt = interior-point optimizer; CppAD = C++ algorithmic differentiation; MEA = monoethanolamine; KKT = Karush–Kuhn–Tucker optimality conditions.

## How each cell was decided

- **Engine.** Only from leaf issues and pull requests, `engine/tests`, `engine/native` and retained analyses, as the root `AGENTS.md` requires. The cited analyses are numerical verification and agreement between implementations, not physical validation, unless a note says otherwise. The published-case reproductions add numerical verification against the papers' model values and physical comparison with literature parameters; most of those data were used to fit the parameters, so the comparisons are not independent validation.
- **Clapeyron.jl, FeOs, teqp, thermopack.** Source read in the local checkouts under `/home/tnnrpolley21/Workspaces/Engineering/Forks-Clones/`. The local FeOs (v0.9.5) and Clapeyron (`6b8a6947`, v0.6.20 plus 100 commits) checkouts are older than the versions measured in the cross-check (FeOs 0.10.1, Clapeyron 0.6.29). Cells that changed between them use the measured version: FeOs v0.10.1 source on GitHub for derivatives, PH/PS flash and fitting datasets; Clapeyron v0.6.29 documentation for implicit derivatives and the generalized flash.
- **CoolProp PC-SAFT backend.** Source read at `aa7c6079` (the last commit touching `PCSAFTBackend.cpp`, 2026-09-27; latest release v8.0.0) and its documentation page.
- **Other tools.** Public documentation and source on GitHub, accessed 2026-10-08. The Aspen Plus cells rest on the version 10.1 manuals found online; later versions may add capability (for example ENRTL-RK amine models, which no retrieved document confirmed).
- **gSAFT and gPROMS.** The gSAFT cells rest on a 2012 abstract. The current Siemens page lists electrolytes among the systems gPROMS Properties handles, with SAFT-γ Mie, but does not say which model handles them or whether the property package includes reactions. The 2012 abstract therefore cannot exclude electrolyte or reactive cells; they stay **?**.
- **G or S.** A cell is G when the tool's documented or source-checked general formulation covers that class, even if a specialized path was the one timed. Clapeyron's generalized flash (non-reactive; PT, PH, PS, UV, QT, QP, TS and VT) makes rows 2, 3 and 16 G, but not LLE, where its use is not documented. teqp's generalized phase equilibrium makes row 16 G. Cantera's `equilibrate` (Gibbs minimization; multiphase TP, HP and SP) and Reaktoro's equilibrium solver make rows 12, 13 and 16 G. Aspen's RGibbs makes row 13 G.
- **Added tool: IDAES generic properties.** A reviewer will ask about it because the Engine claims IDAES coupling and IDAES has its own equation-oriented eNRTL with inherent reactions. Its cells are mostly inferred and should be strengthened before publication.
- **Split row.** Parameter regression is split into neutral pure/binary and electrolyte/reactive because the tools differ sharply between the two.
- **Not added.** GEMS/GEMS3K, a Gibbs-minimization geochemistry code: its solution-model catalogue (TSolMod) lists cubic EOS and aqueous activity models, but no SAFT model and no two-electrolyte-liquid example were found in a short search. Its cells would need the full catalogue (Wagner et al. 2012) to decide.

## Related work: published reactive electrolyte LLE

Ascani, Pabsch, Klinksiek, Gajardo-Parra, Sadowski and Held (Chem. Commun. 2022, 58, 8436–8439, https://doi.org/10.1039/d2cc02943j) solve phase equilibria and reaction (dissociation) equilibria simultaneously with ePC-SAFT advanced. Their examples are CO2 + water + salt (vapor–liquid) and carboxylic acid + water + toluene or methyl isobutyl ketone (liquid–liquid), from which they predict pH. The main text does not describe the solution algorithm or name public software. The Engine now reproduces their acetic acid + water + toluene case from the supporting information ([#374](../2026-ascani-carboxylic-acid-lle/notebook.qmd); row 14).

So ePC-SAFT with simultaneous reaction and liquid–liquid equilibrium is published modeling, not new with the Engine. What the Engine adds, on the evidence here:
- one driver shared with every other equilibrium specification it admits;
- exact derivatives of the solved state, with the orders and restrictions in row 17.

Tsanas, de Hemptinne and Mougin (Chem. Eng. Sci. 2022, 248, 117174, https://doi.org/10.1016/j.ces.2021.117174) present a reactive flash for multiple ion-containing phases: Gibbs minimization with the electric-potential contribution, per-phase electroneutrality and stability analysis on the electrochemical potential. They apply it to water/DEEA/MAPA + CO2 with eNRTL in both liquids (abstract read only). So reactive LLE with ions in both liquids is also published with an activity-coefficient model.

Whether the Engine is public software is not settled: it is GPL-3.0 licensed, but the repository was private on 2026-10-08. Their code was not available for comparison.

## Per-class notes on method and speed

**Pure saturation, PT flash, bubble/dew, envelopes and critical points (rows 1–4).** Every SAFT-family tool has dedicated Newton or successive-substitution solvers here, and Clapeyron also has a generalized flash. The Engine solves the same states as specifications of one Gibbs-minimization driver (Ipopt with exact CppAD Hessians). Values agree to about 1e-9 relative or better; the Engine is 7.7 to 130 times slower on these cases ([retained measurements](#retained-measurements)). That cost is the price of one driver. A 2026-09-23 trace attributes 2.85 ms of a 6.15 ms flash to 333 EOS evaluations and 3.17 ms to driver time (Ipopt, MUMPS, preparation and checks), so EOS speed alone cannot close the gap. The Engine's envelopes cover neutral vapor–liquid only; Clapeyron, teqp and thermopack also cover three-phase or liquid–liquid features. CoolProp's PC-SAFT backend has saturation and bubble/dew (quality) flashes; its PT input classifies the phase from bubble and dew pressures and does not split the feed.

**Stability and neutral LLE/VLLE (rows 5–6).** Clapeyron, FeOs and thermopack have tangent-plane stability searches. The Engine's explicit search covers neutral, reacting and charged trials (charged only with a permittivity form that reaches vacuum permittivity at zero density) and supplies no global certificate. On water/pentanol LLE the Engine is faster than FeOs at 298.15 K and slower at 323.15 K; Clapeyron is about 18 to 19 times faster than the Engine at both.

**Association and electrolyte properties (rows 7–9).** Association is standard in all SAFT-family tools. For NaCl activity, the Engine, FeOs and Clapeyron agree within 1.18e-6 in ln gamma under matched conventions. The Engine has a Born term and four permittivity equations: ion suppression (shared by three input labels), component mixing, density–temperature and dipolar perturbation. Clapeyron and FeOs have Born terms and alternative permittivity models. At fixed NaCl–water states the direct Born term agrees with both to within 7e-16 once the codes' electrostatic-constant ratios are removed, and the Rueben perturbation form agrees with FeOs to 8.2e-12 ([head-to-head](../2026-head-to-head-comparisons/notebook.qmd)). Of the built-in Clapeyron and FeOs models inspected, neither offers the corrected Born term or the ion-suppressed permittivity; the Engine gives both. CoolProp's ion term is Debye–Hückel with the dielectric constant of pure water and no Born term. Geochemical tools (Reaktoro, PHREEQC, Cantera) and Aspen/OLI give activity and osmotic coefficients from activity-coefficient models (Pitzer, HKF, eNRTL, MSE), not from an equation of state. Reaktoro's HKF Born term acts on aqueous solutes through the dielectric properties of pure water.

**Electrolyte phase equilibrium with speciation and reactions (rows 10–14).** This is where the matrix separates the tools.
- *Electrolyte LLE with ions in both liquids (row 11).* Clapeyron (electroneutral Rachford–Rice and Michelsen flashes), Aspen and OLI document it; the Engine solves it in the same driver. On the eight Khudaida feeds, with a matched model (plain Born and Debye–Hückel at σ) and constant permittivity, Clapeyron finds the Engine's states to 4.6e-9 in mole fraction from both its default start and the Engine's guesses. With composition-dependent permittivity the two differ by 5.3e-3, from a source-confirmed Clapeyron 0.6.29 derivative defect on the Debye–Hückel `debye_length` path (`ad.jl:381` zeroes κ's partial whenever the ionic-sum partial vanishes), reported upstream as [Clapeyron.jl #649](https://github.com/ClapeyronThermo/Clapeyron.jl/issues/649); there Clapeyron's default start returned no valid equilibrium ([head-to-head](../2026-head-to-head-comparisons/notebook.qmd)). FeOs has an ePC-SAFT model, and its flash can hold non-volatile components such as ions out of one phase, but no electroneutral two-phase solver with ions in both liquids was found.
- *Reactions.* Clapeyron documents its flashes as non-reactive. FeOs, teqp, thermopack and CoolProp contain no reaction-equilibrium code. Reaktoro, PHREEQC and Cantera solve reactive and phase equilibrium with geochemical activity models, cubic gases and declared phases. Aspen (solution chemistry with aqueous–organic LLE, RGibbs) and OLI (MSE applied to both liquids) document reactive electrolyte LLE with activity models.
- *Row 14 definition.* Reactive electrolyte LLE here means reactions and two liquid phases with ions present. The Engine witnesses have ions confined to the aqueous liquid, or only trace ions in the organic liquid (mark ‡); reactions with more than trace ions in both liquids are not qualified. On that same boundary Reaktoro is partial: a developer-provided example combines HKF aqueous speciation with a separate Peng–Robinson liquid holding no ions, and notes that vapor–liquid–liquid cases need further work. Cantera's multiphase `equilibrate` could pose the same boundary with declared phases, but no example was found, so it stays **?**. Aspen and OLI document ions in both liquids.
- *Engine evidence.* Key-solve work for the nine-species MEA reacting flash and the Khudaida ionic LLE point is under [retained measurements](#retained-measurements). The Li/Na extraction case (two liquids, two reactions, charge balance) closes moieties to 1.42e-14 mol, but its parameters are calibrated pseudo-components with no physical validation. The Ascani 2022 acetic acid + water + toluene case (#374) solves acid dissociation and water autoprotolysis with the liquid split, with reaction log-affinity at most 1.1e-11 over 805 solves. With ions aqueous only, 7 of 19 aqueous pH points lie within 0.005 of their ePC-SAFT advanced values, and differences reach −0.017 at high acid, cause unresolved. Against experimental pH the AARD is 6.95 % for the Engine and 6.90 % for theirs (physical comparison with literature parameters). With ions permitted in both liquids the solve converges, but the returned organic ions are trace (x(Ac⁻) 5.9e-22 to 6.7e-20, x(H⁺) 2.8e-36 to 1.2e-34) and their charge split is undetermined, so single-ion quantities in that liquid are not results (#375, ALG-EQUILIBRIUM) ([notebook](../2026-ascani-carboxylic-acid-lle/notebook.qmd)).
- *Speed.* On row 11 with constant permittivity, Clapeyron takes 2.1–5.0 ms per warm solve against the Engine's 30.6–288 ms, about 6 to 140 times faster; one process per code, 10 samples, not the three-process protocol. No competitor has been timed on rows 10 and 12–14.

**Trace species (row 15).** The Engine preserves certified structural zeros and resolves the cited positive trace cases within stated tolerances (exact support faces, roots below the optimizer composition box). Clapeyron removes zero components before solving. Cantera's VCS solver tracks minor, very small and zeroed species, deletes species below cutoffs and adds them back at the end; no trace qualification was found for it.

**Calorics and PH/PS (row 16).** Common to most tools. The Engine's neutral PH/PS states agree with FeOs's flashes on the same inputs; its reactive calorics (heat of absorption, apparent heat capacity at prescribed T and P) have no checked open-source counterpart with an equation of state. UV and charged PH/PS remain absent in the Engine. CoolProp's PC-SAFT backend rejects PH and PS inputs.

**Derivatives, process coupling and fitting (rows 17–20).**
- *Engine derivatives: order.* Directional actions of solved states of order one to three.
- *Engine derivatives: directions and qualified cases.*
  - Neutral vapor–liquid ternary: orders one to three along T, P, feed and two EOS-parameter directions (k_ij and segment number).
  - Single-phase water/Na⁺/Cl⁻: orders one to three along T, P and feed.
  - Reacting isomerization: orders one to three in T against analytic values.
  - Reacting caloric actions on manufactured branches (#225).
  - Reactive extraction (two liquids, two reactions, ions aqueous): first-order feed responses within 3.5e-6 of finite differences.
- *Engine derivatives: restrictions.* With source-referenced reactions, third order admits only feed directions at fixed T and P, and second order admits at most one parameter or reaction-law direction (`engine/native/equilibrium/implicit_sensitivity.cpp:332,344,347`, #225). At critical candidates only first actions are available; second and third are excluded (#224). No action crosses phase appearance.
- *Engine derivatives: cost.* Three batched exact first responses cost 1.642 ms against 243.314 ms for six finite-difference solves (148 times less), on projection wheel `aa95dcf7…`, run 2 of [reactive-projection-timings.csv](../2026-open-source-equilibrium-comparison/results/reactive-projection-timings.csv). This is incremental derivative cost at a solved state, not complete solve cost. The [reactive-electrolyte-driver notebook](../2026-reactive-electrolyte-driver/notebook.qmd) holds a separate campaign (2.076 ms against 786.728 ms).
- *Agreement and cost against competitors.* On C2 bubble (T, k_ij, association energy and volume) and C4 LLE (T, k_ij), the Engine's exact first actions agree with Clapeyron's implicit derivatives to 4.0e-13 (C2) and 3.4e-10 (C4) relative, and with FeOs's C2 parameter gradients of pressure to 2.2e-10. The Engine's actions on the solved state cost about 1/215 (C2) and 1/90 (C4) of its own central differences ([head-to-head](../2026-head-to-head-comparisons/notebook.qmd)).
- *Competitor derivatives.* FeOs 0.10.1 gives first-order parameter gradients of vapor pressure, liquid density, bubble/dew pressure and a few other pure and binary properties, and is far faster than the Engine on them (0.01043 ms against 13.524 ms for C1). Its Python interface exposes no LLE derivative. Since v0.10.0 its Rust `tp_flash_binary` also differentiates both phase volumes and compositions with respect to T, P, feed or EOS parameters for binaries; the Python wheel does not expose that. Clapeyron 0.6.29 documents arbitrary-order implicit derivatives through IFTDuals (ForwardDiff only), without naming the solvers covered; first order was measured through `bubble_pressure` and `tp_flash2`. Reaktoro gives first-order sensitivities of reactive equilibrium with geochemical models.
- *Process coupling.* IDAES and the gPROMS host are equation-oriented platforms in their own right. The Engine's IDAES coupling is qualified only for neutral properties and a nonreactive two-liquid electrolyte settler.
- *Fitting.* The Engine fits with Ceres least squares and guarded minimax, separate from its equilibrium driver; predictions inside a fit are driver solves. It fits activity, osmotic, solvation, permittivity, ln K, heat-capacity and enthalpy observations with exact implicit Jacobians. Clapeyron fits with a gradient-free optimizer by default. FeOs returns cost vectors and exact gradients for a few pure-fluid and binary properties. Aspen documents electrolyte parameter fitting; its derivative method was not checked.

## Published-case reproductions (#387)

Each case runs a published ePC-SAFT or PC-SAFT paper natively with its printed parameters. "Numerical verification" compares the Engine with the paper's model values, mostly digitized from figures, so agreement is bounded by reading precision. "Physical comparison" compares the Engine with the paper's experimental data using literature parameters. No case fits parameters or changes Engine source. Wording follows each case's reviewed record; numbers are in its notebook, computed from its `results/`.

| Case (merge) | Rows | Result, reviewed wording | Evidence type |
|---|---|---|---|
| [2001-gross](../reproduction/2001-gross/notebook.qmd) (`b2227137`, #396) | 1, 3, 6 | Partial reproduction: saturation states for all 78 compounds; 50 of 53 sampled model points within tolerance (12/12 pure, 38/41 binary). Aggregate AAD tables not reproducible (original data absent). Since #391 cold propylene converges from the default start; its and diethyl ether's relative pressure closure stay outside the criterion. | numerical verification; physical comparison |
| [2002-gross](../reproduction/2002-gross/notebook.qmd) (`997bfd91`, #402) | 1, 3, 6, 7 | 209 of 223 model rows within tolerance, 12 misses undiagnosed; water/1-pentanol three-phase point 0.66 K below the trace; since #391 only three water/1-pentanol dews fail from the default start, all recovered. | numerical verification; physical comparison (fitting data, not independent) |
| [2005-cameretti](../reproduction/2005-cameretti/notebook.qmd) (`cd2b3bc7`, #394) | 3, 8 | Fig. 2 vapor pressure partly reproduced (original ion table, d = σ: 248 of 280 within tolerance); densities and Table 5 not reproduced (paper-side basis and printing issues, inferred). since #391 the default start solves every bubble point. | numerical verification; physical comparison; prediction (activity) |
| [2014-held](../reproduction/2014-held/notebook.qmd) (`220132b5`, #397) | 8, 11 | 7 of 8 osmotic series within 0.005 with a 0.88σ Debye–Hückel closest approach (inferred); Fig. 6a LLE reproduced over the traced range; Fig. 3a and Fig. 5 NaCl not reproduced; Fig. 4a and 4b capability gaps. | numerical verification; physical comparison |
| [2020-bulow](../reproduction/2020-bulow/notebook.qmd) (`c5c95418`, #395) | 8, 9 | Hydration and transfer energies and corrigendum activity curves reproduced with 0.88σ for packing, Born and Debye–Hückel (inferred). | numerical verification; physical comparison |
| [2021-bulow](../reproduction/2021-bulow/notebook.qmd) (`adfad5a0`, #399) | 8, 9 | With unified 0.88σ: all Fig. 4/S1 curves and 10 of 18 transfer bars (all five ethanol predictions) reproduced; with σ, 0 of 18. Methanol transfer energies unresolved. | numerical verification; physical comparison |
| [2022-ascani](../reproduction/2022-ascani/notebook.qmd) (`9896157c`, #392) | 11 | The Engine solves the mixed-salt split with ions in both liquids but does not reproduce the paper's printed phases, which are inconsistent with the printed feed. Experimental comparison pending (Pedraza 2004). | numerical verification (solve checks only) |
| [2023-ascani](../reproduction/2023-ascani/notebook.qmd) (`9feeb790`, #393) | 5, 13 | Reaction and liquid–liquid equilibrium solved together and natively; K_a consistent with the paper (exploratory); tie lines deviate 8.29 % (Senina's combined statistic). Two phase-discovery defects, fixed under #390 (PR #408). | numerical verification; physical comparison; prediction (system 2) |

**Ion-diameter convention (inferred).** Three Held-group papers reproduce only when the ion diameter in the packing, Born and Debye–Hückel terms is 0.88σ, although they print σ for the Born and Debye–Hückel terms and 0.88σ only for packing (Held 2014 Eq. 21):
- Held 2014: 7 of 8 osmotic series within 0.005 in Φ with 0.88σ, against series maxima of 0.016–0.019 with σ.
- Bülow 2020: Born with d = σ moves Li⁺ hydration from −555.0 to −489.3 kJ/mol (paper −555.7).
- Bülow 2021: 10 of 18 transfer bars with unified 0.88σ, 0 of 18 with σ.

The convention is inferred from the numbers, not printed; no case adopts it as an Engine default. Cameretti 2005, which prints no ionic hard-sphere diameter, comes closest with d = σ. With 0.88σ the [Held 2014 chloride verdict](../2026-held-2014-chlorides/notebook.qmd) fails for all four salts, including LiCl and NaCl, which pass with σ (#398).

**Capability gaps found.** Each is recorded in its case and not built:
- ion pairing (Bjerrum), Bülow 2021;
- solid–liquid equilibrium, Bülow 2021 and Held 2014 Fig. 4b;
- ionic-liquid solvents, Bülow 2021: `component-permittivity-mixing` requires neutral solvents; this admission restriction is removable and was not changed;
- molar-mass parameter correlations, Gross and Sadowski 2001 Appendix B and Fig. 5: they cannot be declared;
- native osmotic coefficient for multi-solvent mixtures, Held 2014 Fig. 4a: the native Φ admits one pure solvent, so amino acid + KCl + water has only a labelled Python diagnostic.

**Start robustness and phase discovery.** The reproductions found two Engine defect groups, both since closed:
- #391: default starts now solve every Cameretti bubble point and all but three Gross 2002 water/1-pentanol dews, which explicit starts or the cooling-onset operation recover ([Cameretti](../reproduction/2005-cameretti/README.md), [Gross 2002](../reproduction/2002-gross/README.md)). Cold propylene (Gross 2001) converges from the default start; its relative pressure closure and that of diethyl ether stay outside the criterion ([Gross 2001](../reproduction/2001-gross/README.md)).
- #390: finite-split discovery returns the lower-G two-liquid state on all four esterification cases and refuses rather than return one liquid silently (PR #408).

## Retained measurements

The Engine C1–C4 rows are fresh timings on main 026b3031 from the [head-to-head comparisons](../2026-head-to-head-comparisons/notebook.qmd) ([timing summary](../2026-head-to-head-comparisons/results/h2h-engine-timing-summary.csv)). The competitor columns, and the Engine methane/ethane rows, are the retained historical measurements, not rerun. Times are milliseconds, warm process, model construction excluded; they compare the exposed calls, not equal algorithms.

Each timing is retained by its source:
- C1 propane saturation 300 K, C2 ethanol/butane bubble, C3 NaCl activity (reused reference), C4 water/pentanol LLE at 298.15 and 323.15 K, and C1 and C2 with a parameter gradient: competitors in the [cross-check](../2026-competitor-cross-check/notebook.qmd), Engine in the [timing summary](../2026-head-to-head-comparisons/results/h2h-engine-timing-summary.csv).
- Methane/ethane TP flash 180 K, 1 MPa: [open-source comparison](../2026-open-source-equilibrium-comparison/README.md).
- Methane/ethane bubble 180 K: [bubble point](../2026-binary-methane-ethane-bubble-point/README.md).
- Khudaida electrolyte LLE, 8 feeds, constant permittivity: [head-to-head](../2026-head-to-head-comparisons/notebook.qmd).

Provenance and limits, by row group:
- *Cross-check rows.* AMD Ryzen 5 5500, one CPU, one thread per pool, medians of three process medians. The Engine values are the 2026-10-08 timings of main 026b3031 under the cross-check protocol, with every timed return retained and checked against the retained Engine values; the rows shown are 17–31 % below the original campaign. Process medians spread by up to 49 % of their median (C3; 34 % on the C1–C4 equilibrium blocks) at a one-minute host load of 9.2–9.6, so speed changes between Engine versions are indicative only. The competitor values are the cross-check's retained 2026-10-03 measurements. Clapeyron's cold compile and first-call work was 604–696 s per Julia process, against 0.230–0.334 s for the Engine. Clapeyron's timed paths were `saturation_pressure` (C1), `bubble_pressure` (C2) and `MichelsenTPFlash` (C4), all specialized.
- *Methane/ethane flash.* The Engine row is wheel `aa95dcf7…` (2026-09-23), 7/16/8 work. The FeOs, teqp and thermopack rows are historical and unpaired. teqp used supplied phase guesses plus the lever rule, and Clapeyron's flash timing hit the memory cap.
- *Methane/ethane bubble.* All three rows come from one historical run.
- *Khudaida electrolyte LLE.* Matched model (plain Born and Debye–Hückel at σ), ranges over the eight feeds; median of 10 warm samples in one process per code. Clapeyron's range covers both its default start and the Engine's guesses. The Engine's two slowest feeds needed 88 and 52 iterations.

Value agreement on the same cases: propane saturation across five codes within 1.10e-9 relative in pressure; methane/ethane flash within about 1.1e-9 in mole fraction; C3 NaCl ln gamma within 1.18e-6; Khudaida electrolyte LLE with Clapeyron within 4.6e-9 in mole fraction (constant permittivity); binary critical points with FeOs and teqp within 1e-8 targets. FeOs acetone (C5) misses its 1e-8 pressure target (1.21e-7).

**Engine key solves (cold work counts, iterations/evaluations/Hessian calls, main 49aac619, `tools/time-key-solves`).** MEA reacting flash at 318.15 K, 5918.5286 Pa: 19/35/20 (43/57/44 before #370). Ionic LLE (Khudaida Figure 2 point 7): 10/15/11. The last retained absolute cold times are 96.9–114.4 ms (MEA) and about 55 ms (ionic LLE) on pre-#370 wheels ([paired campaign](../2026-competitor-cross-check/README.md#paired-ionic-lle-and-mea-cost)). No competitor was run on these two cases.

## Claims this matrix supports for a manuscript

1. No equivalent built-in combination was identified in the pinned sources surveyed: the ePC-SAFT equation of state, with chemical-reaction and phase equilibrium of charged species solved simultaneously in one calculation (electrolyte VLE with speciation, reactive flash, reactive LLE with ions confined to the aqueous liquid or only trace ions in the organic liquid), and exact derivatives of the solved state (first-order feed responses on the reactive extraction case; orders one to three on the cases listed in row 17). Evidence: rows 10, 13, 14 and 17; Clapeyron flash documentation (non-reactive); FeOs, teqp, thermopack and CoolProp source. Limits: **?** cells are unknown; Reaktoro, PHREEQC, Cantera, Aspen and OLI solve reactive phase equilibrium with other models; the Held and Sadowski group has published the modeling combination without public code identified, and Tsanas et al. (2022) a reactive flash with ions in both liquids using eNRTL.
2. Among the SAFT-family tools checked, chemical-reaction equilibrium was found only in the Engine; the gSAFT and Multiflash cells are unknown. Evidence: rows 12–13.
3. Within the Engine, one Gibbs-minimization driver with exact Hessians solves every admitted equilibrium specification in the matrix (rows 1–7 and 10–16), from pure saturation and critical points to reactive electrolyte LLE, with no class-specific equilibrium solver. Its neutral results agree with four independent codes to about 1e-9. Density roots, property evaluations (rows 8–9) and fitting (rows 19–20) use separate routines. Shared formulations also exist elsewhere (Clapeyron's generalized flash, teqp, Reaktoro, Cantera, Aspen RGibbs); what is identified only in the Engine is one driver spanning neutral, charged and reacting specifications. Evidence: rows 1–7 and 10–16; retained agreement table.
4. The Engine returns exact derivatives of solved electrolyte and reacting equilibrium states. Order, direction and restrictions differ by case (row 17): first-order feed responses on the reactive extraction case; orders one to three on the ionic single-phase, reacting isomerization and neutral parameter cases; restricted higher orders with source-referenced reactions; first order only at critical candidates. Three batched exact first responses cost 148 times less than six finite-difference solves (1.642 ms against 243.314 ms, incremental derivative cost). Evidence: row 17; `reactive-projection-timings.csv` run 2.
5. None of the open-source tools checked was identified, in the pinned sources surveyed, as fitting electrolyte and reaction-constant parameters against equilibrium observations with exact implicit Jacobians, as the Engine does. Evidence: row 20. Limit: the IDAES cell is inferred and unchecked.
6. The Engine preserves certified structural zeros and resolves the cited positive trace cases within stated tolerances. Evidence: row 15.
7. The Engine is slower than specialized solvers on small neutral equilibrium problems: 1.7 to 130 times on the retained cases (faster than FeOs on one LLE case), and about 1300 times for a C1 value-plus-gradient. Evidence: retained measurements (fresh Engine timings on main 026b3031 against retained competitor timings). Limit: process medians spread by up to 49 % on a loaded host.
8. In the cases compared, the Engine's exact first-order equilibrium derivatives agree with Clapeyron's implicit derivatives and FeOs's automatic-differentiation gradients to within 3.4e-10 relative (C2 bubble along T, k_ij and association parameters; C4 LLE along T and k_ij). This is numerical verification between implementations. Evidence: row 17; [head-to-head sensitivities](../2026-head-to-head-comparisons/results/h2h-sensitivities.csv).

## Claims it does not support

- That the Engine is the only tool for reactive electrolyte LLE: Aspen Plus and OLI document it with activity-coefficient models, Reaktoro has an ion-confined example, the Held and Sadowski group has published it with ePC-SAFT, and Tsanas et al. (2022) with eNRTL and ions in both liquids.
- Reactions combined with more than trace ions in both liquids: the Engine witnesses have ions aqueous only, or organic ions at mole fractions of 1e-19 or less with an undetermined charge split (#375).
- That a shared solver formulation is unique to the Engine.
- That the Engine is public software: the repository was private on 2026-10-08.
- That the Engine is fast on electrolyte LLE: on constant-permittivity electrolyte LLE, Clapeyron is about 6 to 140 times faster than the Engine.
- Any speed claim on reactive classes: no competitor was timed on them.
- Physical validation: the Engine rows are numerical verification. The Li/Na extraction case uses calibrated pseudo-components. The published-case physical comparisons use literature parameters fitted to much of the same data.
- Global phase stability: the Engine search gives no certificate.
- That a **?** cell means absence: it means not shown in the sources checked.
- Commercial-tool rankings: those cells rest on old or partial documentation.
- Electrolyte LLE agreement with a competitor using the corrected Born term or composition-dependent permittivity: only the plain-Born, constant-permittivity variant agrees; the composition-dependent case differs by 5.3e-3 because of the Clapeyron derivative defect.
- Mixed-salt agreement with other codes: still not measured. Two reproductions solve mixed salts natively, neither against another code: Cameretti 2005 Table 5 bubble pressures (closer to experiment than the printed model column, but the printed model values are not reproduced) and Ascani 2022 water/1-butanol/NaCl/KCl LLE (solved, but the paper's printed phases are not reproduced and are inconsistent with its feed).
- Reproduction of the Held-group electrolyte papers with ion diameters as printed: they reproduce only with the inferred 0.88σ convention.

## Head-to-head measurements that would strengthen the matrix most

Items 1–4 were run in the [head-to-head comparisons](../2026-head-to-head-comparisons/notebook.qmd) (PR #378); their results are in rows 9, 11 and 17 and the retained measurements. Ordered by value to the manuscript. Costs assume the locked comparator environments of the [cross-check](../2026-competitor-cross-check/README.md) (FeOs 0.10.1, Julia 1.13.1 with Clapeyron 0.6.29) still exist in the cache; each Julia process costs about 10–12 minutes of compilation.

1. **Electrolyte LLE, Engine against Clapeyron's electrolyte flash, on the eight Khudaida feeds.** This is the only open-source competitor for row 11, so it decides the strongest contested cell, for both values and time. Cost: about 1–2 days, mostly convention matching (screening diameter, permittivity, association size factor) as in C3; runtime under an hour.
2. **Fresh Engine timing of C1–C4 on current main with returned values retained.** The current Engine timings are historical and their returns were not retained. Cost: about 2 hours with the existing scripts; competitors need not rerun.
3. **Born term and concentration-dependent permittivity at fixed state, against FeOs (advanced) and Clapeyron.** Row 9 has no numerical comparison. Cost: about half a day, reusing the C3 harness.
4. **Exact equilibrium sensitivities, Engine against Clapeyron's implicit derivatives and FeOs gradients, on C2 bubble and C4 LLE.** This tests row 17 on values and cost. Cost: about 1 day.
5. **Charged tangent-plane test on a Khudaida feed, against Clapeyron's electrolyte trial plan.** Row 5 for ions. Cost: about 1 day.
6. **Mixed-salt activity (for example NaCl + KCl), against Clapeyron.** This also qualifies the Engine's unqualified mixed-salt cell. Cost: about 1 day.

A reacting comparison against Reaktoro or Cantera would compare solvers, not models, because neither can evaluate ePC-SAFT; it is not recommended. An Aspen ENRTL-RK MEA comparison needs a license.

## Limits

- The matrix is a dated inventory (2026-10-08), not an exhaustive survey. Each cell's notes say what was checked.
- Local competitor checkouts lag the measured versions for FeOs and Clapeyron; the cells that changed use the measured versions, and other cells may still differ in current releases.
- Commercial cells rest on version 10.1 Aspen manuals, a 2004 Aspen Polymers slide set, vendor blog and training pages, a 2012 gSAFT abstract and the current Siemens gPROMS page.
- The IDAES column is mostly inferred.
- Timings come from one host and compare exposed calls rather than algorithms. The Engine C1–C4 timings are fresh but scattered up to 49 % on a loaded host; competitor timings are historical and unpaired.
- Engine capability rows are bounded by the cases named in their evidence.
- The matrix above is a source inventory; its quoted numbers come from the retained analyses it links. The higher-order comparison notebook below owns the newly retained derivative operands.

## Named-release higher-order rerun (#680)

The inventory did not declare an executable study command or retain the
second- and third-order check operands during #680. The paper placeholder was
therefore unavailable in that rerun. The
[execution limitation and wheel identity](results/main-8298c602/unavailable.json)
record that historical disposition and remain unchanged.

## Retained higher-order comparisons (#687)

The [notebook](notebook.qmd) reads the native second- and third-order operands
from [the new results folder](results/main-8298c602-higher-order/summary.json).
It retains neutral vapor–liquid ternary, aqueous electrolyte (with and without
thermal records), and selected analytic isomerization fraction/amount comparisons.
The complete native selectors passed. These are numerical verification results
on the declared smooth branches; delivered review and investigator promotion
remain pending. No Engine, native-test, Meson or native-run code changed.

The dated inventory used a separate release checkout at
`8298c6029ebb8b467929badc4a223c580ec591ab`, with at most six build jobs and
`EPCSAFT_HSL_LIBRARY`, `LD_LIBRARY_PATH` and `LD_PRELOAD` unset. Its setup and
capture recipe is withdrawn from this frozen survey.

The producer transports existing JSON records and computes discrepancies and
maxima from their operands. All thermodynamic calculations and numerical
criteria remain in the existing C++ checks. The summary retains executable and
linked-library hashes, source/build paths, compiler/options, the single build/list
command and all capture return codes. A failed repeat removes the previous
complete summary and leaves its diagnostic and full logs for inspection.
