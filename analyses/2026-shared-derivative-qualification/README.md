# Shared derivatives and density closure of selected models (#26)

**Question.** Does an assembled, explicitly selected ePC-SAFT model keep every
internal-state (association), dielectric, temperature, density, composition and
active-parameter dependency that its state and pressure-closed quantities need,
and is neutral-solvent density adequate before salt-density error is blamed on
ionic physics?

**Answer (bounded).** Yes on the tested smooth liquid branches. All four native
permittivity forms, including the density-dependent DensityTemperature and
DipolarPerturbation forms that #27 admitted, carry ∂ε_r/∂ρ and ∂ε_r/∂T through
pressure, density closure, caloric actions of orders 1–3 on both charts and
residual energy contractions of orders 3–4. Pure water from the Figiel, Held and
MEA packets is within 0.046 % of IAPWS-95 at 0.1 MPa over 278–368 K. One
admission defect was found and fixed: a neutral selection from an electrolyte
packet was rejected.

## Evidence

| Claim | Owner | Result |
|---|---|---|
| Z = 1 + ρ ∂a_res/∂ρ with the ε_r(ρ) chain, four forms | `native.properties.issue26` (`engine/tests/native/properties_p3_f9_checks.cpp`) | 4/4, rel. error ≤ 2.4e-9 on Z−1 |
| Declared liquid root persists across the ion-bearing vapor-root fold; sole roots unlabelled; closure | same | 48/48 (4 forms × 12 pressures, 2 kPa–5 MPa) |
| Residual h, s, c_p actions, orders 1–3, T,ρ,x,θ and T,P,x,θ charts, θ = (k_ij water–Na⁺, association ε/k) | same | 72/72 on the P4 criterion, error/max(1,\|ref\|) ≤ 5e-6 / 5e-5 / 5e-4; worst scaled 1.6e-9 / 3.0e-6 / 4.4e-4 (order 3: 1.6e-3 relative to \|ref\| = 0.28) |
| E3, repeated-direction E3 and E4 of the extensive residual energy in (T,V,n,θ) | `native.properties.issue26e` | 15/15; relative error ≤ 9.3e-11. Criterion 1e-6. The reference-uncertainty estimate (the difference of the Richardson (h,h/2) and (h/2,h/4) extrapolations, asymptotic and conservative for the finer one, not a bound) must be ≤ 1e-7 and is ≤ 2.1e-10 |
| Real packets: pressure-chart residual h, s, c_p orders 1–3 in T,P,x | `results/pressure-actions.csv`, floor inputs in `results/reference-ladder.csv` | 36/36 pass on 5e-6 / 5e-5 / 5e-4 × the term scale S; worst error/S 4.8e-10 / 2.3e-6 / 1.7e-4; every reference uncertainty ≤ 0.52 × its criterion |
| Real packets: pressure closure and branch at the intended liquid state | `results/pressure-closure.csv` | reduced residual ≤ 2.7e-13; ∂P/∂ρ > 0; two stable roots |
| Neutral-solvent density vs NIST SRD 69 at 0.1 MPa | `results/solvent-density.csv`, `inputs/nist-srd69-solvent-isobars.tsv` | water ≤ 4.6e-4 (−7.4e-5 at 298.15 K); methanol +0.35 % to +0.73 % |
| Neutral selection equals the ion-free electrolyte limit | `engine/tests/eos/test_state.py::test_neutral_selection_is_the_ion_free_electrolyte_limit` | P, h_res, μ_w/RT and the 1-bar liquid root agree to 1e-12 |

Native synthetic inputs: `properties_support::aqueous_electrolyte()` (2B water,
Na⁺, Cl⁻, DH + Born) at 298.15 K, x = (0.96, 0.02, 0.02). ComponentMixing uses
ε_ion = 8, DensityTemperature uses Yang2024 water d_V = 0.3777 dm³/mol and
d_T = 1403 K, and DipolarPerturbation uses the Rueben2024 PT2 records
(full EOS §6). The action ladder runs at 2 MPa, at least 1 MPa from every
vapor-root fold. The real packets are Figiel (Born, ion-fraction suppression),
Held 2008 (DH, solvent-only), Ascani (component mixing, water–1-butanol
cross-association) and MEA v4 (temperature-polynomial ε_r, MEA/MEAH⁺/carbamate,
Wolbach–Sandler MEA–water cross volume with water σ(T)).

Reproduce with `.venv/bin/python analyses/2026-shared-derivative-qualification/scripts/run.py`
(≈2 s, installed wheel). The native checks run in `bash tools/check-all.sh`.
`results/summary.json` holds the extension hash, the packet hashes, counts and
claim limits.

## Findings

- **Admission defect (fixed).** `Parameters.from_json(figiel…, components=("water",))`
  failed with `invalid ionic_dispersion`: surface admission declared the packet's
  electrolyte choices irrespective of the selection. Full EOS §12 says the
  selection fixes topology. The builder now leaves the electrolyte, dispersion and
  permittivity forms undeclared for an uncharged selection. This is exact, because DH
  and Born vanish identically without ions.
- **Vapor-root fold moves with ε_r(ρ).** At 298.15 K and x_ion = 0.04, the
  ion-bearing vapor-like root disappears between 0.5 and 1 MPa for IonSuppression
  and ComponentMixing, between 20 and 50 kPa for DensityTemperature, and between
  95 and 100 kPa for Dipolar. As ε_r → 1 at vapor density, the ionic free energy
  there becomes prohibitive. At 1 bar the requested liquid root is then the sole
  stable root and is returned unlabelled, as ALG-DENSITY specifies. A difference
  stencil straddling the fold is correctly refused by the branch-admission
  ladder. This is branch topology, not a derivative defect.
- **The pure-water baseline is sound.** The three water parameterizations deviate
  from IAPWS-95 by −3.7e-4 to +4.6e-4 at 0.1 MPa, 278–368 K. That bounds the
  pure-solvent error a salt-density comparison would carry, but without solution
  references or a sensitivity study it does not apportion a salt-density
  discrepancy. Ionic terms are one candidate, alongside solvent–ion dispersion (k_ij),
  packing diameters and partial-volume behavior remain candidates. No retained
  aqueous NaCl solution-density reference exists (the only retained salt
  densities are VCl₃ rows outside every packet), so no salt-density diagnosis is
  made here. Methanol is 0.55 % dense at 298.15 K. Ethanol has no NIST SRD 69
  reference and is not assessed.
- **Near-zero third-order actions are cancelling sums.** A directional action is
  multilinear, so it is the sum of the 3^k exact sub-direction actions obtained by
  splitting each direction into its temperature, pressure and composition parts.
  At order 3 these thermal, mechanical and compositional terms cancel: MEA h has
  |ref| = 1.0e-3 J/mol against a term scale S = Σ|terms| = 1.9e-2 J/mol, and Figiel
  c_p has 6.1e-6 against 2.2e-4 J/(mol K). A criterion relative to |ref| therefore
  asks for more accuracy than the terms themselves carry (summation condition
  S/|ref| up to 37). The criterion is tolerance × S, as for the reaction affinity
  (criterion × Σ|ν_i action_i|). S ≥ |ref| always; for non-cancelling actions S ≈ |ref|.
- **The finite-difference reference was fixed, not the tolerance.** Previously each
  node kept its root-closure miss (reduced ≈ 2.5e-13), which third-order stencils
  amplify by 1/(2h)³ to ~1e-4 J/mol; three near-zero actions were then left
  tolerance gated. Each node now takes one Newton step at fixed (T, x) to its target
  pressure, leaving a rounding-level miss, and the ladder h…8h adds the Richardson
  value (4D(h) − D(2h))/3 with its |R(h) − R(2h)| truncation estimate. The reference
  with the smallest propagated uncertainty (node miss × ∂obs/∂ρ, 8ε|obs| per node,
  truncation) is used, and a pass also requires that uncertainty to be below the
  criterion. The former gated rows now pass: Figiel c_p error/S 3.3e-5, Held c_p
  3.7e-5, MEA h 4.3e-5 (criterion 5e-4), with error/uncertainty 0.3–0.7, so the
  remaining error is at the reference's own noise level.

## Availability matrix (smooth declared branches)

Orders are directional-action orders. "E" is the energy order of the
residual-Helmholtz contraction.

| Model / permittivity | Association | Observable | Chart | Orders qualified | Evidence |
|---|---|---|---|---|---|
| Neutral PC-SAFT | none | h,u,s,c_p,c_v (ideal/residual/total) | T,ρ,x and T,P,x | 1–3 | `native.properties.p4`, `tests/eos/test_actions.py` |
| Neutral PC-SAFT | cross-assoc. binary | h,u,s,c_p,c_v with ε/k, κ directions | T,ρ,x,θ and T,P,x,θ | 1–3 | `native.properties.p4` |
| Neutral PC-SAFT | cross-assoc. binary | residual energy | T,V,n,θ | E3, E4 | `native.properties.issue26e` |
| Electrolyte, all four ε_r forms | 2B water | Z / pressure | T,ρ,x | value + ∂/∂ρ identity | `native.properties.issue26` |
| Electrolyte, all four ε_r forms | 2B water | residual h,s,c_p with k_ij, ε/k | T,ρ,x,θ and T,P,x,θ | 1–3 | `native.properties.issue26` |
| Electrolyte, all four ε_r forms | 2B water | residual energy | T,V,n,θ | E3 (both paths), E4 | `native.properties.issue26e` |
| Electrolyte, admitted real packets | 2B / cross-assoc. | residual h,s,c_p | T,P,x | 1–3 (term-scale criterion) | `results/pressure-actions.csv` |
| Electrolyte, IonSuppression | 2B water | residual h,s,c_p; total h | T,ρ,x at ρ = 1000 | 1–3 | `native.properties.p4` (supplemental) |
| Neutral, associating, reactive, electrolyte fixtures of E5 | fixture-specific | equilibrium A1–A3 | solved state | 1–3 on those fixtures only | `native.equilibrium.e5-*` (slow evidence; not re-run here) |

**Unavailable or unqualified, and kept explicit.**
- Active electrolyte records (Born and DH diameters, ε_r coefficients, C_ε, f_k) are fixed inputs, not
  ParameterFamily coordinates. They have no parameter derivatives, and none are claimed. Ion
  σ_i is a separate record from the packing, DH and Born diameters, so an active ion σ does not move them.
- Ionic insertion at zero ionic strength: the DH leading power is 3/2, so insertion
  Hessians are unbounded there (`IonicInsertion` status; full EOS §13).
- Trace and zero-support association insertion directions remain unresolved. The archived
  finite-depth harness (tag `archive/association-archive-harness-mixed`) left 28 of 2,588
  actions unresolved. That harness evaluated a superseded finite-depth candidate and
  was not replayed. The current implicit path is qualified above on interior supports only.
- The native order-3 rows pass only the scaled P4 criterion (1.6e-3 relative at worst).
- Order-3 real-packet agreement is bounded by the reference noise (error/uncertainty
  up to 0.7), not by a demonstrated Engine error; relative to |ref| the former near-zero
  rows agree to 1.1e-4 (Held c_p) to 1.2e-3 (Figiel c_p).
- Ascani order-3 h, which passed the old relative criterion, is 9.3e-4 relative to
  |ref| on the new reference (error/S 1.7e-4, S/|ref| = 5.5).
- Active-parameter directions are qualified on the synthetic native fixtures only,
  not on the real Figiel, Ascani or MEA packets.
- Anchored pressure requests are not re-qualified here. The charged and
  associating liquid sequences in `analyses/2026-fixed-pressure-anchor` are historical numerical
  verification on that record’s wheel, not a current anchored-request claim.
- Site fractions are internal states. Their third-order chain is exercised through the
  E3/E4 energy contractions; no direct X‴ check exists.
- Near a vapor-root fold, pressure-chart differences cannot be admitted. The action
  itself remains available on the liquid root. The branch scan samples 12 pressures;
  it does not prove global branch persistence.
- A neutral selection reports `bulk_relative_permittivity = 0`: no dielectric form is
  declared. The ion-bearing model at zero ion amount reports the solvent value (78.09
  for Figiel water). Their equality holds for the thermodynamic quantities tested
  (P, h_res, μ/RT, liquid root), not for this diagnostic output.

## Claim limits

Finite states and directions on smooth declared branches. Exact differentiation is
not physical accuracy. The Figiel, Held and Ascani source evidence is 298.15 K and 1 bar;
the solvent sweep beyond it is extrapolation of the parameters. Synthetic native
inputs qualify equations and chains, not fitted-model accuracy.
