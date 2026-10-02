# Charged-trial domain addendum readiness

Reviewer: independent GPT-6.1 Sol, xhigh; requested 2026-10-01.
Task: `0881f4ac-fdf6-4480-b7a1-8a7193284ac5`.
This is design readiness, not delivered capability review.

**Verdict: Block for resuming reacting stability/discovery qualification.** The addendum’s conclusion and stop boundary are supported. No computable, source-backed charged-liquid admission criterion is ready to implement from the inspected sources, existing density rules, or selected MEA v3 inputs. Accurate retention of that unresolved question does not establish readiness.

This verdict covers only the charged-trial domain addendum in [#223](https://github.com/tannerpolley/ePC-SAFT/issues/223). I read the full current body and all six comments, and confirmed that the body remained unchanged at the end of inspection.

1. **Existing density admission accepts this point mechanically.** `ALG-DENSITY` leaves a genuine sole stable root unlabelled. The declared-phase check selects the cold root, anchors on its density, and checks pressure slope, closure and displacement. It imposes no additional electrolyte applicability boundary. The retained point has \(\rho=2.321466537657506\ \mathrm{mol/m^3}\), positive \(dP/d\rho=2504.2517\ \mathrm{Pa/(mol/m^3)}\), and the same sole `Unspecified` root for cold Liquid and Vapor requests. See [ALG-DENSITY](https://github.com/tannerpolley/ePC-SAFT/blob/443d8289/engine/docs/science/algorithms.md#alg-density), [density.cpp, lines 382–406](https://github.com/tannerpolley/ePC-SAFT/blob/443d8289/engine/native/properties/density.cpp#L382-L406), and [shared admission, lines 582–629](https://github.com/tannerpolley/ePC-SAFT/blob/443d8289/engine/native/equilibrium/phase_boundary.cpp#L582-L629).

   Testing the anchored `Liquid` label would not resolve this: the anchor routine returns the requested label after establishing local uniqueness. Requiring the **unanchored** root to carry `Liquid` would exclude sole dense roots as well. That would define a narrower numerical search domain without deriving a physical validity condition. The [formulation, lines 137–140](https://github.com/tannerpolley/ePC-SAFT/blob/443d8289/engine/docs/equilibrium/formulation.md#L137-L140) explicitly disclaims universal phase classification by the selector.

2. **The inspected sources explain the limitation but supply no exclusion inequality.** Bülow–Ascani–Held Part I, §2.1, Eqs. (2)–(4), neglects dielectric/Born-radius pressure dependence and assigns zero Born pressure contribution. Ascani–Held’s thermodynamic framework, pp. 3–4, Eqs. (10)–(13), uses continuum solvation and composition mixing for liquid–liquid applications. These passages do not derive a density or branch boundary applicable to this MEA trial. [Bülow–Ascani–Held 2021](https://doi.org/10.1016/j.fluid.2021.112967), [Ascani–Held 2021](https://d-nb.info/1232858226/34).

   Perdomo2025 §4.5, Eqs. (73)–(76), uses density-dependent solvent permittivity in its SAFT-γ Mie implementation. It supplies an alternative formulation, not an admission condition for the selected solvent-only MEA model. The existing [full-EOS diagnosis, lines 514–521](https://github.com/tannerpolley/ePC-SAFT/blob/443d8289/engine/docs/science/full-eos.md#L514-L521) likewise withholds authority for silent dielectric substitution. [Perdomo2025](https://doi.org/10.1016/j.compchemeng.2024.108977).

3. **The negative calculation remains consequential evidence.** At 318.15 K and 5918.5286 Pa, the retained direct tangent-plane distance is −4.291232943714732 and the donor slope is −4.291232943717396. The finite split gives
   \[
   \Delta G/RT=-1.9951025594702504\times10^{-5}\ \mathrm{mol},
   \]
   versus the native estimate \(2.7239525076341137\times10^{-7}\ \mathrm{mol}\), about a 73-fold margin. Donor pressure reclosure preserves −1.99510237757131×10⁻⁵ mol; smaller splits preserve descent. All nine inventory errors at the original split are zero in the retained arithmetic, with charge residuals below 4.8×10⁻¹⁷. See the [retained diagnosis](~/.t3/worktrees/ePC-SAFT/t3code-223-build-tpd-stability-test-and-phase-35ebff8f/analyses/2026-tpd-phase-discovery/results/diagnostics/phase-domain/diagnosis.md:14).

   The optimizer residual of 6.28×10⁻⁶ prevents a stationary-trial identity claim. It does not invalidate the separately reconstructed feasible descent permitted by adopted D9. The solvent-only model’s \(\varepsilon_r=71.4898\) and Born contribution −322.0579 RT at this density explain why physical applicability needs resolution. These are numerical evaluations with estimated uncertainty; physical validation, rigorous sign bounds and a converged added phase remain unestablished.

The **smallest required scientific decision** belongs to the owner: establish and adopt a computable physical-domain condition with primary-source applicability to these inputs, or authorize a dielectric/model-input redesign if such a condition cannot be established. The MEA packet’s [current domains](https://github.com/tannerpolley/ePC-SAFT/blob/443d8289/data/packets/mea-co2-h2o-nine-species-estimation/3/packet.toml#L49-L63) provide T/P ranges and candidate-extrapolation qualifications, without a physical density boundary. A cutoff, failed-trial filter or changed reference expectation cannot fill that gap.

Before resumption, refresh readiness against the selected condition or redesign and its numerical falsifier. Subsequent qualification must preserve this exact miss and refresh the existing case-7 below/above-bubble and pinned-snapshot checks on the identified wheel, with the original inputs and tolerances unless the owner explicitly changes them. Source checks, retention, neutral/case-6 qualification and acceptance item 10 remain authorized independent work.

Inspection used branch `issue-223-phase-discovery`, HEAD `443d8289`, its working tree, and retained diagnosis wheel `384d608…`. I verified the original diagnostic hashes and copied probe source. No builds, installations, Engine calls, tests, file edits or Git/GitHub mutations were performed. Delivered review remains separate and pending.
