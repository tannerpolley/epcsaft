# Sources

All sources were accessed on 2026-10-08. Local paths are relative to
`/home/tnnrpolley21/Workspaces/Engineering/` for `Forks-Clones/` and to the
repository root for Engine evidence.

## Engine evidence (repository, commit 92ce206a; the head-to-head and #374 analyses from main 1c209c48)

- Leaf issues and pull requests: #9, #27, #46, #47, #172, #219, #220, #223, #224, #225, #231, #243, #269, #300, #322, #350, #352, #353, #356, #358, #374, #375; PRs #370, #377, #378 (`gh issue view`, `gh pr view`).
- Native source: `engine/native/equilibrium/implicit_sensitivity.cpp` (lines 332, 344, 347: higher-order restrictions with source-referenced reactions), `engine/native/core/parameters.hpp` (line 17: the four permittivity equations), `engine/native/surface/parameter_builder.cpp` (line 188: input labels mapped to them).
- Tests: `engine/tests/eos/test_state.py`, `test_activity.py`, `test_permittivity_forms.py`; `engine/tests/equilibrium/test_problem.py`, `test_continuation.py`, `test_critical_envelope.py`, `test_caloric_specification.py`, `test_density_diagnostics.py`, `test_native_transport.py`, `test_neutral_reference.py`, `test_trace_reliability.py`, `test_declared_starts.py`; `engine/tests/regression/test_fit.py`, `test_reaction_constants.py`.
- Retained analyses: `analyses/2026-competitor-cross-check`, `2026-open-source-equilibrium-comparison`, `2026-neutral-pcsaft-equilibrium-comparison`, `2026-binary-methane-ethane-bubble-point`, `2026-binary-methane-ethane-tp-flash`, `2026-tpd-phase-discovery`, `2026-supplied-tangent-plane`, `2026-reactive-electrolyte-driver`, `2026-khudaida-idaes-settler`, `2026-idaes-asl-adapter`, `2026-figiel-divalent-chlorides`, `2026-divalent-chloride-activity`, `2026-dielectric-alternatives`, `2026-electrolyte-parameter-gradients`, `2026-reaction-constant-fitting`, `2026-reactive-caloric-actions`, `2026-pure-saturation-regression`, `2026-association-native-comparison`, `2026-equilibrium-polish-masking`, `2026-head-to-head-comparisons` (fresh Engine C1–C4 timings on main 026b3031; Khudaida electrolyte LLE, solved-state sensitivities and Born/permittivity against FeOs 0.10.1 and Clapeyron 0.6.29), `2026-ascani-carboxylic-acid-lle` (#374), `reproduction/2020-cleeton`.

## Open-source competitors, local source

| Tool | Checkout | Identity |
|---|---|---|
| Clapeyron.jl | `Forks-Clones/Clapeyron.jl` | `6b8a6947` (2026-08-07), `git describe` v0.6.20-100 |
| FeOs | `Forks-Clones/feos` | `a508463c`, tag v0.9.5 (2026-04-14) |
| teqp | `Forks-Clones/teqp` | `5f62a6f5`, v0.23.2 (2026-07-21) |
| thermopack | `Forks-Clones/thermopack` | `d68c794c` (2026-02-26) |
| CoolProp (PC-SAFT backend) | GitHub, not cloned | `aa7c6079` (2026-09-27), last commit touching `src/Backends/PCSAFT/PCSAFTBackend.cpp`; latest release v8.0.0 |

The retained cross-check measured FeOs 0.10.1 and Clapeyron 0.6.29; its FeOs locators use `https://github.com/feos-org/feos/blob/v0.10.1/`. The FeOs derivative, PH/PS and fitting-dataset cells use v0.10.1 source on GitHub (tag commit `fa37b2cd`; `crates/feos-core/src/phase_equilibria/tp_flash.rs` line 48, `crates/feos-core/src/ad/properties/`) and the installed 0.10.1 wheel in `~/.cache/epcsaft-comparators/222/feos` for the Python surface.

## Documentation, source on GitHub and papers

- Clapeyron flash methods: https://clapeyronthermo.github.io/Clapeyron.jl/dev/properties/flash/
- Clapeyron implicit differentiation: https://clapeyronthermo.github.io/Clapeyron.jl/dev/implicit_differentiation/ and, for the measured version, https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/docs/src/implicit_differentiation.md
- FeOs releases: https://github.com/feos-org/feos/releases/tag/v0.10.0 (differentiable binary flash, PH/PS flash, datasets) and https://github.com/feos-org/feos/releases/tag/v0.10.1 (fixes only)
- CoolProp PC-SAFT backend: https://coolprop.org/coolprop/PCSAFT.html; source https://github.com/CoolProp/CoolProp/blob/aa7c6079/src/Backends/PCSAFT/PCSAFTBackend.cpp
- pcsaft (Zmeri): https://pcsaft.readthedocs.io/en/latest/ and https://github.com/zmeri/PC-SAFT
- Reaktoro activity models: https://reaktoro.org/tutorials/basics/specifying-activity-models.html
- Reaktoro source (main): https://github.com/reaktoro/reaktoro (README.md; Reaktoro/Equilibrium/EquilibriumSolver.cpp, EquilibriumSpecs.hpp, EquilibriumSensitivity.hpp)
- Reaktoro HKF Born functions: https://github.com/reaktoro/reaktoro/blob/main/Reaktoro/Models/StandardThermoModels/StandardThermoModelHKF.cpp (lines 87, 97)
- Reaktoro aqueous speciation with a separate liquid phase (developer answer): https://github.com/reaktoro/reaktoro/discussions/373
- Reaktoro custom activity models: https://nbviewer.org/github/reaktoro/reaktoro-jupyter/blob/master/tutorial/eq.custom-activity-models.ipynb
- Reaktoro reactive-transport coupling: https://arxiv.org/pdf/2009.01194
- PHREEQC version 3: https://www.usgs.gov/software/phreeqc-version-3; Parkhurst and Appelo (2013), USGS Techniques and Methods 6-A43, https://doi.org/10.3133/tm6A43
- Cantera phase models: https://cantera.org/stable/reference/thermo/phase-thermo.html
- Cantera `equilibrate`: https://cantera.org/stable/python/thermo.html
- Cantera VCS species states (minor, small, zeroed, deleted): https://cantera.org/3.2/cxx/da/d80/vcs__defs_8h.html
- IDAES generic property framework (main): https://github.com/IDAES/idaes-pse/tree/main/idaes/models/properties/modular_properties (eos/enrtl.py, eos/ceos.py, base/generic_property.py, phase_equil/smooth_VLE.py, phase_equil/smooth_VLE_2.py)
- WaterTAP true and apparent species: https://watertap.readthedocs.io/en/1.0.0/how_to_guides/how_to_use_apparent_and_true_chemical_species.html
- Aspen Plus 10.1 Physical Property Methods and Models: https://sites.chemengr.ucsb.edu/~ceweb/courses/che184b/aspenplus/PhysPropMeth&Models.pdf (pages 1-4 to 1-5, 1-13, 1-16, 2-15, 3-60 to 3-62, 5-1 to 5-2, 5-9)
- Aspen Plus 10.1 Getting Started Modeling Processes with Electrolytes: https://sites.chemengr.ucsb.edu/~ceweb/courses/che184b/aspenplus/GettingStartedElectrolytes.pdf
- Aspen Polymers (POLYPCSF) slides: https://www.slideshare.net/slideshow/aspen-polymersunitopsv8-2usr/42140231
- OLI MSE and MSE-SRK: https://www.olisystems.com/resources/blog/msesrk-thermodynamic-framework-in-oli/
- OLI Engine in Aspen Plus: https://wiki.olisystems.com/wiki/OLI_Engine_in_Aspen_Plus
- MSE extraction example (RSC 2026): https://pubs.rsc.org/en/content/articlehtml/2026/cp/d5cp03992d
- gPROMS Properties (Siemens product page; lists electrolytes and SAFT-γ Mie): https://www.siemens.com/en-gb/products/gproms/process/
- gSAFT (AIChE 2012 abstract): https://proceedings.aiche.org/conferences/aiche-annual-meeting/2012/proceeding/paper/591h-gsaft-application-saft-g-mie-group-contribution-eos-oilgas-industry-academic-research
- SAFT-gamma Mie electrolytes and PSE rights (Imperial CPSE newsletter, May 2015): https://t4awstest.cc.ic.ac.uk/media/imperial-college/research-centres-and-groups/cpse/annual-report/newsletter-and-technical-articles/Technical-Article---Newsletter-Issue--12---May-2015.pdf
- Ascani, Pabsch, Klinksiek, Gajardo-Parra, Sadowski, Held, Prediction of pH in multiphase multicomponent systems with ePC-SAFT advanced, Chem. Commun. 58 (2022) 8436–8439, https://doi.org/10.1039/d2cc02943j (open-access copy https://eldorado.tu-dortmund.de/bitstream/2003/41310/1/d2cc02943j.pdf; main text read for the matrix; the supporting information was read and its acetic acid case reproduced in #374, Zotero AGPHR9MB)
- Tsanas, de Hemptinne, Mougin, Calculation of phase and chemical equilibrium for multiple ion-containing phases including stability analysis, Chem. Eng. Sci. 248 (2022) 117174, https://doi.org/10.1016/j.ces.2021.117174 (Zotero IU9YI9VZ; abstract read only)
- Clapeyron 0.6.29 `debye_length` derivative defect: `src/models/Electrolytes/electrolytes.jl:149` and `src/modules/solvers/ad.jl:381`, reproduced in the head-to-head analysis; not yet reported upstream
- GEMS3K and TSolMod: https://gems.web.psi.ch/GEMS3K/doc/html/gems3k-about.html, https://gems.web.psi.ch/TSolMod
- KBC Multiflash: https://www.colan.org/process-modeling-component/multiflash/; training outline https://www.kbc.global/uploads/files/Events/Multiflash_3_day_training_content.pdf; release 7.4 https://www.kbc.global/about/news/kbc-launches-multiflash-7.4-software-to-accelerate-decarbonization-efforts/

## Not accessible

- Aspen Plus input language guide (University of Alberta copy) returned HTTP 404.
- The RSC HTML page of Ascani et al. (2022) returned HTTP 403; the open-access repository copy was read.
- Reaktoro equilibrium-sensitivity tutorial and Cantera equilibrium reference pages returned HTTP 404; the GitHub source and the Python `ThermoPhase.equilibrate` reference page were used instead.
- Clapeyron.jl issue #649 (filed 2026-10-09): https://github.com/ClapeyronThermo/Clapeyron.jl/issues/649. Minimal reproduction on 0.6.29: ForwardDiff gives 0 for dκ/dn_water against −2.83e8 by central difference when ϵ_r depends on composition; guarding only sqrt(I) restores agreement.
