# Chemical-equilibrium validation ladder

This campaign validates the installed public
`epcsaft_equilibrium.chemical_equilibrium` route in four increasing levels of
difficulty:

1. NIST-JANAF `N2O4 <=> 2 NO2` at 500 K over three pressures. This is a
   two-species ideal-gas benchmark with a closed-form extent and a nonzero
   reaction mole-number change, so it checks the 1 bar source-to-concentration
   reference conversion.
2. NIST-JANAF water-gas shift at 600, 900, and 1100 K. This is a four-species
   ideal-gas benchmark with temperature-dependent Shomate records and a
   closed-form extent.
3. The homogeneous gas-only restriction of Belov and Aristova's 2023
   `test1.dat` carbon/oxygen system. Eight species span approximately 47
   decades in amount and are compared with an independent 80-decimal oracle.
4. IAPWS water self-ionization at 298.15 K and 1 bar through the exact
   installed Provider model and source-standard-state transform.

The first three cases validate chemical-equilibrium mathematics and numerical
robustness under an explicitly ideal-gas model. They are not ePC-SAFT
predictive evidence. The IAPWS case exercises the current nonideal
Provider-backed route, but it remains one state and one reaction. All results
remain local fixed-temperature, fixed-pressure homogeneous equilibria with
`globality_status=not_guaranteed`.

The aggregate campaign owner is
`results/chemical-equilibrium-validation.json`. The two detailed Chase result
owners are `results/n2o4-no2-equilibrium.json` and
`results/water-gas-shift-equilibrium.json`; the source-standard-state and
trace-resolution owners remain in the `2019-iapws` and `2023-belov` analysis
leaves. The single campaign entry point is the historical
`scripts/validation/chemical_equilibrium_validation_ladder.py`; it verified
the exact wheel hashes and installed `RECORD` members before importing either
production package. No campaign wrapper or replay artifact is an active owner.

The separate `campaigns/2018-iapws/cases/carbonate-complexation/` source packet is not an
installed-artifact campaign. Its carbonate and calcium-complexation cases
require species-complete Provider records and approved source-reference
transformations, so they remain future provenance rather than being made
executable with substitute chemistry.

## Architecture baseline

- executable campaign: 546 lines in one script;
- persistent tests: 74 lines in one parameterized test module;
- campaign documentation, source records, and generated JSON: 645 lines in
  ten files;
- future carbonate/complexation provenance: 254 lines in eight non-executable
  files;
- generated policy: the historical campaign generated the four detailed
  result owners and the aggregate summary;
- public production imports: `epcsaft` and `epcsaft_equilibrium`;
- direct dependencies: the Python standard library and the two exact installed
  wheels; and
- native modules or targets: none.

The 546-line campaign is the only file above 500 lines. It owns four distinct
source transformations plus wheel-`RECORD` verification; splitting those into
a validation-owned framework would add a second abstraction without another
campaign contract.
