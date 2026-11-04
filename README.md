# ePC-SAFT

[![release](https://img.shields.io/github/v/release/tannerpolley/epcsaft)](https://github.com/tannerpolley/epcsaft/releases)
[![documentation](https://img.shields.io/badge/docs-github--pages-blue)](https://tannerpolley.github.io/epcsaft/)
[![licence](https://img.shields.io/badge/licence-GPLv3-blue)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23634228-blue)](https://doi.org/10.5281/zenodo.23634228)

A C++ implementation of the PC-SAFT and electrolyte PC-SAFT (ePC-SAFT)
equations of state with a Python interface. It computes fluid properties,
phase and reaction equilibrium, exact derivatives of solved states, and
parameter estimates. Every density root, equilibrium iteration and fit runs in
the native Engine; Python declares the calculation and reads the results.

```python
from epcsaft import Mixture, Parameters

# Water and ion parameters of Held, Cameretti and Sadowski (2008),
# read from data/packets in the source archive
parameters = Parameters.from_packets(
    "data/packets",
    [("held-cameretti-sadowski-2008/4", ["water", "sodium-cation", "chloride-anion"])],
)
brine = Mixture(parameters)
reference = brine.reference(
    298.15, 1.0e5,                          # T / K, P / Pa
    solvent_fractions=(1.0, 0.0, 0.0),      # pure water solvent
    ions=("sodium-cation", "chloride-anion"),
    stoichiometry=(1, 1),
)
for molality in (0.1, 1.0, 2.0):            # mol NaCl per kg water
    activity = reference.activity(molality)
    print(f"{molality:4.1f} mol/kg  gamma± = {activity.mean_ionic_molality_scale:.4f}"
          f"  phi = {activity.osmotic_coefficient:.4f}")
```

```text
 0.1 mol/kg  gamma± = 0.7631  phi = 0.9251
 1.0 mol/kg  gamma± = 0.6553  phi = 0.9528
 2.0 mol/kg  gamma± = 0.7046  phi = 1.0335
```

These are the model's mean ionic activity and osmotic coefficients of aqueous
NaCl at 298.15 K and 1 bar; the example makes no comparison with measurements.

## Models

- **PC-SAFT** hard-chain and dispersion terms
  ([Gross and Sadowski, 2001](https://doi.org/10.1021/ie0003887)).
- **Association** after Wertheim
  ([Chapman et al., 1990](https://doi.org/10.1021/ie00104a021);
  [Gross and Sadowski, 2002](https://doi.org/10.1021/ie010954d)).
- **Polar terms**, optional: dipole–dipole
  ([Gross and Vrabec, 2006](https://doi.org/10.1002/aic.10683)) and
  quadrupole–quadrupole ([Gross, 2005](https://doi.org/10.1002/aic.10502)).
- **Electrolytes (ePC-SAFT)**: Debye–Hückel term
  ([Cameretti, Sadowski and Mollerup, 2005](https://doi.org/10.1021/ie0488142);
  [Held, Cameretti and Sadowski, 2008](https://doi.org/10.1016/j.fluid.2008.06.010))
  with cation–anion dispersion
  ([Held et al., 2014](https://doi.org/10.1016/j.cherd.2014.05.017)).
- **Permittivity and Born term**, optional: composition-dependent permittivity
  ([Bülow et al., 2019](https://doi.org/10.1016/j.fluid.2019.03.010)) and Born
  solvation ([Bülow, Ascani and Held, 2021](https://doi.org/10.1016/j.fluid.2021.112989);
  [Figiel, Yu and Held, 2025](https://doi.org/10.1021/acs.iecr.5c00475)).

The [equation registry](https://tannerpolley.github.io/epcsaft/theory/equations.html)
gives every implemented equation with its source equation number.

## Capabilities

- **Properties:** density on a declared branch, pressure, fugacity and activity
  coefficients, molality-scale mean ionic activity and osmotic coefficients,
  residual properties, and total calorics when ideal-gas records are supplied.
- **Phase and reaction equilibrium:** flashes, bubble and dew points,
  saturation, liquid–liquid splits, chemical reactions and pure solids;
  phase-envelope tracing; tangent-plane stability tests and finite phase
  discovery (not a global stability proof).
- **Derivatives:** exact first and second derivatives of properties and solved
  equilibrium states with respect to conditions and parameters.
- **Parameter fitting:** least squares, errors-in-variables and guarded
  minimax, with local covariance.
- **Process models:** Pyomo/IDAES external functions with exact derivatives.

The [limits page](https://tannerpolley.github.io/epcsaft/limits.html) states
what each capability does not cover. The examples are numerical verification
unless they say they compare with measurements.

## Installation

Install release 1.3.0 for Linux and Windows (x86-64) from the
`tannerwpolley` conda channel:

```bash
conda create -n epcsaft -c conda-forge -c tannerwpolley python=3.13 epcsaft=1.3.0
```

The conda package includes the FeOs Gross–Sadowski 2001 and 2002 packets.
Download the matching [release source archive](https://github.com/tannerpolley/epcsaft/releases)
for other packets and example inputs; load packets from its `data/packets`
directory with `Parameters.from_packets`. Each packet states its use terms. To build from source, see
[installation](https://tannerpolley.github.io/epcsaft/installation.html).

## Documentation

The [documentation](https://tannerpolley.github.io/epcsaft/) has a quickstart,
a user guide, the theory and the Python reference. The
[examples](https://tannerpolley.github.io/epcsaft/examples/index.html) are
executed notebooks: water, methane, methane/ethane, aqueous NaCl, halite
solubility, ammonia synthesis, water/isobutanol phase discovery, propane
parameter estimation and an IDAES settler.

## Citing

Please cite the software with the metadata in [CITATION.cff](CITATION.cff)
(DOI [10.5281/zenodo.23634228](https://doi.org/10.5281/zenodo.23634228)), and
cite the model and parameter sources you use.

## Licence

GPLv3 with an additional permission for linking with EPL-licensed Ipopt
([LICENSE](LICENSE)); retained third-party data have their own terms
([NOTICE](NOTICE)).

## Contributing

See the [contributing guide](https://tannerpolley.github.io/epcsaft/contributing.html).
