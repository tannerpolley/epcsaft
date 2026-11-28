# ePC-SAFT

[![release](https://img.shields.io/github/v/release/tannerpolley/epcsaft)](https://github.com/tannerpolley/epcsaft/releases)
[![documentation](https://img.shields.io/badge/docs-github--pages-blue)](https://tannerpolley.github.io/epcsaft/)
[![licence](https://img.shields.io/badge/licence-GPLv3-blue)](LICENSE)

Calculate fluid properties, phase and reaction equilibrium, exact derivatives
and fitted parameters with PC-SAFT and electrolyte PC-SAFT. Python declares
inputs and reads results from the C++ Engine.

## Install

```bash
conda create -n epcsaft -c conda-forge -c tannerwpolley python=3.13 epcsaft=0.1.2
conda activate epcsaft
```

See the [installation guide](https://tannerpolley.github.io/epcsaft/installation.html)
for supported platforms and source installation.

## Aqueous NaCl activity

<!-- nacl-example-start -->
From an empty directory, download and extract the Held, Cameretti and
Sadowski (2008) parameter packet from the **0.1.2** source release. The conda
package includes two neutral-fluid packets.

```bash
curl -fLO https://github.com/tannerpolley/epcsaft/releases/download/v0.1.2/epcsaft-0.1.2.tar.gz
tar -xzf epcsaft-0.1.2.tar.gz --strip-components=1 epcsaft-0.1.2/data/packets/held-cameretti-sadowski-2008/4
```

In Python, calculate the molality-scale mean ionic activity coefficient
`gamma±` and the osmotic coefficient `phi` at 298.15 K and 1 bar.
Molality is in mol NaCl per kg water.

```doctest
>>> from epcsaft import Mixture, Parameters
>>> ions = ("sodium-cation", "chloride-anion")
>>> parameters = Parameters.from_packets(
...     "data/packets", [("held-cameretti-sadowski-2008/4", ["water", *ions])])
>>> brine = Mixture(parameters)
>>> reference = brine.reference(
...     298.15, 100_000.0, solvent_fractions=(1.0, 0.0, 0.0),
...     ions=ions, stoichiometry=(1, 1))
>>> for molality in (0.1, 1.0, 2.0):
...     activity = reference.activity(molality)
...     print(f"{molality:.1f} mol/kg  gamma± = {activity.mean_ionic_molality_scale:.4f}"
...           f"  phi = {activity.osmotic_coefficient:.4f}")
0.1 mol/kg  gamma± = 0.7631  phi = 0.9251
1.0 mol/kg  gamma± = 0.6553  phi = 0.9528
2.0 mol/kg  gamma± = 0.7046  phi = 1.0335

```

These coefficients are model predictions.
<!-- nacl-example-end -->

## Documentation

- [First calculations](https://tannerpolley.github.io/epcsaft/quickstart.html) covers pure saturation, mixture flash and electrolyte activity.
- [Parameter guide](https://tannerpolley.github.io/epcsaft/guide/parameters-and-models.html) explains how to choose components and source records.
- [Examples](https://tannerpolley.github.io/epcsaft/examples/index.html) include executed notebooks for state properties, equilibrium and parameter estimation.
- The [equation registry](https://tannerpolley.github.io/epcsaft/theory/equations.html) records the source equations. The [limits page](https://tannerpolley.github.io/epcsaft/limits.html) describes the scope of supported calculations.

Cite the software using [CITATION.cff](CITATION.cff) and cite the model and
parameter sources you use. See [LICENSE](LICENSE) and [NOTICE](NOTICE)
for software and third-party data terms.
