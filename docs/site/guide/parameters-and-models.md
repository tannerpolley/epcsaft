# Parameters and models

A parameter set declares the physical model as well as its numbers. Select
species, parameter source and version, model families, association sites and
unlike pairs before constructing a model. The Engine does not select physics
from an application name. The
[water notebook](../examples/parameters-and-water-state.ipynb) demonstrates a
published source selection; the [parameter algorithm](../theory/algorithms.md#alg-parameters)
and [EOS conventions](../theory/eos.md) define the units and model terms.

## Select source records

`Parameters.from_json` reads a parameter document with an optional component
selection. `Parameters.from_packets` combines named component selections from
versioned packets below a supplied directory. `parameter_document` exposes the
combined document before native construction. These calls are documented in
[parameter inputs](../api/parameters.rst).

Search the [1.3.2 component index](https://github.com/tannerpolley/epcsaft/blob/v1.3.2/data/docs/parameter-index.md)
by component name, then use the packet ID and version as `packet/version`.
Copy the exact `component_id` values from the index or the packet's
`parameters/parameters.json` (`components[].component_id`). A missing name
raises `ParameterError` and identifies the packet. Keep the selected packet versions,
component order, source locators and caller overrides with the calculation.
A combined document retains materialization source information. A successful
load establishes mathematical admission, not the validity of extrapolation.

| Quantity | Input unit |
|---|---|
| Segment count, unlike `k_ij` and `l_ij` | dimensionless |
| Segment and ionic diameters | Å |
| Dispersion and association energy divided by Boltzmann's constant | K |
| Association volume | dimensionless |
| Molar mass | kg/mol |
| State pressure and molar density | Pa and mol/m³ |

Source conversions must be explicit. The parameter reader and native builder
check the declared units; do not silently treat an Å coefficient as metres.

## Combine model choices deliberately

Packet combination refuses conflicting base, electrolyte, permittivity or
polar choices. A component without association sites can join a general-site
association model, but arbitrary model conflicts are not resolved by taking
the first source. Selecting two copies of the same component is also refused.

Sites and their interaction edges are explicit. A missing cross-association
edge is not permission to invent one. Declare the intended edge or the
source's deliberate omission. Unlike-component `k_ij` and `l_ij` overrides
refer to two selected, distinct components. A constant override replaces that
pair's corresponding temperature law as well as its constant coefficient.
Keep a caller override labelled as an assumption or a fitted value.

See the [ammonia parameter combination](../examples/index.md#ammonia-parameter-combination)
for an explicit combined declaration.

## Construct and reuse the model

Use `Pure(parameters)` for exactly one component and `Mixture(parameters)` for
two or more. `Pure` fixes the composition at one and offers a saturation
problem declaration. Both copy their complete physical inputs into the native
model. Editing a Python input afterwards does not update a prepared model.
Construct a new model for a different parameter set.

Pass thermal records in component order when the requested outputs consume
them. Declare active parameters when requesting parameter actions; selecting
an active coordinate changes what is differentiated, not the physical source
record. See [states](states-and-properties.md),
[parameter estimation](parameter-estimation.md), and the
[model reference](../api/eos.rst).

## Choose electrolyte and polar terms by their equations

The EOS includes the declared hard-chain, dispersion and association terms,
and optional electrostatic and polar contributions. A dipole-dipole residual
term and a bulk-permittivity choice are distinct declarations. The
[EOS theory](../theory/eos.md) defines the admitted families and exclusions.
See [aqueous NaCl inputs](../examples/index.md#aqueous-nacl-inputs)
for the Held source selection.

Charged ordinary states must be electroneutral. Mathematical admission does
not establish parameter accuracy at a temperature, concentration or solvent
composition. Read the source's fit domain and retain unknown applicability
as unknown; the [limits](../limits.md) explain this distinction.

## Installed and downloaded packets

Conda installs `feos-gross-sadowski-2001/1` for nonassociating fluids and
`feos-gross-sadowski-2002/1` for associating fluids. The packet root is
`Path(epcsaft.__file__).resolve().parent / "data" / "packets"`.
Source builds with the accessor expose this same path as
`epcsaft.installed_packet_root()`. The directory can be absent when the
installation contains no bundled packets.

```{doctest}
>>> from pathlib import Path
>>> import epcsaft
>>> assert epcsaft.installed_packet_root() == Path(epcsaft.__file__).resolve().parent / "data" / "packets"
>>> packets = Path("data/packets")  # extracted source or checkout
>>> epcsaft.Parameters.from_packets(packets, [("feos-gross-sadowski-2001/1", ["propan"])])
Traceback (most recent call last):
    ...
epcsaft._core.ParameterError: feos-gross-sadowski-2001/1 has no component propan; see components[].component_id in data/packets/feos-gross-sadowski-2001/1/parameters/parameters.json
```

For other packets, download the [matching source release](https://github.com/tannerpolley/epcsaft/releases)
and pass its `data/packets` directory to `Parameters.from_packets`. The
[NaCl quickstart](../quickstart.md#aqueous-nacl-activity) extracts one packet
from the 1.3.2 archive. Each packet's `packet.toml` and any
`parameters/LICENSE-NOTICE.md` state its use terms.

The source archive omits the Esper 2023 pure-component (`feos-esper-2023/1`) and
Rehner 2023 binary (`feos-rehner-2023/1`) packets because their use basis is
CC BY-NC 4.0, which permits noncommercial use only. Rebuild them from their
original FeOs distribution,
[feos-org/feos `parameters/pcsaft` at commit `af84e3f9`](https://github.com/feos-org/feos/tree/af84e3f9da2fac2d102d449cd07cbbb48237617b/parameters/pcsaft),
by running this from the extracted source root:

```bash
bash data/tools/restore_feos.sh
```

It downloads `gross2001.json`, `gross2002.json`, `esper2023.json` and
`rehner2023_binary.json` at that commit, refuses any file whose SHA-256 differs
from the pinned value, writes the four `feos-*/1` packets under `data/packets`
and regenerates `data/docs/parameter-index.md`. The Esper and Rehner packet
tests skip until this step has run. Cite the parameter sources:

- Esper, T.; Bauer, G.; Rehner, P.; Gross, J. PCP-SAFT Parameters of Pure
  Substances Using Large Experimental Databases. *Ind. Eng. Chem. Res.* **2023**,
  DOI [10.1021/acs.iecr.3c02255](https://doi.org/10.1021/acs.iecr.3c02255).
- Rehner, P.; Bardow, A.; Gross, J. Modeling Mixtures with PCP-SAFT: Insights
  from Large-Scale Parametrization and Group-Contribution Method for Binary
  Interaction Parameters. *Int. J. Thermophys.* **2023**, DOI
  [10.1007/s10765-023-03290-3](https://doi.org/10.1007/s10765-023-03290-3).

Selection and mathematical
admission do not establish a fit range or validate an extrapolation.

The polar model may select a quadrupole–quadrupole contribution. Only one
polar choice is active, and dipole–quadrupole mixing is not supported. Gross's
mixing rule is unsuitable for real mixtures with more than one quadrupolar
component. The [EOS theory](../theory/eos.md) describes the contribution.

## Caloric data availability

Both installed FeOs packets supply residual EOS parameters and omit ideal-gas
thermal records. The downloadable Held, Cameretti and Sadowski NaCl packet
also omits them. Total enthalpy requires an ideal-gas thermal record for each
component. The [propane quickstart](../quickstart.md#pure-propane-saturation)
shows the returned `MissingRecord` and `None` for total enthalpy.

The downloadable `nasa-glenn-2002-hda-species/1` packet supplies NASA9
ideal-gas records for hydrogen, methane, benzene and toluene. Load those
records separately from the EOS parameters. The
[methane caloric example](../examples/methane-caloric-properties.ipynb)
shows how to attach the methane record and read total and residual properties.

## Composition domain report

For packet-loaded parameters, call `mixture.outside_domains(x)` with one phase's
mole fractions in component order (including a phase returned by equilibrium).
Each returned dictionary names `packet`, `domain_id`, `quantity`, `limit`,
`bound`, `value`, `unit`, and `solute` for molality. Quantity and limit are native
`DomainQuantity` and `DomainLimit` enum values. Comparisons are strict.

Only declared bounds are reported. Temperature and pressure are outside this
report. A missing solvent or absent bounded solute skips its comparison. Each
packet contributes only bounds for its own selected solutes; ionic strength
still includes all ions in the phase. A bare parameter JSON uses its own domains;
fitted and directly constructed parameters have no packet bounds. Applications
choose how to respond to a report; calculations are unaffected.
