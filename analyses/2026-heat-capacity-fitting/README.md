# Total isobaric heat-capacity fitting

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The results used
a Python solver before the native calculation rule; the ideal-Cp projection and
coefficient/comparison writers were removed, so the current scripts cannot reproduce
them.

Pre-rule record: the retained results were produced by a Python solver before the native calculation rule (#234, 2026-10-01); they are history, not reproducible by the current scripts. The ideal-Cp polynomial projection and its coefficient and comparison writers were removed from `scripts/reference.py`.

The analysis asks what the exact first parameter derivatives establish and how
Held, Cameretti and Sadowski's 2008 2B water model can match density, total molar
isobaric heat capacity and saturation pressure jointly. Parameter adoption is
outside this numerical demonstration. The [notebook](notebook.qmd) reads retained
numbers and does not run the model.

The scripts read the packet-v3 parameter document from
`data/packets/held-cameretti-sadowski-2008/3` (#341 removed the byte-identical leaf copy);
only water is selected. `inputs/source.json` records its source path, SHA-256,
packet fingerprint, IAPWS release URL and Zotero key. The constant part of
sigma(T) is fitted; its exponential terms remain fixed. The native checks carry
the same Table 1 water subset. The source reported electrolyte conditions do not
establish validation of this pure-water demonstration over 280–440 K.

`reference.py` verifies chemicals 1.5.2 against IAPWS R6-95(2018), Table 3 and
Table 6 for cp, all Table 7 single-phase rows, Table 8 at 275 and 450 K, and the
ten retained NIST SRD 69 densities. Table 7's low-pressure pressure-rounding
caveat is retained; this run passes even its printed half-digit criterion. The
retained ideal-gas centered polynomial approximates the release's Eq. (5) over 270–450 K.
The only fitted quantities are the five residual-EOS coordinates. Density and
cp targets are compressed liquid at 1 MPa, with saturation pressure at each
280–440 K temperature in 20 K steps. Residuals are equally weighted LogRatio.
The companion fit omits cp rows and uses the same chart and default controls.

The original non-editable-wheel calculation recipe is withdrawn; the retained
methods and output locations below identify the historical calculation.

The native `regression_checks cp-derivatives` and `cp-recovery` output goes to
`build/issue205-derivatives.log` and `build/issue205-recovery.log`; `scripts/retain.py`
converts the comparisons to CSV. `scripts/parameter-sensitivities.py` retains the diagnosis of the
original signal-floor conflict. For only water dispersion energy at 298.15 K and
the binary cross-site association energy, the owner approved a stricter relative
ladder with resolved refinement estimates. All other signal floors, states,
steps, value identities and density-lift checks remain as designed in #205.
`plot.py` renders retained results without running thermodynamics.

References: Held, Cameretti and Sadowski (2008), *Fluid Phase Equilibria* 270,
87–96, Table 1 and Eq. (7), [DOI](https://doi.org/10.1016/j.fluid.2008.06.010);
IAPWS (2018), *R6-95(2018)*, [official release](https://iapws.org/technical-guidance/release/IAPWS-95),
Zotero item `Y3FU9CXJ`. The nine targets are evaluated reference values, not
independent measurements. Conditional covariance and local rank do not establish
unique parameters or predictive accuracy outside these inputs.
