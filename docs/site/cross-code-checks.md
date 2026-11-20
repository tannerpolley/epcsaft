# Cross-code check values

{download}`Download the cross-code check values <data/cross-code-check-values.csv>`.
The table contains retained pure-propane saturation and critical-equation
values, methane/ethane critical and bubble/dew values, electrolyte liquid compositions, solved-state derivatives, and fixed-state
Born values. These comparisons check numerical implementations at the stated
conditions; they are not measurements or physical validation.

Each row names the state, quantity and unit, the ePC-SAFT value, the other code
and its version, that code's value, an agreement bound and the conventions to
match. Values are copied as retained. Bounds identify whether differences are absolute
or relative and whether the reference value is the denominator. The propane
saturation pressure bound is the reported spread; the density, critical and
bubble/dew bounds are comparison criteria. The electrolyte bound is the
reported maximum difference in mole fraction. Derivative bounds are relative
to ePC-SAFT unless an absolute term is stated. For Born values, the bound
applies after matching constants, while the listed values retain the constants
each code actually used. These bounds are not uncertainty estimates. Source revisions marked as modified are not uniquely reproducible
from the commit alone.

Matching the stated parameters and model terms is necessary to reproduce these
values. The electrolyte rows use constant permittivity and the direct Born
term; the codes used slightly different electrostatic constants. Electrostatic
constants must also be matched for closer agreement in the Born term.

The critical point is a solution of the critical equations, not a stability
certificate or a measured critical point.
