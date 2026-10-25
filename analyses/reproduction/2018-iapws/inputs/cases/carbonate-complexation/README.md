# Chemical-equilibrium preintegration sources

This packet prepares two independent reference families for future
fixed-temperature, fixed-pressure homogeneous chemical validation:

- an ideal-dilute carbonate ladder using the 25 °C concentration constants in
  IAPWS TGD9-18; and
- calcium carbonate/bicarbonate complex records from the `wateq4f.dat`
  database distributed with PHREEQC 3.9.0.

The IAPWS result table contains three closed, electroneutral analytical
benchmarks at total inorganic carbon `0.001 mol/L`. Sodium is an inert charge
carrier. The three sodium totals select carbonic-acid-, bicarbonate-, and
carbonate-dominant regimes without prescribing pH.

The PHREEQC records preserve reaction orientation, base-10 constants,
temperature expressions, and the exact database identity. Published PHREEQC
Example 1 values establish that both calcium complexes are material in a broad
seawater calculation. They are cross-model context, not expected ePC-SAFT
amounts.

The retained PHREEQC input also runs four closed aqueous solutions through the
official `3.9.0-17591` source release and its bundled `wateq4f.dat`. The result
checks the three carbonate-dominance regimes and confirms nonzero
`CaHCO3+`/`CaCO3` complex formation while closing carbon, sodium, calcium, and
charge balances. PHREEQC reports these quantities as molality. The release
archive identifies itself as 3.9.0, but its executable banner still says
`PHREEQC_3.8.9`; the source contract preserves both identities.

The packet deliberately does not:

- transform either source into the Provider chemical-reference basis;
- invent Provider components or neutral-reference records;
- use PHREEQC activity coefficients as ePC-SAFT activity coefficients;
- import a private Equilibrium binding;
- include calcite precipitation or another phase-appearance problem; or
- claim that an installed chemical-equilibrium calculation has run.

This packet does not execute the public chemical-equilibrium route because no
approved Provider-basis records exist for these
carbonate and calcium-complexation reactions. A future campaign may consume
the packet once those records and transforms are source-complete; the PHREEQC
result remains an independent cross-model check rather than an acceptance
target.
