# Incomplete HCl parameter estimation

The five-coordinate native Ceres fits all converged. Their returned spectra have
full numerical rank, but every fit ended on a numerical bound. The lowest-cost
pair puts the independent Born diameter exactly at its lower trial bound.
The covariance owner consequently withholds its inverse-information matrix.
This is a failure of the declared interior-identification requirement, not
measurement evidence for a proton Born diameter.

The first declared reduction ties Born and Debye–Hückel diameters. Its middle
start converged at an interior endpoint with available inverse information.
Both the free-coordinate correlation and the normalized local excursions exceed
their predeclared limits, as retained in `identification.csv` and
`local-resolution.csv`. The other allocated starts exhausted their iteration
limits; they were not wall-time stops. Their objective histories, endpoints,
exact columns and any rejected trial causes are retained in `fit-values.csv`
and `trial-failures.csv`.

These facts separate two causes. The completed five-coordinate stage fails
identification at its selected endpoint. The tied-diameter stage is incomplete
because its allocated starts did not all converge. It also contains a converged
endpoint that fails practical identification. The retained histories do not
establish that a longer run would converge, that the final reductions would
pass, or that the model cannot represent the data.

No extra iterations, restart, fit method, weight, bound or restriction is
permitted by the frozen declaration. The design explicitly returns an incomplete
stage to the lead rather than advancing its reduction sequence. A numerical
root-cause correction or a larger work allocation requires disposition and a
new frozen declaration; none was attempted here. This is an unfinished study,
not a supported negative conclusion about physical model adequacy.

The earlier subdirectories retain driver declaration errors. They all have zero
iterations and evaluations: the first used component indices in native pair
slots, the next reversed the canonical reciprocal order, and the third omitted
the liquid branch. The current producer uses flattened reciprocal slots in
canonical order, explicitly declares the liquid branch and validates the full
request before optimization. The declaration errors do not consume numerical
start budgets or support any scientific inference.

Raw correlation, excursion and distance fields in `identification.csv` are
meaningful only when their native prerequisites are available. Zero placeholders
in rows whose covariance was withheld are unavailable quantities, not zero
correlation or perfect resolution. The notebook presents them as unavailable.
No identified endpoint was selected for binary accuracy or mixture
qualification. There is no qualified parameter packet and no claim of mixture
accuracy. Reserved source rows remain unused by the fitting objective.
