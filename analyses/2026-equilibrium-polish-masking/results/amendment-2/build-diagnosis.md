# Manufactured Newton admission check compilation

The first candidate build stopped when `specification_checks.cpp` included
`gibbs_minimizer.hpp`: that header requires `IpTNLP.hpp`, while the specification
check target declares no Ipopt include dependency. The shared native solver
library compiled before this failure, but no installed candidate qualification
was produced.

The correction moves only the shared admission function declaration into the
existing `problem.hpp`, which both the solver and the manufactured check already
consume. The function implementation remains in `gibbs_minimizer.cpp`.
No equation, numerical method, tolerance, physical input or build dependency changes.
`candidate-build-diagnostic.log` retains the exact compiler stop and verified
reference-wheel restoration.
