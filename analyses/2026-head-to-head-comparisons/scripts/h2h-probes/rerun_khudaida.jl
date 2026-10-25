H2H_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
delete!(H2H_RESULT, "khudaida"); include(joinpath(H2H_ROOT, "analyses/2026-head-to-head-comparisons/scripts/h2h_clapeyron_tasks.jl"))
