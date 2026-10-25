ENV["EPCSAFT_COMPARISON_T"] = "300.0"
ENV["CLAPEYRON_OUTPUT"] = joinpath(@__DIR__, "..", "results", "main-8298c602", "clapeyron.json")
include(joinpath(@__DIR__, "collect_clapeyron.jl"))
delete!(ENV, "EPCSAFT_COMPARISON_T")
delete!(ENV, "CLAPEYRON_OUTPUT")
