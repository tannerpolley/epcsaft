# Reuse the existing one-process FIFO server; relocate declarations to this leaf.
const LEAF681 = abspath(joinpath(@__DIR__, ".."))
const REPO681 = abspath(joinpath(LEAF681, "../.."))
const SOURCE681 = joinpath(REPO681, "analyses/2026-audit-timing-campaign/scripts/clapeyron_server.jl")
text681 = read(SOURCE681, String)
text681 = replace(text681,
    "include(joinpath(@__DIR__, \"../../2026-competitor-cross-check/scripts/clapeyron.jl\"))" =>
    "include(\"$(REPO681)/analyses/2026-competitor-cross-check/scripts/clapeyron.jl\")",
    "JSON.parsefile(joinpath(@__DIR__, \"../inputs/inputs.json\"))" =>
    "JSON.parsefile(\"$(LEAF681)/inputs/inputs.json\")")
text681 = replace(text681, r"const DEF533 = .*" => "const DEF533 = \"$(LEAF681)/scripts/definitions.jl\"")
include_string(Main, text681, SOURCE681)
