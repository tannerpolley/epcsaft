using Clapeyron
using JSON

T = parse(Float64, get(ENV, "EPCSAFT_COMPARISON_T", ARGS[1]))
model = PCSAFT("propane")
p, v_l, v_v = saturation_pressure(model, T, ChemPotVSaturation())
result_json = JSON.json(Dict(
    "pressure_Pa" => p,
    "rho_liquid_mol_m3" => 1 / v_l,
    "rho_vapor_mol_m3" => 1 / v_v,
    "phase_pressure_rel_error" => nothing,
))

println(result_json)
if haskey(ENV, "CLAPEYRON_OUTPUT")
    write(ENV["CLAPEYRON_OUTPUT"], result_json * "\n")
end
