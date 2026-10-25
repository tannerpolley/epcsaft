# Reuse the cross-check's source-backed PC-SAFT/ePC-SAFT parameter conversions.
include(joinpath(@__DIR__, "../../2026-competitor-cross-check/scripts/clapeyron.jl"))

const PRIMARY_MODE = "unique states (one pass)"
const REPEATED_MODE = "repeated identical input (rewards caching)"
const ACCURACY_LIMIT = 1e-10
const ABSENT_VAPOR_NOTE = "Requested vapor branch absent; volume(... phase=:vapor) returned NaN"

function request_for(case)
    Dict("case"=>case["model"], "diameter"=>"sigma", "electrostatic"=>"matched",
        "kcl"=>-0.25, "dipole_convention"=>"matched")
end

function prepare_model(document, case, state)
    spec = deepcopy(document["comparator_inputs"]["models"][case["model"]])
    if case["model"] == "water"
        spec["components"][1]["sigma"] = Float64(state["sigma_angstrom"])
    end
    request = request_for(case)
    model, _ = build_model(spec, request, document["comparator_inputs"]["constants"])
    return model
end

function prepare_runtime(model, state)
    return (model, Float64(state["P_Pa"]), Float64(state["T_K"]), Float64.(state["x"]), Symbol(state["branch"]))
end

function density_sample(mode, sample_index, state, elapsed_us, rho, status, note="")
    oracle = Float64(state["oracle_rho_mol_m3"])
    delta = isnothing(rho) ? nothing : abs(rho / oracle - 1)
    return Dict(
        "mode"=>mode,
        "sample_index"=>sample_index,
        "state_index"=>Int(state["index"]),
        "elapsed_us"=>elapsed_us,
        "rho_mol_m3"=>rho,
        "oracle_rho_mol_m3"=>oracle,
        "oracle_branch"=>state["oracle_branch"],
        "oracle_stable_root_count"=>Int(state["oracle_stable_root_count"]),
        "relative_density_difference"=>delta,
        "accuracy_within_1e10"=>(isnothing(delta) ? nothing : delta <= ACCURACY_LIMIT),
        "status"=>status,
        "note"=>note,
    )
end

function warm_volume(runtime)
    model, p, T, z, phase = runtime
    volume(model, p, T, z; phase=phase)
end

function timed_volume(runtime, mode, sample_index, state)
    model, p, T, z, phase = runtime
    begun = time_ns()
    local v
    try
        v = volume(model, p, T, z; phase=phase)
    catch error
        elapsed = Float64(time_ns() - begun) / 1000
        return density_sample(mode, sample_index, state, elapsed, nothing, "solve_failed",
            string(typeof(error), ": ", sprint(showerror, error)))
    end
    elapsed = Float64(time_ns() - begun) / 1000
    local vf
    try
        vf = Float64(v)
    catch error
        return density_sample(mode, sample_index, state, elapsed, nothing, "solve_failed",
            string(typeof(error), ": ", sprint(showerror, error)))
    end
    if !isfinite(vf)
        if isnan(vf) && state["branch"] == "vapor"
            return density_sample(mode, sample_index, state, elapsed, nothing, "unrepresentable", ABSENT_VAPOR_NOTE)
        end
        return density_sample(mode, sample_index, state, elapsed, nothing, "solve_failed", string("invalid volume ", vf))
    end
    rho = inv(vf)
    if !isfinite(rho) || rho <= 0
        return density_sample(mode, sample_index, state, elapsed, nothing, "solve_failed", string("invalid density ", rho))
    end
    return density_sample(mode, sample_index, state, elapsed, rho, "represented")
end

function warmups(runtime, errors)
    for _ in 1:5
        if isnothing(runtime)
            push!(errors, "model preparation failed")
            continue
        end
        try
            warm_volume(runtime)
        catch error
            push!(errors, string(typeof(error), ": ", sprint(showerror, error)))
        end
    end
end

function measure_case(document, case)
    states = case["states"]
    runtimes = Union{Nothing,Tuple}[]
    preparation_errors = String[]
    models = Dict{Any,Any}()
    for state in states
        try
            sigma = case["model"] == "water" ? Float64(state["sigma_angstrom"]) : nothing
            if !haskey(models, sigma)
                models[sigma] = prepare_model(document, case, state)
            end
            push!(runtimes, prepare_runtime(models[sigma], state))
            push!(preparation_errors, "")
        catch error
            push!(runtimes, nothing)
            push!(preparation_errors, string(typeof(error), ": ", sprint(showerror, error)))
        end
    end
    repeated = nothing
    repeated_error = ""
    try
        repeated_model = prepare_model(document, case, states[1])
        repeated = prepare_runtime(repeated_model, states[1])
    catch error
        repeated_error = string(typeof(error), ": ", sprint(showerror, error))
    end

    warmup_errors = Dict(PRIMARY_MODE=>String[], REPEATED_MODE=>String[])
    warmups(runtimes[1], warmup_errors[PRIMARY_MODE])
    warmups(repeated, warmup_errors[REPEATED_MODE])

    samples = Any[]
    for i in eachindex(states)
        state = states[i]
        runtime = runtimes[i]
        if isnothing(runtime)
            push!(samples, density_sample(PRIMARY_MODE, 0, state, nothing, nothing, "solve_failed", preparation_errors[i]))
        else
            push!(samples, timed_volume(runtime, PRIMARY_MODE, 0, state))
        end
    end
    for sample_index in 0:(Int(case["N"]) - 1)
        if isnothing(repeated)
            push!(samples, density_sample(REPEATED_MODE, sample_index, states[1], nothing, nothing, "solve_failed", repeated_error))
        else
            push!(samples, timed_volume(repeated, REPEATED_MODE, sample_index, states[1]))
        end
    end
    return Dict("id"=>case["id"], "samples"=>samples, "warmup_errors"=>warmup_errors)
end

length(ARGS) == 4 && ARGS[1] == "--worker" || error("launch through competitors.py")
replica = parse(Int, ARGS[2])
interpreter = ARGS[3]
effective_nice = parse(Int, ARGS[4])
identity = Dict(
    "type"=>"ready",
    "replica"=>replica,
    "pid"=>getpid(),
    "version"=>string(pkgversion(Clapeyron)),
    "effective_nice"=>effective_nice,
    "interpreter"=>interpreter,
)
println(JSON.json(identity))
flush(stdout)

for line in eachline(stdin)
    request = nothing
    try
        request = JSON.parse(line)
        result = measure_case(request["document"], request["case"])
        println(JSON.json(result))
    catch error
        case_id = isnothing(request) ? "" : get(request["case"], "id", "")
        println(JSON.json(Dict("id"=>case_id, "error_type"=>string(typeof(error)), "error"=>sprint(showerror, error))))
    end
    flush(stdout)
end
