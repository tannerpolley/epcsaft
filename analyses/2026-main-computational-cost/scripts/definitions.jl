# Reuse the cross-check's source-backed PC-SAFT/ePC-SAFT parameter conversions.
include("/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT/.worktrees/issue-533-full-speed-comparison-with-clapeyron-and/analyses/2026-competitor-cross-check/scripts/clapeyron.jl")

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


  # Clapeyron promotes number types through every submodel
function khudaida_model(spec, permittivity = "salt-free-suppressed")
    cs = spec["components"]; names = String[c["id"] for c in cs]; n = length(cs); T = Float64(spec["T"])
    Z = Int[Int(c["z"]) for c in cs]
    sigT(c) = c["id"] == "water" ? c["sigma"] + CP.Δσh20(T)*1e10 : c["sigma"]  # cross-volume size factor at T (Engine rule)
    epsab = Dict{Tuple{Tuple{String,String},Tuple{String,String}},Float64}(); kap = copy(epsab)
    for i in 1:n, j in i:n
        ci, cj = cs[i], cs[j]
        ci["kappa_ab"]*cj["kappa_ab"] == 0 && continue
        si, sj = sigT(ci), sigT(cj)
        kij = sqrt(ci["kappa_ab"]*cj["kappa_ab"])*(sqrt(si*sj)/((si+sj)/2))^3
        eij = (ci["epsilon_ab"]+cj["epsilon_ab"])/2
        for (a,b) in (("H","e"),("e","H"))
            epsab[((names[i],a),(names[j],b))] = eij; kap[((names[i],a),(names[j],b))] = kij
        end
    end
    sites = [c["kappa_ab"] > 0 ? 1 : 0 for c in cs]
    k = [Float64(spec["k"][i][j]) for i in 1:n, j in 1:n]
    pars = (Mw=Float64[c["Mw"] for c in cs], segment=Float64[c["m"] for c in cs], sigma=Float64[c["sigma"] for c in cs],
            epsilon=Float64[c["epsilon"] for c in cs], n_H=sites, n_e=sites, k=k, epsilon_assoc=epsab, bondvol=kap)
    neutral = pharmaPCSAFT(names; userlocations=pars, assoc_options=ASSOC)
    for i in 1:n, j in 1:n
        Z[i]*Z[j] > 0 && (neutral.params.epsilon.values[i,j] = 0.)  # same-sign ion pairs, as the Engine family
    end
    rsp = permittivity == "constant" ? ConstRSP(78.09) : SaltFreeSuppressedRSP(names, (epsilon=CP.SingleParam("epsilon", names, Float64[something(c["relative_permittivity"], 0.) for c in cs]),
        Mw=CP.SingleParam("Mw", names, Float64[c["Mw"] for c in cs]), charge=CP.SingleParam("charge", names, Float64.(Z))), Float64(spec["ion_fraction_suppression"]))
    solvents, ions = names[Z .== 0], names[Z .!= 0]
    ion = DHBorn(solvents, ions; RSPmodel=rsp, userlocations=(; sigma_born=Float64[something(c["born_diameter"], 0.) for c in cs]))
    return CP.ePCSAFT(names, Z, BasicIdeal(names), neutral, ion, ["10.1016/j.cherd.2014.05.017"])
end

function lle_residuals(model, P, T, f)
    Z = model.charge; μ = [CP.VT_chemical_potential(model, v, T, x)./(CP.Rgas(model)*T) for (v,x) in zip(f.volumes, f.compositions)]
    neutral = maximum(abs(μ[1][i]-μ[2][i]) for i in eachindex(Z) if Z[i] == 0)
    salt = abs(sum(μ[1][Z .!= 0]) - sum(μ[2][Z .!= 0]))
    pres = maximum(abs(pressure(model, v, T, x)/P - 1) for (v,x) in zip(f.volumes, f.compositions))
    charge = maximum(abs(sum(Z .* x)) for x in f.compositions)
    return Dict("neutral_mu_over_RT"=>neutral, "salt_mu_over_RT"=>salt, "pressure_relative"=>pres, "phase_charge"=>charge)
end


function with_parameter(model, which, θ)
    m = CP.promote_model(typeof(θ), model)
    if which == "k_ij"
        K = zeros(typeof(θ), 2, 2); K[1,2] = K[2,1] = θ; CP.set_k!(m, K)
    else  # ethanol self-association: the only nonzero association entries in C2
        p = which == "association_energy_over_k" ? m.params.epsilon_assoc : m.params.bondvol
        v = CP.raw_values(p); v[.!iszero.(CP.raw_values(which == "association_energy_over_k" ? model.params.epsilon_assoc : model.params.bondvol))] .= θ
    end
    return m
end
