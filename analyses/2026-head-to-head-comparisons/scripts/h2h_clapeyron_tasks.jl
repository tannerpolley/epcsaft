# Clapeyron 0.6.29 side of the 2026-10-08 head-to-head; included by h2h_clapeyron_server.jl.
# Reads the Engine-declared inputs (build/h2h/engine-h2h.json), writes build/h2h/clapeyron-h2h.json.
H2H_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
H2H_IN = JSON.parsefile(joinpath(H2H_ROOT, "build/h2h/engine-h2h.json"))
H2H_CASES = JSON.parsefile(joinpath(H2H_ROOT, "analyses/2026-competitor-cross-check/inputs/cases.json"))
H2H_OUT = joinpath(H2H_ROOT, "build/h2h/clapeyron-h2h.json")
isdefined(Main, :H2H_RESULT) || (H2H_RESULT = isfile(H2H_OUT) ? JSON.parsefile(H2H_OUT; dicttype=Dict{String,Any}) : Dict{String,Any}())  # resume finished tasks

function h2h_timed(f; warm = 3, samples = 10)
    t = time_ns(); value = f(); first = (time_ns()-t)/1e6
    for _ in 1:warm; f(); end
    ns = Float64[]
    for _ in 1:samples; t = time_ns(); f(); push!(ns, (time_ns()-t)/1e6); end
    sort!(ns)
    return value, Dict("first_ms"=>first, "median_ms"=>(ns[(samples+1)÷2]+ns[samples÷2+1])/2, "min_ms"=>ns[1], "max_ms"=>ns[end], "samples"=>samples)
end

# ---- 2. Khudaida electrolyte LLE, matched variants ----
Base.eltype(::SaltFreeSuppressedRSP) = Float64  # Clapeyron promotes number types through every submodel
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

function run_khudaida()
    spec = H2H_IN["clapeyron_khudaida"]; t0 = time_ns(); out = Dict{String,Any}()
    T, P = Float64(spec["T"]), Float64(spec["P"])
    for (variant, permittivity) in (("matched", "salt-free-suppressed"), ("constant-permittivity", "constant"))
        t1 = time_ns(); model = khudaida_model(spec, permittivity); build_s = (time_ns()-t1)/1e9; rows = []
        for (point, feed) in sort(collect(spec["feeds"]); by = p -> parse(Int, p[1]))
            z = Float64.(feed); row = Dict{String,Any}("point"=>parse(Int, point))
            method = MichelsenTPFlash(equilibrium=:lle, K_tol=1e-12, ss_iters=1000, nacc=5, second_order=true)
            guess = MichelsenTPFlash(equilibrium=:lle, K_tol=1e-12, ss_iters=1000, nacc=5, second_order=true,
                x0=[0.80,0.05,0.05,0.05,0.05], y0=[0.20,0.10,0.6998,0.0001,0.0001])
            for (label, m) in (("default-initialization", method), ("engine-guesses", guess))
                try
                    f, cost = h2h_timed(() -> CP.tp_flash2(model, P, T, z, m))
                    # Raw return and its diagnostics only; head_to_head.py classifies validity and phase count.
                    ix = sortperm([x[1] for x in f.compositions]; rev=true)
                    row[label] = Dict("x"=>f.compositions[ix], "rho"=>1 ./ f.volumes[ix],
                        "beta"=>f.fractions[ix[1]]/sum(f.fractions), "residuals"=>lle_residuals(model, P, T, f), "cost"=>cost)
                catch e
                    row[label] = Dict("error"=>sprint(showerror, e))
                end
            end
            push!(rows, row); println("khudaida ", variant, " point ", point); flush(stdout)
        end
        # DH chemical-potential check at the Engine's solved phase states: automatic vs central difference.
        eng = H2H_IN["khudaida"][variant][1]; checks = []
        for (x, rho) in zip(eng["x"], eng["rho"])
            x = Float64.(x); V = 1/Float64(rho)
            fdh(n) = sum(n)*CP.a_dh(model.ionmodel, V, T, n, CP.iondata(model, V, T, n))
            push!(checks, Dict("x"=>x, "automatic"=>ForwardDiff.gradient(fdh, x),
                "central_difference"=>[(fdh(x .+ 1e-6 .* (1:5 .== i)) - fdh(x .- 1e-6 .* (1:5 .== i)))/2e-6 for i in 1:5]))
        end
        out[variant] = Dict("rows"=>rows, "model_build_s"=>build_s, "dh_mu_check_point1"=>checks)
    end
    out["wall_s"] = (time_ns()-t0)/1e9
    H2H_RESULT["khudaida"] = out
end

# ---- 3. Implicit (IFT) derivatives through bubble_pressure and tp_flash2 ----
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

function run_sensitivities()
    out = Dict{String,Any}(); t0 = time_ns()
    for (ident, coords) in (("C2-0.5", ["k_ij","association_energy_over_k","association_volume"]), ("C4-298.15", ["k_ij"]))
        q = only(r for r in H2H_CASES["requests"] if r["id"] == ident); spec = H2H_CASES["models"][q["case"]]
        model, _ = build_model(spec, q, H2H_CASES["constants"]); T0 = Float64(q["T"]); x = Float64.(q["x"]); c0 = spec["components"][1]
        lle = MichelsenTPFlash(equilibrium=:lle, K_tol=1e-12, ss_iters=1000, nacc=5, second_order=true)
        obs(m, T) = if q["case"] == "C2"
            r = bubble_pressure(m, T, x, BUB); [r[1], r[4][1]]
        else
            f = CP.tp_flash2(m, Float64(q["P"]), T, x, lle)
            ix = sortperm([ForwardDiff.value(c[1]) for c in f.compositions]; rev=true)
            [f.compositions[ix[1]][1], f.compositions[ix[2]][1], f.fractions[ix[1]]/sum(f.fractions)]
        end
        origin = Dict("k_ij"=>Float64(spec["k"]), "association_energy_over_k"=>Float64(get(c0,"epsilon_ab",0.)), "association_volume"=>Float64(get(c0,"kappa_ab",0.)))
        rows = []; central, primal_cost = h2h_timed(() -> obs(model, T0))
        d, cost = h2h_timed(() -> ForwardDiff.derivative(t -> obs(model, t), T0))
        push!(rows, Dict("direction"=>"T", "derivative"=>d, "cost"=>cost))
        for c in coords
            d, cost = h2h_timed(() -> ForwardDiff.derivative(θ -> obs(with_parameter(model, c, θ), T0), origin[c]))
            push!(rows, Dict("direction"=>c, "derivative"=>d, "cost"=>cost))
        end
        g, gcost = h2h_timed(() -> ForwardDiff.jacobian(v -> vcat(obs(with_parameter(model, coords[1], v[1]), v[end]))[:], [origin[coords[1]], T0]))
        out[ident] = Dict("central"=>central, "primal_cost"=>primal_cost, "rows"=>rows, "two_direction_jacobian_cost"=>gcost)
        println("sensitivity ", ident); flush(stdout)
    end
    out["wall_s"] = (time_ns()-t0)/1e9
    H2H_RESULT["sensitivities"] = out
end

# ---- 4. Born term and permittivity at Engine-declared fixed states ----
function born_model(form)
    cs = H2H_IN["born_comparators"]["components"]; names = String[c["id"] for c in cs]; Z = Int[Int(c["z"]) for c in cs]
    k = [0. -0.3 -0.3; -0.3 0. 0.8; -0.3 0.8 0.]
    epsab = Dict(((names[1],"H"),(names[1],"e"))=>2425.7); kap = Dict(((names[1],"H"),(names[1],"e"))=>0.04509)
    pars = (Mw=[18.0153, 22.98977, 35.45], segment=Float64[c["m"] for c in cs], sigma=Float64[c["sigma"] for c in cs],
            epsilon=Float64[c["epsilon"] for c in cs], n_H=[1,0,0], n_e=[1,0,0], k=k, epsilon_assoc=epsab, bondvol=kap)
    neutral = pharmaPCSAFT(names; userlocations=pars, assoc_options=ASSOC)
    rsp = form == "mole-fraction-mixing" ? LinMixRSP(names[1:1], names[2:3]; userlocations=(; dielectric=Float64[c["relative_permittivity"] for c in cs])) :
        Schreckenberg(names[1:1], names[2:3]; userlocations=(; d_T=Float64[something(c["dielectric_temperature"], 0.) for c in cs], d_V=Float64[something(c["dielectric_volume_m3_per_mol"], 0.) for c in cs]))
    ion = DHBorn(names[1:1], names[2:3]; RSPmodel=rsp, userlocations=(; sigma_born=Float64[something(c["born_diameter"], 0.) for c in cs]))
    return CP.ePCSAFT(names, Z, BasicIdeal(names), neutral, ion, ["10.1016/j.cherd.2014.05.017"])
end

function run_born()
    rows = []; t0 = time_ns()
    for form in ("mole-fraction-mixing", "schreckenberg")
        model = born_model(form)
        for st in H2H_IN["born_states"]
            V, T, x = 1/Float64(st["rho"]), Float64(st["T"]), Float64.(st["x"])
            d = CP.iondata(model, V, T, x)
            push!(rows, Dict("form"=>form, "m"=>st["m"], "permittivity"=>CP.dielectric_constant(model, V, T, x),
                "a_born"=>CP.a_born(model.ionmodel, V, T, x, d), "a_dh_at_sigma"=>CP.a_dh(model.ionmodel, V, T, x, d)))
        end
    end
    H2H_RESULT["born"] = rows; H2H_RESULT["born_wall_s"] = (time_ns()-t0)/1e9
    H2H_RESULT["constants"] = Dict("e"=>CP.e_c, "epsilon0"=>CP.ϵ_0, "kB"=>CP.k_B, "NA"=>CP.N_A)
end

for (name, f) in (("born", run_born), ("sensitivities", run_sensitivities), ("khudaida", run_khudaida))
    haskey(H2H_RESULT, name) && continue
    t = time_ns(); delete!(H2H_RESULT, name*"_error")
    try
        f()
    catch e
        H2H_RESULT[name*"_error"] = sprint(showerror, e, catch_backtrace()); println(stderr, name, ": ", sprint(showerror, e))
    end
    H2H_RESULT[name*"_task_wall_s"] = (time_ns()-t)/1e9
    write(H2H_OUT, JSON.json(H2H_RESULT; allownan=true))
end
H2H_RESULT["julia"] = Dict("version"=>string(VERSION), "clapeyron"=>string(pkgversion(Clapeyron)), "compile_enabled"=>Base.JLOptions().compile_enabled, "opt_level"=>Base.JLOptions().opt_level)
write(H2H_OUT, JSON.json(H2H_RESULT; allownan=true))
println("tasks written ", H2H_OUT); flush(stdout)
