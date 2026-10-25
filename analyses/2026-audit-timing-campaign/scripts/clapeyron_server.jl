# One persistent Clapeyron 0.6.29 process for the whole #578 campaign.
# Reads request-file paths from the FIFO named by ARGS[1]; writes each reply JSON atomically to request["out"].
# Conversions: the #222 cross-check (clapeyron.jl), #533 density timing (definitions.jl) and the Opus P01 model.
include(joinpath(@__DIR__, "../../2026-competitor-cross-check/scripts/clapeyron.jl"))
const DEF533 = "/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT/build/issue-533/scripts/definitions.jl"
include_string(Main, replace(read(DEF533, String), r"^include\(.*\)$"m => ""), DEF533)
const INPUTS = JSON.parsefile(joinpath(@__DIR__, "../inputs/inputs.json"))
const ROOT = abspath(joinpath(@__DIR__, "../../.."))
const CROSS = JSON.parsefile(joinpath(ROOT, "analyses/2026-competitor-cross-check/inputs/cases.json"))
const DENSITY = JSON.parsefile(joinpath(ROOT, "analyses/2026-density-root-benchmark/inputs/cases.json"))
const T1 = Dict(r["id"] => r for r in INPUTS["t1"])
const RG = CP.R̄

# Opus probes/p01/clapeyron_models.jl:model, with the frozen component constants.
function audit_model(key; kij=0.0)
    comps = [INPUTS["components"][n] for n in INPUTS["models"][key]]
    n = length(comps); names = String[c["cid"] for c in comps]
    k = zeros(n, n); n == 2 && (k[1, 2] = k[2, 1] = kij)
    nH = [haskey(c, "assoc") ? 1 : 0 for c in comps]
    ea = Dict{Tuple{Tuple{String,String},Tuple{String,String}},Float64}(); bv = copy(ea)
    for c in comps
        haskey(c, "assoc") || continue
        ea[((c["cid"], "e"), (c["cid"], "H"))] = c["assoc"][1]; bv[((c["cid"], "e"), (c["cid"], "H"))] = c["assoc"][2]
    end
    base = (; Mw=Float64[c["M"] * 1000 for c in comps], segment=Float64[c["m"] for c in comps],
            sigma=Float64[c["sigma"] for c in comps], epsilon=Float64[c["eps"] for c in comps], k, n_H=nH, n_e=nH)
    base = isempty(ea) ? merge(base, (; epsilon_assoc=nothing, bondvol=nothing)) : merge(base, (; epsilon_assoc=ea, bondvol=bv))
    any(haskey(c, "mu") for c in comps) && return PCPSAFT(names; userlocations=merge(base, (; dipole=[get(c, "mu", 0.0) for c in comps])))
    return PCSAFT(names; userlocations=base)
end
t1_model(key) = key in ("gross-methane-ethane", "methane-ethane") ? audit_model("methane-ethane") : audit_model(key)

const PREPARED = Dict{String,Any}()

function make(req)
    op = req["op"]
    if op == "t1-bundle"
        r = T1[req["state"]]; m = t1_model(r["model"]); T = Float64(r["T"]); V = 1 / Float64(r["rho"]); x = Float64.(r["x"])
        return () -> begin
            p = pressure(m, V, T, x); Z = p * V / (RG * T)
            mu = CP.VT_chemical_potential_res(m, V, T, x) ./ (RG * T)
            cv = -T * ForwardDiff.derivative(t -> ForwardDiff.derivative(s -> CP.eos_res(m, V, s, x), t), T)
            Dict("P" => p, "h_res" => CP.VT_enthalpy_res(m, V, T, x), "cv_res" => cv,
                 "lnphi" => Z > 0 ? mu .- log(Z) : nothing, "a_res" => CP.eos_res(m, V, T, x) / (RG * T))
        end
    elseif op == "t1-derivs"
        r = T1[req["state"]]; m = t1_model(r["model"]); y = Float64[r["T"], r["rho"]]
        f(v) = CP.eos_res(m, 1 / v[2], v[1], [1.0]) / (RG * v[1])
        order = req["order"]
        return () -> order == "Value" ? Dict("value" => f(y)) :
                     order == "Gradient" ? Dict("value" => f(y), "gradient" => ForwardDiff.gradient(f, y)) :
                     Dict("value" => f(y), "gradient" => ForwardDiff.gradient(f, y), "hessian" => vec(ForwardDiff.hessian(f, y)))
    elseif op == "t2-sequence"
        case = only(c for c in DENSITY["cases"] if c["id"] == req["case"])
        return () -> begin
            out = measure_case(DENSITY, case)
            unique = [s for s in out["samples"] if s["mode"] == PRIMARY_MODE]
            Dict("seconds" => sum(something(s["elapsed_us"], 0.0) for s in unique) / 1e6, "per_call" => length(unique), "value" => unique)
        end
    elseif op == "t2-loop"
        loop = only(c for c in INPUTS["t2"]["ethane_loop"] if c["id"] == req["case"]); m = audit_model("ethane")
        T = Float64(loop["T"]); phase = Symbol(loop["branch"])
        return () -> begin
            rows = []; ns = 0
            for P in Float64.(loop["pressures"])
                t = time_ns(); v = volume(m, P, T, [1.0]; phase=phase); ns += time_ns() - t
                push!(rows, isfinite(v) ? Dict("rho" => 1 / v, "P" => pressure(m, v, T, [1.0])) : Dict("error" => "volume $v"))
            end
            Dict("seconds" => ns / 1e9, "per_call" => length(rows), "value" => rows)
        end
    elseif op == "t3"
        name = req["case"]
        if name == "methane-bubble"
            spec = DENSITY["comparator_inputs"]["models"]["methane-ethane"]
            q = Dict("case" => "methane-ethane", "x" => [0.5, 0.5], "T" => 180.0, "operation" => "bubble")
        else
            q = only(r for r in CROSS["requests"] if r["id"] == name); spec = CROSS["models"][q["case"]]
        end
        m, _ = build_model(spec, q, CROSS["constants"])
        return () -> call_block(m, q)
    elseif op == "t4"
        sys = req["system"]; z = [0.5, 0.5]
        m = sys == "methane_ethane" ? audit_model("methane-ethane") : audit_model("methanol-hexane"; kij=INPUTS["t4"]["methanol_hexane"]["kij"])
        T, P = Float64(req["T"]), Float64(req["P"])
        method = req["method"] == "michelsen" ? MichelsenTPFlash() : nothing
        return () -> begin
            f = isnothing(method) ? CP.tp_flash2(m, P, T, z) : CP.tp_flash2(m, P, T, z, method)
            Dict("compositions" => f.compositions, "volumes" => f.volumes, "fractions" => f.fractions,
                 "phases" => length(f.compositions))
        end
    elseif op == "t6-curve"
        q0 = only(r for r in CROSS["requests"] if r["id"] == INPUTS["t6"]["request"] * "-1")
        m, _ = build_model(CROSS["models"]["C3"], q0, CROSS["constants"])
        ms = Float64.(INPUTS["t6"]["molalities"])
        return () -> [call_block(m, merge(q0, Dict("m" => mi))) for mi in ms]
    end
    error("unknown op $op")
end

finite(x) = x isa AbstractFloat ? (isfinite(x) ? x : nothing) : x isa AbstractDict ? Dict(k => finite(v) for (k, v) in x) :
            x isa AbstractVector ? finite.(x) : x isa Tuple ? finite.(collect(x)) : x

function serve(req)
    key = JSON.json(Dict(k => v for (k, v) in req if !(k in ("samples", "warmups", "out", "block"))))
    prepare = nothing
    if !haskey(PREPARED, key)
        t = time_ns(); PREPARED[key] = make(req); prepare = (time_ns() - t) / 1e9
    end
    f = PREPARED[key]
    for _ in 1:get(req, "warmups", 0); Base.invokelatest(f); end
    samples = []
    for _ in 1:get(req, "samples", 1)
        t = time_ns()
        value = try
            Base.invokelatest(f)
        catch e
            Dict("error" => sprint(showerror, e))
        end
        elapsed = (time_ns() - t) / 1e9
        push!(samples, value isa AbstractDict && haskey(value, "seconds") ? value : Dict("value" => value, "seconds" => elapsed))
    end
    Dict("samples" => samples, "prepare_seconds" => prepare)
end

println("READY ", JSON.json(Dict("pid" => getpid(), "version" => string(pkgversion(Clapeyron)), "julia" => string(VERSION),
                                 "threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads())))
flush(stdout)
while true
    for line in eachline(ARGS[1])
        path = strip(line); path == "quit" && exit(); isempty(path) && continue
        req = JSON.parsefile(path)
        reply = try
            serve(req)
        catch e
            Dict("error" => sprint(showerror, e, catch_backtrace()))
        end
        open(io -> JSON.print(io, finite(json_numbers(reply))), req["out"] * ".tmp", "w")
        mv(req["out"] * ".tmp", req["out"]; force=true)
        println("DONE ", path); flush(stdout)
    end
end
