# One long-lived Julia process for the 2026-10-08 head-to-head: compile Clapeyron once, then include each
# task file whose path arrives on the FIFO in ARGS[1] ("quit" ends it).
include(joinpath(@__DIR__, "..", "..", "2026-competitor-cross-check", "scripts", "clapeyron.jl"))  # build_model, ASSOC, BUB, CP; its own entry points stay inactive here
println("server import+include ", IMPORT_NS/1e9, " s"); flush(stdout)

# Engine packet permittivity (Figiel Eq. 11-12): salt-free mass-weighted solvent permittivity
# suppressed by the total ion mole fraction. Clapeyron ships no such RSP model; this user subtype of its
# RSPModel interface supplies only the scalar relation, and Clapeyron differentiates and solves with it.
struct SaltFreeSuppressedRSP <: CP.RSPModel
    components::Vector{String}
    params::NamedTuple{(:epsilon, :Mw, :charge), Tuple{CP.SingleParam{Float64}, CP.SingleParam{Float64}, CP.SingleParam{Float64}}}
    C::Float64
end
function CP.dielectric_constant(m::SaltFreeSuppressedRSP, V, T, z, Z = nothing)
    q = m.params.charge.values; eps = m.params.epsilon.values; mw = m.params.Mw.values
    num = zero(eltype(z)); den = zero(eltype(z)); ion = zero(eltype(z))
    for i in eachindex(z)
        if iszero(q[i]); num += z[i]*mw[i]*eps[i]; den += z[i]*mw[i]; else; ion += z[i]; end
    end
    return (num/den)/(1 + m.C*ion/sum(z))
end
CP.is_splittable(::SaltFreeSuppressedRSP) = false

while true
    for line in eachline(ARGS[1])
        path = strip(line)
        path == "quit" && exit()
        isempty(path) && continue
        t = time_ns()
        try
            Base.include(Main, path)
        catch e
            println(stderr, sprint(showerror, e, catch_backtrace()))
        end
        println("done ", path, " ", (time_ns()-t)/1e9, " s"); flush(stdout); flush(stderr)
    end
end
