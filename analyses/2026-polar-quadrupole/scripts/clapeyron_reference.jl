# Clapeyron 0.6.29 QPCPSAFT comparator rows for the Gross (2005) quadrupole term (#574), in one Julia process:
#   julia --startup-file=no --project=analyses/2026-competitor-cross-check/inputs analyses/2026-polar-quadrupole/scripts/clapeyron_reference.jl
# Gross (2005) Table A1 sets, k_ij = 0. Clapeyron converts D A with exactly 1e-19 J A^5; the Engine uses the CODATA 2018
# vacuum permittivity, so each moment is multiplied by sqrt(r), r = Engine (D A)^2/(4 pi eps0) / 1e-19, to give
# identical reduced moments. Clapeyron applies no chain-length cap to the quadrupole. Span-Wagner CO2 is SingleFluid.
# Saturation and bubble solves use the 1e-12 tolerances of analyses/2026-competitor-cross-check (the defaults stop near 1e-9).
using Clapeyron, Printf
const C = Clapeyron
debye = 1e-21 / 299792458.0
const r = debye^2 / (4pi * 8.8541878128e-12 * 1e-30) / 1e-19
const s = sqrt(r)
const SAT = ChemPotVSaturation(atol=1e-12, rtol=1e-12, max_iters=1000)
const BUB = ChemPotBubblePressure(atol=1e-12, rtol=1e-12, max_iters=1000)
model(names, Mw, m, sigma, eps, Q) = QPCPSAFT(names; userlocations=(; Mw, segment=m, sigma, epsilon=eps,
    k=zeros(length(names), length(names)), dipole=zeros(length(names)), quadrupole=s .* Q,
    n_H=zeros(Int, length(names)), n_e=zeros(Int, length(names)), epsilon_assoc=nothing, bondvol=nothing))
co2 = model(["carbon-dioxide"], [44.01], [1.5131], [3.1869], [163.33], [4.4])
n2 = (["carbon-dioxide", "nitrogen"], [44.01, 28.01], [1.5131, 1.1504], [3.1869, 3.3848], [163.33, 91.40], [4.4, 1.43])
co2n2 = model(n2...)
benzene = model(["benzene"], [78.114], [2.2463], [3.7852], [296.24], [5.5907])
out = joinpath(@__DIR__, "..", "results")
@printf(stderr, "moment ratio r - 1 = %.6e\n", r - 1)
open(joinpath(out, "clapeyron-term-rows.csv"), "w") do io
    println(io, "species,x,T_K,rho_mol_m3,a_qq_clapeyron")
    for (names, mdl, x, T, rho) in (("carbon-dioxide", co2, [1.0], 250.0, 1000.0), ("carbon-dioxide", co2, [1.0], 250.0, 22000.0),
                                    ("carbon-dioxide+nitrogen", co2n2, [0.7, 0.3], 240.0, 15000.0),
                                    ("benzene", benzene, [1.0], 400.0, 9000.0))
        @printf(io, "%s,%s,%.17g,%.17g,%.17g\n", names, join(x, " "), T, rho, C.a_qq(mdl, 1 / rho, T, x))
    end
end
sw = SingleFluid("carbon dioxide")
open(joinpath(out, "clapeyron-saturation-rows.csv"), "w") do io
    println(io, "species,T_K,p_sat_clapeyron_Pa,p_sat_span_wagner_Pa")
    for T in 220.0:5.0:300.0
        @printf(io, "carbon-dioxide,%.17g,%.17g,%.17g\n", T, saturation_pressure(co2, T, SAT)[1], saturation_pressure(sw, T)[1])
    end
    for T in 290.0:10.0:540.0
        @printf(io, "benzene,%.17g,%.17g,\n", T, saturation_pressure(benzene, T, SAT)[1])
    end
end
open(joinpath(out, "clapeyron-bubble-rows.csv"), "w") do io
    println(io, "T_K,x_nitrogen,p_bubble_clapeyron_Pa,y_nitrogen_clapeyron")
    for T in (240.0, 270.0)
        p, _, _, y = bubble_pressure(co2n2, T, [0.95, 0.05], BUB)
        @printf(io, "%.17g,0.05,%.17g,%.17g\n", T, p, y[2])
    end
end
