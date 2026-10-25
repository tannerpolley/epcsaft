let spec = H2H_IN["clapeyron_khudaida"], model = khudaida_model(spec), eng = H2H_IN["khudaida"]["matched"][1]
    z = Float64.(spec["feeds"]["1"]); T = 293.15; P = 1e5
    m = MichelsenTPFlash(equilibrium=:lle, K_tol=1e-12, ss_iters=1000, nacc=5, second_order=true, x0=[0.80,0.05,0.05,0.05,0.05], y0=[0.20,0.10,0.6998,0.0001,0.0001])
    f = CP.tp_flash2(model, P, T, z, m)
    println("probe4 clapeyron-solution ", lle_residuals(model, P, T, f), " phases ", f.compositions)
    xs = [Float64.(x) for x in eng["x"]]; vs = [volume(model, P, T, x; phase=:liquid, vol0=1/r) for (x,r) in zip(xs, eng["rho"])]
    fe = CP.FlashResult(xs, [eng["beta"], 1-eng["beta"]], vs, CP.FlashData(P, T))
    println("probe4 engine-solution-in-clapeyron ", lle_residuals(model, P, T, fe))
    μ = [CP.VT_chemical_potential(model, v, T, x)./(CP.Rgas(model)*T) for (v,x) in zip(f.volumes, f.compositions)]
    println("probe4 clap mu diff ", μ[1] .- μ[2])
    println("probe4 gibbs clap ", sum(f.fractions[i]*dot(f.compositions[i], μ[i]) for i in 1:2))
    μe = [CP.VT_chemical_potential(model, v, T, x)./(CP.Rgas(model)*T) for (v,x) in zip(vs, xs)]
    println("probe4 engine mu diff ", μe[1] .- μe[2])
    println("probe4 gibbs engine ", eng["beta"]*dot(xs[1], μe[1]) + (1-eng["beta"])*dot(xs[2], μe[2]))
end
