Base.eltype(::SaltFreeSuppressedRSP) = Float64
let spec = H2H_IN["clapeyron_khudaida"], model = khudaida_model(spec)
    z = Float64.(spec["feeds"]["1"])
    try
        f = CP.tp_flash2(model, 1e5, 293.15, z, MichelsenTPFlash(equilibrium=:lle, K_tol=1e-12, ss_iters=1000, nacc=5, second_order=true))
        println("probe ok ", f.compositions)
    catch e
        println("probe error ", sprint(showerror, e, catch_backtrace())[1:min(end,3000)])
    end
end
