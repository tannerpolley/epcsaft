let spec = H2H_IN["clapeyron_khudaida"], model = khudaida_model(spec)
    z = Float64.(spec["feeds"]["1"])
    for m in (MichelsenTPFlash(equilibrium=:lle, K_tol=1e-12, ss_iters=1000, nacc=5, second_order=true, x0=[0.80,0.05,0.05,0.05,0.05], y0=[0.20,0.10,0.6998,0.0001,0.0001]),
              MichelsenTPFlash(equilibrium=:lle, K_tol=1e-12, ss_iters=1000, nacc=5, second_order=false, x0=[0.80,0.05,0.05,0.05,0.05], y0=[0.20,0.10,0.6998,0.0001,0.0001]),
              MichelsenTPFlash(equilibrium=:lle, K_tol=1e-12, ss_iters=1000, nacc=5, second_order=false))
        try
            t = time_ns(); f = CP.tp_flash2(model, 1e5, 293.15, z, m)
            println("probe2 ", (time_ns()-t)/1e9, " s ", f.compositions, " ", f.fractions)
        catch e
            println("probe2 error ", sprint(showerror, e, catch_backtrace())[1:min(end,2000)])
        end
    end
end
