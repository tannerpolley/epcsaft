let model = khudaida_model(H2H_IN["clapeyron_khudaida"]), rows = JSON.parsefile(joinpath(H2H_ROOT, "build/h2h/probe-engine-terms.json"))
    for r in rows
        x = Float64.(r["x"]); V = 1/Float64(r["rho"]); T = 293.15
        nm = model.neutralmodel; d = CP.data(nm, V, T, x); id = CP.iondata(model, V, T, x)
        c = Dict("hard_chain"=>CP.a_hc(nm,V,T,x,d), "dispersion"=>CP.a_disp(nm,V,T,x,d), "association"=>CP.a_assoc(nm,V,T,x,d),
                 "debye_huckel"=>CP.a_dh(model.ionmodel,V,T,x,id), "born"=>CP.a_born(model.ionmodel,V,T,x,id),
                 "bulk_relative_permittivity"=>id[3], "residual_helmholtz"=>CP.a_res(model,V,T,x), "pressure"=>pressure(model,V,T,x))
        for k in sort(collect(keys(c)))
            println("probe3 ", k, " clap=", c[k], " engine=", r[k], " diff=", c[k]-r[k])
        end
        println("probe3 dh_sigma ", id[2], " hs d ", CP.d(nm,V,T,x))
    end
end
