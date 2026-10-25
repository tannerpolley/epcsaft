let model = khudaida_model(H2H_IN["clapeyron_khudaida"]), rows = JSON.parsefile(joinpath(H2H_ROOT, "build/h2h/probe-perturbed.json"))
    for r in rows
        x = Float64.(r["x"]); V = 1/Float64(r["rho"]); id = CP.iondata(model, V, 293.15, x)
        println("probe6 eps ", id[3], " vs ", r["eps"], " | dh ", CP.a_dh(model.ionmodel, V, 293.15, x, id), " vs ", r["debye_huckel"], " | born ", CP.a_born(model.ionmodel, V, 293.15, x, id), " vs ", r["born"])
    end
end
