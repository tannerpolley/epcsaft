let model = khudaida_model(H2H_IN["clapeyron_khudaida"]), rows = JSON.parsefile(joinpath(H2H_ROOT, "build/h2h/probe-engine-terms.json")), out = []
    nm = model.neutralmodel; T = 293.15
    terms = Dict("hard_chain"=>(V,z)->CP.a_hc(nm,V,T,z,CP.data(nm,V,T,z)), "dispersion"=>(V,z)->CP.a_disp(nm,V,T,z,CP.data(nm,V,T,z)),
        "association"=>(V,z)->CP.a_assoc(nm,V,T,z,CP.data(nm,V,T,z)),
        "debye_huckel"=>(V,z)->CP.a_dh(model.ionmodel,V,T,z,CP.iondata(model,V,T,z)), "born"=>(V,z)->CP.a_born(model.ionmodel,V,T,z,CP.iondata(model,V,T,z)))
    for r in rows
        x = Float64.(r["x"]); V = 1/Float64(r["rho"]); d = Dict{String,Any}()
        for (k,f) in terms
            d[k] = ForwardDiff.gradient(z -> sum(z)*f(V,z), x)
        end
        push!(out, d)
    end
    open(io -> JSON.print(io, out), joinpath(H2H_ROOT, "build/h2h/probe-clapeyron-mu-terms.json"), "w")
end
