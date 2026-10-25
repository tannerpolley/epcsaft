let model = khudaida_model(H2H_IN["clapeyron_khudaida"]), r = JSON.parsefile(joinpath(H2H_ROOT, "build/h2h/probe-engine-terms.json"))[1]
    x = Float64.(r["x"]); V = 1/Float64(r["rho"]); T = 293.15
    fdh(z) = sum(z)*CP.a_dh(model.ionmodel, V, T, z, CP.iondata(model, V, T, z))
    feps(z) = CP.iondata(model, V, T, z)[3]
    ad = ForwardDiff.gradient(fdh, x); ade = ForwardDiff.gradient(feps, x)
    fd = [(fdh(x .+ 1e-6 .* (1:5 .== i)) - fdh(x .- 1e-6 .* (1:5 .== i)))/2e-6 for i in 1:5]
    fde = [(feps(x .+ 1e-6 .* (1:5 .== i)) - feps(x .- 1e-6 .* (1:5 .== i)))/2e-6 for i in 1:5]
    println("probe7 dh AD ", ad); println("probe7 dh FD ", fd); println("probe7 eps AD ", ade); println("probe7 eps FD ", fde)
    κ(z) = CP.debye_length(V, T, z, CP.iondata(model, V, T, z)[3], model.charge)
    println("probe7 kappa AD ", ForwardDiff.gradient(κ, x))
end
