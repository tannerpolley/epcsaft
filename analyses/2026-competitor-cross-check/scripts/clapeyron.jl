const IMPORT_START = time_ns()
using JSON, LinearAlgebra
BLAS.set_num_threads(1)

# full-eos §10: NaCl, MW in kg/mol; ionic reference and pure water are distinct.
function contract(q)
    l, r, p = q["lnphi"], q["reference_lnphi"], q["pure_lnphi"]
    gx = exp(((l[2]-r[2])+(l[3]-r[3]))/2)
    pw = p isa AbstractVector ? p[1] : p
    return Dict("gamma_x"=>gx, "gamma"=>gx*q["xwater"],
                "phi"=>-(log(q["xwater"])+l[1]-pw)/(2*q["m"]*q["MW"]))
end
if !isempty(ARGS) && ARGS[1] == "--contract"
    for line in eachline(stdin)
        try
            println(JSON.json(contract(JSON.parse(line))))
        catch e
            println(JSON.json(Dict("error"=>sprint(showerror,e))))
        end
        flush(stdout)
    end
    exit()
end
using Clapeyron, ForwardDiff
const IMPORT_NS = time_ns()-IMPORT_START
const CP = Clapeyron
const ASSOC = AssocOptions(atol=1e-13,rtol=1e-13,max_iters=1000,combining=:nocombining)
const SAT = ChemPotVSaturation(atol=1e-12,rtol=1e-12,max_iters=1000)
const BUB = ChemPotBubblePressure(atol=1e-12,rtol=1e-12,max_iters=1000)

function build_model(spec, q, constants; assoc=ASSOC)
    cs = spec["components"]
    names = String[c["id"] for c in cs]
    n = length(cs)
    epsab = Dict{Tuple{Tuple{String,String},Tuple{String,String}},Float64}()
    kap = copy(epsab)
    for name in names
        epsab[((name,"H"),(name,"e"))] = 0.
        kap[((name,"H"),(name,"e"))] = 0.
    end
    factors = []
    for i in 1:n, j in i:n
        ci,cj = cs[i],cs[j]
        ki,kj = get(ci,"kappa_ab",0.),get(cj,"kappa_ab",0.)
        ki*kj == 0 && continue
        sij = (ci["sigma"]+cj["sigma"])/2
        kij = sqrt(ki*kj)*(sqrt(ci["sigma"]*cj["sigma"])/sij)^3
        eij = (ci["epsilon_ab"]+cj["epsilon_ab"])/2
        for (a,b) in (("H","e"),("e","H"))
            epsab[((names[i],a),(names[j],b))] = eij
            kap[((names[i],a),(names[j],b))] = kij
        end
        source = sqrt(ki*kj)*(ci["sigma"]*cj["sigma"])^1.5
        push!(factors,Dict("i"=>i,"j"=>j,"epsilon_ab"=>eij,"kappa_ab"=>kij,
            "entered_volume_factor"=>sij^3*kij,"source_volume_factor"=>source,
            "relative_difference"=>(sij^3*kij-source)/source))
    end
    k = zeros(n,n)
    if n == 2
        k[1,2] = k[2,1] = spec["k"]
    end
    for (i,j,v) in get(spec,"pairs",[])
        k[Int(i)+1,Int(j)+1] = k[Int(j)+1,Int(i)+1] = v
    end
    if q["case"] == "C3"
        k[1,3] = k[3,1] = q["kcl"]
    end
    sites = [get(c,"kappa_ab",0.) > 0 ? 1 : 0 for c in cs]
    pars = (Mw=Float64[1000*c["molar_mass"] for c in cs], segment=Float64[c["m"] for c in cs],
        sigma=Float64[c["sigma"] for c in cs], epsilon=Float64[c["epsilon"] for c in cs],
        n_H=sites,n_e=sites,k=k,epsilon_assoc=epsab,bondvol=kap)
    if q["case"] != "C3"
        model = q["case"] == "C5" ? PCPSAFT(names; userlocations=merge(pars,
            (dipole=Float64[c["mu"]*sqrt(c["nmu"])*(get(q,"dipole_convention","")=="matched" ?
                sqrt((1e-21/299792458.)^2/(4pi*constants["Engine"]["epsilon0"]*1e-30*1e-19)) : 1.) for c in cs],)),assoc_options=assoc) :
            PCSAFT(names;userlocations=pars,assoc_options=assoc)
        return model,factors
    end
    # ePCSAFT.jl constructor: suppress like charges after mixed energies form.
    neutral = pharmaPCSAFT(names;userlocations=pars,assoc_options=assoc)
    charges = Int[get(c,"z",0) for c in cs]
    for i in 1:n, j in 1:n
        if charges[i]*charges[j] > 0
            neutral.params.epsilon.values[i,j] = 0.
        end
    end
    ce,cc = constants["Engine"],constants["Clapeyron"]
    er = q["electrostatic"] == "matched" ?
        78.09*(cc["e"]/ce["e"])^2*(ce["epsilon0"]/cc["epsilon0"])*(ce["kB"]/cc["kB"]) : 78.09
    ion = (q["diameter"] == "sigma" ? DH : hsdDH)(names[1:1],names[2:3];
        charge=charges,RSPmodel=ConstRSP(er))
    model = CP.ePCSAFT(names,charges,BasicIdeal(names),neutral,ion,["10.1016/j.cherd.2014.05.017"])
    return model,factors
end

composition(model,m) = molality_to_composition(model,CP.auto_binary_salts(model),[Float64(m)],[1.])
function fixed_values(model,T,rho,x)
    v = inv(rho)
    a = CP.a_res(model,v,T,x)
    p = pressure(model,v,T,x)
    Z = p*v/(CP.Rgas(model)*T)
    mu = CP.VT_chemical_potential_res(model,v,T,x)/(CP.Rgas(model)*T)
    return Dict("rho"=>rho,"x"=>x,"pressure"=>p,"ar"=>a,"Z"=>Z,
                "lnphi"=>Z>0 ? mu .- log(Z) : nothing,"mu_res_over_RT"=>mu)
end
function phase_values(model,T,rho,x,label)
    d = fixed_values(model,T,rho,x)
    d["phase"] = label
    v = inv(rho)
    p,dpdv = CP.p∂p∂V(model,v,T,x)
    d["pressure_density_derivative"] = -dpdv/rho^2
    neutral = model isa CP.ESElectrolyteModel ? model.neutralmodel : model
    data = CP.data(neutral,v,T,x)
    d["contributions"] = Dict("hard_chain"=>CP.a_hc(neutral,v,T,x,data),
        "dispersion"=>CP.a_disp(neutral,v,T,x,data),"association"=>CP.a_assoc(neutral,v,T,x,data),
        "debye_huckel"=>model isa CP.ESElectrolyteModel ? CP.a_res(model.ionmodel,v,T,x,CP.iondata(model,v,T,x)) : 0.,
        "dipole"=>neutral isa CP.PCPSAFTModel ? CP.a_polar(neutral,v,T,x,data) : 0.,"born"=>0.)
    if model isa CP.ESElectrolyteModel
        z,diam,er = CP.iondata(model,v,T,x)
        kap = CP.debye_length(v,T,x,er,z)
        d["kappa_d"] = kap .* diam
        d["screening_diameter"] = diam
        d["relative_permittivity"] = er
        d["bjerrum_length"] = CP.e_c^2/(4pi*CP.ϵ_0*er*CP.k_B*T)
        d["sigma_water"] = neutral.params.sigma.values[1,1] + CP.Δσh20(T)
    end
    return d
end

# Timed extraction excludes qualification contributions and source diagnostics.
function call_block(model,q)
    T,op = Float64(q["T"]),String(q["operation"])
    if op == "fixed"
        return fixed_values(model,T,q["rho"],Float64.(q["x"]))
    elseif op == "saturation"
        p,vl,vv = Base.invokelatest(saturation_pressure,model,T,SAT)
        return Dict("P"=>p,"densities"=>[1/vl,1/vv],"compositions"=>[[1.],[1.]])
    elseif op == "bubble"
        x = Float64.(q["x"])
        p,vl,vv,y = Base.invokelatest(bubble_pressure,model,T,x,BUB)
        return Dict("P"=>p,"densities"=>[1/vl,1/vv],"compositions"=>[x,collect(y)])
    elseif op == "lle"
        z = Float64.(q["x"])
        method = MichelsenTPFlash(equilibrium=:lle,K_tol=1e-12,ss_iters=1000,nacc=5,second_order=true)
        f = Base.invokelatest(CP.tp_flash2,model,Float64(q["P"]),T,z,method)
        ix = sortperm([x[1] for x in f.compositions];rev=true)
        return Dict("P"=>q["P"],"densities"=>1 ./ f.volumes[ix],
            "compositions"=>f.compositions[ix],"beta"=>f.fractions[ix[1]]/sum(f.fractions))
    elseif op == "activity"
        g = Base.invokelatest(mean_ionic_activity_coefficient,model,Float64(q["P"]),T,[Float64(q["m"])])
        phi = Base.invokelatest(osmotic_coefficient,model,Float64(q["P"]),T,[Float64(q["m"])])
        return Dict("gamma"=>only(g),"phi"=>only(phi),"P"=>q["P"])
    elseif op == "reference"
        x = composition(model,q["m"])
        v = volume(model,q["P"],T,x;phase=:liquid,threaded=false)
        return Dict("P"=>q["P"],"densities"=>[1/v],"compositions"=>[x])
    end
    error("Unsupported operation: $op")
end
function qualify(model,q,d)
    op,T = q["operation"],q["T"]
    out = Dict("id"=>q["id"],"case"=>q["case"],"T"=>T)
    merge!(out,d)
    if op == "activity"
        x = composition(model,q["m"])
        rho = 1/volume(model,q["P"],T,x;phase=:liquid,threaded=false)
        out["phases"] = [phase_values(model,T,rho,x,"liquid")]
        refs = []
        for m in (1e-20,0.)
            xr = composition(model,m)
            rr = 1/volume(model,q["P"],T,xr;phase=:liquid,threaded=false)
            push!(refs,merge(phase_values(model,T,rr,xr,"liquid"),Dict("m"=>m)))
        end
        out["reference_phases"] = refs
        out["reference_m"] = 1e-20
    elseif op == "fixed"
        out["phases"] = [phase_values(model,T,q["rho"],Float64.(q["x"]),get(q,"phase","liquid"))]
    else
        labels = op == "lle" ? ["water-rich","water-poor"] :
            op == "reference" ? ["liquid"] : ["liquid","vapor"]
        out["phases"] = [phase_values(model,T,r,x,l) for (r,x,l) in zip(d["densities"],d["compositions"],labels)]
    end
    for ph in out["phases"]
        if ph["lnphi"] === nothing
            ph["lnphi_absence"] = "Nonpositive pressure/Z at this fixed state"
        end
        all(isfinite,[ph["rho"],ph["ar"],ph["Z"],ph["pressure"]]) || error("Nonfinite phase result")
    end
    return out
end

# DH.jl published rational coefficients, lifted from their actual Float64 values.
function f_analytic(y)
    if abs(y) <= 0.5
        s,t,k = one(y)/3,one(y),0
        while true
            k += 1
            t *= -y
            term = t/(k+3)
            s += term
            k >= 3 && abs(term) < big"1e-60" && return s
        end
    end
    return (log1p(y)-y+y*y/2)/(y*y*y)
end
function f_source(y)
    if y <= 1.5
        p = BigFloat.((0.3333333333333333,0.7222222222346917,0.5314393939617089,0.15151515152716674,0.013227513229354138))
        q = BigFloat.((1.,2.9166666667040753,3.181818181913183,1.5909090909939427,0.3535353535662185,0.026515151518857805))
        return evalpoly(y,p)/evalpoly(y,q)
    end
    return f_analytic(y) # Above the rational branch the two equations are identical.
end
f_feos(y) = inv(y^3)*(log(1+y)-2*(1+y)+(1+y)^2/2+1.5)
function source_functions(ys)
    rows = []
    for y in ys, (name,f) in (("analytic",f_analytic),("library_expression",f_source))
        d1 = ForwardDiff.derivative(f,y)
        d2 = ForwardDiff.derivative(t->ForwardDiff.derivative(f,t),y)
        push!(rows,Dict("y"=>string(y),"function"=>name,"f"=>string(f(y)),"f_prime"=>string(d1),"f_second"=>string(d2)))
    end
    return rows
end
function source_correction(model,T,rho,x)
    _,diam,er = CP.iondata(model,inv(rho),T,x)
    if iszero(x[2]+x[3])
        return Dict("delta_a"=>big"0","delta_Z"=>big"0","delta_mu"=>zeros(BigFloat,3),
            "delta_lnphi_fixed"=>zeros(BigFloat,3),"delta_P"=>big"0",
            "delta_P_density_derivative"=>big"0","kappa_d"=>zeros(BigFloat,3),
            "source_function"=>source_functions([big"0"]))
    end
    lb = BigFloat(CP.e_c)^2/(4big(pi)*BigFloat(CP.ϵ_0)*BigFloat(er)*BigFloat(CP.k_B)*T)
    na,rt = BigFloat(CP.N_A),BigFloat(CP.Rgas(model))*T
    function da(u)
        r,n = u[1],u[2:end]
        xx = n/sum(n)
        kap = sqrt(4big(pi)*lb*na*r*(xx[2]+xx[3]))
        return -lb*kap*sum(xx[i]*(f_source(kap*diam[i])-f_analytic(kap*diam[i])) for i in 2:3)
    end
    u = vcat(rho,x)
    g = ForwardDiff.gradient(da,u)
    h = ForwardDiff.hessian(da,u)
    dz = rho*g[1]
    dm = da(u) .+ dz .+ g[2:end] .- dot(x,g[2:end])
    raw = fixed_values(model,T,rho,x)
    dl = dm .- log(raw["Z"]/(raw["Z"]-dz))
    dp = rho*rt*dz
    slope = rt*(2rho*g[1]+rho^2*h[1,1])
    kap = sqrt(4big(pi)*lb*na*rho*(x[2]+x[3]))
    return Dict("delta_a"=>da(u),"delta_Z"=>dz,"delta_mu"=>dm,"delta_lnphi_fixed"=>dl,
        "delta_P"=>dp,"delta_P_density_derivative"=>slope,"kappa_d"=>kap .* diam,
        "source_function"=>source_functions(vcat(kap .* diam[2:3],[big"1.5"-big"1e-12",big"1.5"+big"1e-12"])))
end
function diagnose_phase(model,q,ph,bits)
    setprecision(BigFloat,bits) do
        T,rho,x = BigFloat(q["T"]),BigFloat(ph["rho"]),BigFloat.(ph["x"])
        x /= sum(x)
        # First close the raw library root in the same precision; preserve its shift.
        supplied = rho
        if q["operation"] != "fixed"
            for _ in 1:2
                st = phase_values(model,T,rho,x,ph["phase"])
                rho -= (st["pressure"]-BigFloat(q["P"]))/st["pressure_density_derivative"]
            end
        end
        rawroot = rho
        initial = fixed_values(model,T,rho,x)
        rows = []
        for step in 0:(q["operation"] == "fixed" ? 0 : 2)
            raw = phase_values(model,T,rho,x,ph["phase"])
            cor = source_correction(model,T,rho,x)
            analytic = raw["lnphi"] - cor["delta_lnphi_fixed"]
            displacement = -(raw["pressure"]-cor["delta_P"]-BigFloat(q["P"]))/
                (raw["pressure_density_derivative"]-cor["delta_P_density_derivative"])
            iszero(cor["delta_P"]) && (displacement = zero(rho))
            push!(rows,Dict("precision_bits"=>bits,"refinement"=>step,"rho"=>string(rho),
                "raw"=>raw,"correction"=>cor,"analytic_lnphi"=>analytic,
                "complete_signed_lnphi_error"=>initial["lnphi"]-analytic,"raw_root_density_shift"=>string(rawroot-supplied),
                "density_displacement"=>string(displacement)))
            rho += displacement
        end
        return rows
    end
end

function diagnose_request(model,q,r,mw)
    phases = vcat(r["phases"],get(r,"reference_phases",[]))
    bybits,rows = Dict(),[]
    for bits in (256,512)
        allph = [Base.invokelatest(diagnose_phase,model,q,ph,bits) for ph in phases]
        bybits[bits] = allph
        for (i,rs) in enumerate(allph), s in rs
            push!(rows,Dict("id"=>q["id"],"case"=>"C3","phase_index"=>i,"diagnostic"=>s))
        end
        if q["operation"] == "activity"
            setprecision(BigFloat,bits) do
                xw = allph[1][1]["raw"]["x"][1]
                rawlogs = [p[1]["raw"]["lnphi"] for p in allph]
                function coefficients(logs)
                    contract(Dict("lnphi"=>logs[1],"reference_lnphi"=>logs[2],
                        "pure_lnphi"=>logs[3],"m"=>BigFloat(q["m"]),"MW"=>BigFloat(mw),"xwater"=>xw))
                end
                original = coefficients(rawlogs)
                for step in 1:3
                    logs = [p[step]["analytic_lnphi"] for p in allph]
                    corrected = coefficients(logs)
                    push!(rows,Dict("id"=>q["id"],"case"=>"C3","precision_bits"=>bits,
                        "refinement"=>step-1,"corrected_gamma"=>corrected["gamma"],"corrected_phi"=>corrected["phi"],
                        "delta_ln_gamma"=>log(original["gamma"])-log(corrected["gamma"]),
                        "delta_phi"=>original["phi"]-corrected["phi"],
                        "native_helper_minus_refined_gamma"=>r["gamma"]-original["gamma"],
                        "native_helper_minus_refined_phi"=>r["phi"]-original["phi"]))
                end
            end
        end
    end
    for (aa,bb) in zip(bybits[256],bybits[512]), (a,b) in zip(aa,bb)
        for key in ("delta_a","delta_Z","delta_lnphi_fixed")
            va,vb = a["correction"][key],b["correction"][key]
            err,mag = maximum(abs.(va.-vb)),maximum(abs.(vb))
            limit = mag == 0 ? big"1e-13" : min(big"1e-13",mag*big"1e-3")
            err <= limit || error("256/512-bit correction not stable: $key")
        end
        va,vb = a["complete_signed_lnphi_error"],b["complete_signed_lnphi_error"]
        err,mag = maximum(abs.(va.-vb)),maximum(abs.(vb))
        err <= (mag==0 ? big"1e-13" : min(big"1e-13",mag*big"1e-3")) || error("256/512-bit complete correction not stable")
    end
    return rows
end

# Table values use Float64 transport; precision checks precede conversion.
function json_numbers(x)
    x isa BigFloat && return Float64(x)
    x isa AbstractDict && return Dict(k=>json_numbers(v) for (k,v) in x)
    x isa AbstractVector && return json_numbers.(x)
    return x
end
function association_probe(input,native)
    rows=[]
    fixed=[q for q in input["requests"] if q["case"]=="C4" && q["operation"]=="fixed"]
    for (tag,tol,implicit) in (("baseline",1e-13,false),("implicit",1e-13,true),("tighter",1e-30,true))
        options=AssocOptions(atol=tol,rtol=tol,max_iters=1000,combining=:nocombining,implicit_ad=implicit)
        model,_=build_model(input["models"]["C4"],first(fixed),input["constants"];assoc=options)
        @assert model.assoc_options.implicit_ad==implicit
        for q in fixed, bits in (53,256)
            d=setprecision(BigFloat,256) do
                scalar=bits==53 ? Float64 : BigFloat
                T,rho,x=scalar(q["T"]),scalar(q["rho"]),scalar.(q["x"])
                x/=sum(x)
                values=phase_values(model,T,rho,x,get(q,"phase","liquid"))
                K=CP.assoc_site_matrix(model,inv(rho),T,x)
                J,_=CP.compress_assoc_matrix(K)
                values["matrix_dimension"]=size(K,1)
                values["compressed_dimension"]=size(J,1)
                values["compression_selected"]=CP.__maybe_compress(K)
                values
            end
            push!(rows,Dict("id"=>q["id"]*":$tag:$bits","precision_bits"=>bits,"tolerance"=>tol,"implicit_ad"=>implicit,"native"=>d))
        end
        tag=="baseline" || continue
        for r in native["records"]
            r["case"]=="C4" && haskey(r,"beta") || continue
            q=only(q for q in input["requests"] if q["id"]==r["id"])
            for p in r["phases"]
                T,P,x=Float64(q["T"]),Float64(q["P"]),Float64.(p["x"])
                v=volume(model,P,T,x;phase=:liquid,vol0=inv(p["rho"]))
                push!(rows,Dict("id"=>r["id"]*":volume-reseed:"*p["phase"],"native"=>phase_values(model,T,inv(v),x,p["phase"])))
            end
        end
    end
    return rows
end
function timing_row(q,tier,kind,rep,sample,elapsed,work)
    Dict("id"=>q["id"],"tier"=>tier,"phase_work"=>work,"selection"=>get(q,"diameter",q["case"]),
        "warm"=>kind,"rep"=>rep,"sample"=>sample,"elapsed_ns"=>elapsed,
        "julia_compile_enabled"=>Base.JLOptions().compile_enabled,"julia_opt_level"=>Base.JLOptions().opt_level)
end
Base.@noinline function main()
    diag = !isempty(ARGS) && ARGS[1] == "--diagnose"
    args = diag ? ARGS[2:end] : ARGS
    input = JSON.parsefile(args[1])
    rep = length(args)>2 ? parse(Int,args[3]) : 0
    records,timings,absent = [],[],["Clapeyron solver iteration/residual and association mass-action residual are not exposed by these returned values", "D6: Clapeyron gradients excluded; no gradient availability claim"]
    models = Dict(); warm_elapsed_ns = UInt64(0)
    for q in input["requests"]
        if q["operation"] == "gradient" || get(q,"reference_m",0.)>0
            push!(absent,q["id"]*": diagnostic finite FeOs reference or adopted D6 exclusion")
            continue
        end
        rep>0 && !(get(q,"timing","") in ("T1","T2")) && continue
        diag && q["case"] != "C3" && continue
        println(stderr,"Clapeyron request ",q["id"])
        try
            key = (q["case"],get(q,"diameter",""),get(q,"electrostatic",""),get(q,"kcl",0.),get(q,"dipole_convention",""))
            if !haskey(models,key)
                start = time_ns()
                models[key] = Base.invokelatest(build_model,input["models"][q["case"]],q,input["constants"])
                push!(timings,timing_row(q,"model","cold",rep,0,time_ns()-start,"model construction"))
            end
            model,factors = models[key]
            start = time_ns()
            d = Base.invokelatest(call_block,model,q)
            push!(timings,timing_row(q,get(q,"timing",""),"cold",rep,0,time_ns()-start,"first case-output call"))
            r = Base.invokelatest(qualify,model,q,d)
            r["association_parameter_check"] = factors
            if diag
                mw = input["models"]["C3"]["components"][1]["molar_mass"]
                append!(records,Base.invokelatest(diagnose_request,model,q,r,mw))
            end
            push!(records,r)
            tier = get(q,"timing","")
            if rep>0 && tier in ("T1","T2") && !diag
                warm_start = time_ns()
                for _ in 1:5
                    Base.invokelatest(call_block,model,q)
                end
                for sample in 1:25
                    start = time_ns()
                    Base.invokelatest(call_block,model,q)
                    elapsed=time_ns()-start
                    elapsed<=180e9 || error("Warm case exceeded 180 seconds")
                    push!(timings,timing_row(q,tier,"warm",rep,sample,elapsed,
                        q["operation"] == "activity" ? "both helpers, repeated solution/reference solves" : q["operation"]))
                end
                warm_elapsed_ns += time_ns()-warm_start
            end
        catch e
            push!(records,Dict("id"=>q["id"],"case"=>q["case"],"error"=>sprint(showerror,e,catch_backtrace())))
            println(stderr,q["id"],": ",sprint(showerror,e))
        end
        open(io->JSON.print(io,json_numbers(Dict("records"=>records,"timings"=>timings,"absent_diagnostics"=>absent))),args[2]*".tmp","w")
    end
    push!(timings,timing_row(Dict("id"=>"Julia-import","case"=>"all"),"import","cold",rep,0,IMPORT_NS,"JSON, LinearAlgebra, Clapeyron, ForwardDiff imports"))
    tmp = args[2]*".tmp"
    open(io->JSON.print(io,json_numbers(Dict("records"=>records,"timings"=>timings,"warm_elapsed_ns"=>warm_elapsed_ns,"absent_diagnostics"=>absent))),tmp,"w")
    mv(tmp,args[2];force=true)
end
if !isempty(ARGS) && ARGS[1] == "--association-probe"
    rows=association_probe(JSON.parsefile(ARGS[2]),JSON.parsefile(ARGS[3]))
    open(io->JSON.print(io,json_numbers(rows)),ARGS[4],"w")
elseif !isempty(ARGS) && ARGS[1] == "--ladder"
    input=JSON.parsefile(ARGS[2]); rows=[]
    for (lib,file) in (("FeOs",ARGS[3]),("Clapeyron",ARGS[4]))
        native=Dict(r["id"]=>r for r in JSON.parsefile(file)["records"])
        for diameter in (lib=="FeOs" ? ("0.88sigma",) : ("sigma","0.88sigma"))
            solution=native["C3-$diameter-1"]["phases"][1]
            zero=native["C3-ref-$diameter-0"]["phases"][1]["lnphi"]
            payload=Dict("lnphi"=>solution["lnphi"],"pure_lnphi"=>zero,"reference_lnphi"=>zero,
                "m"=>1.,"MW"=>input["models"]["C3"]["components"][1]["molar_mass"],"xwater"=>solution["x"][1])
            g0=contract(payload)["gamma"]; fit_m=[];fit_y=[]
            for m in (1e-4,1e-6,1e-8,1e-10,1e-12)
                reference=only(r for (id,r) in native if startswith(id,"C3-ref-$diameter-") && get(r,"phases",[])[1]["x"][2]>0 && isapprox(r["phases"][1]["x"][2]/r["phases"][1]["x"][1]/payload["MW"],m;rtol=1e-12))
                payload["reference_lnphi"]=reference["phases"][1]["lnphi"]
                shift=log(contract(payload)["gamma"])-log(g0)
                push!(rows,Dict("id"=>"$lib:$diameter:$m","molality"=>m,"finite_minus_zero_ln_gamma"=>shift))
                if 1e-10<=m<=1e-6;push!(fit_m,sqrt(m));push!(fit_y,shift);end
            end
            design=hcat(ones(length(fit_m)),Float64.(fit_m));coeff=design\Float64.(fit_y)
            push!(rows,Dict("id"=>"$lib:$diameter:fit","intercept"=>coeff[1],"sqrt_m_slope"=>coeff[2],"max_fit_residual"=>maximum(abs.(design*coeff-Float64.(fit_y)))))
        end
    end
    open(io->JSON.print(io,rows),ARGS[5],"w")
elseif !isempty(ARGS) && ARGS[1] == "--screening"
    for line in eachline(stdin)
        q = JSON.parse(line)
        c = q["constants"]
        for bits in (256,512)
            setprecision(BigFloat,bits) do
                lb = BigFloat(c["e"])^2/(4big(pi)*BigFloat(c["epsilon0"])*BigFloat(q["epsilon_r"])*BigFloat(c["kB"])*BigFloat(q["T"]))
                k = sqrt(4big(pi)*lb*BigFloat(c["NA"])*BigFloat(q["rho"])*sum(BigFloat.(q["x"])[2:3]))
                ys = k .* BigFloat.(q["diameters"])
                d = Dict("id"=>q["id"],"precision_bits"=>bits,"bjerrum_length"=>string(lb),
                    "kappa_d"=>string.(ys),"source_function"=>source_functions(ys),"producer"=>"independent source expression at supplied native state")
                d["feos_float64_source"] = [iszero(y) ? nothing : Dict("f"=>f_feos(Float64(y)),
                    "f_prime"=>ForwardDiff.derivative(f_feos,Float64(y)),
                    "f_second"=>ForwardDiff.derivative(t->ForwardDiff.derivative(f_feos,t),Float64(y))) for y in ys]
                println(JSON.json(d));flush(stdout)
            end
        end
    end
elseif abspath(PROGRAM_FILE) == (@__FILE__)
    main()
end
