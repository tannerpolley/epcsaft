# FeOs-minus-Engine source-digit diagnostic, Gross2001 Eqs.4-10/full-eos §4.
using JSON, ForwardDiff, DelimitedFiles
inp=JSON.parsefile(ARGS[1]);data=JSON.parsefile(ARGS[2]);table=readdlm(ARGS[3],',',String;skipstart=1)
coeff=Dict((l,a)=>zeros(7) for l in ("Engine","FeOs") for a in ("A0","A1","A2","B0","B1","B2"))
for row in eachrow(table);coeff[(row[1],row[2])][parse(Int,row[3])+1]=parse(Float64,row[4]);end
dip=Dict((l,a)=>zeros(3,5) for l in ("Engine","FeOs") for a in ("A","B","C"))
for row in eachrow(readdlm(joinpath(dirname(ARGS[3]),"dipole-coefficients.csv"),',',String;skipstart=1));dip[(row[1],row[2])][parse(Int,row[3])+1,parse(Int,row[4])+1]=parse(Float64,row[5]);end
requests=Dict(q["id"]=>q for q in inp["requests"]);rows=[]
for rec in data["records"]
 q=requests[rec["id"]];q["operation"]=="fixed" || continue
 cs=inp["models"][q["case"]]["components"];n=length(cs);T=q["T"];m=Float64[c["m"] for c in cs];sig=Float64[c["sigma"] for c in cs];eps=Float64[c["epsilon"] for c in cs]
 if q["case"]=="C3";sig[1]+=10.11exp(-.01775T)-1.417exp(-.01146T);end
 diam=sig.*(1 .- .12exp.(-3eps./T));k=zeros(n,n);z=[get(c,"z",0) for c in cs]
 if n==2;k[1,2]=k[2,1]=inp["models"][q["case"]]["k"];end
 if q["case"]=="C3"
  for (i,j,v) in inp["models"]["C3"]["pairs"];k[Int(i)+1,Int(j)+1]=k[Int(j)+1,Int(i)+1]=v;end
  k[1,3]=k[3,1]=q["kcl"];diam[2:3]=.88sig[2:3]
 end
 function delta(u)
  rho=u[1];x=u[2:end]/sum(u[2:end]);rn=rho*6.02214076e-7;eta=pi/6*rn*sum(x.*m.*diam.^3);mm=sum(x.*m)
  f1=(mm-1)/mm;f2=f1*(mm-2)/mm
  da=coeff[("FeOs","A0")]-coeff[("Engine","A0")]+f1*(coeff[("FeOs","A1")]-coeff[("Engine","A1")])+f2*(coeff[("FeOs","A2")]-coeff[("Engine","A2")])
  db=coeff[("FeOs","B0")]-coeff[("Engine","B0")]+f1*(coeff[("FeOs","B1")]-coeff[("Engine","B1")])+f2*(coeff[("FeOs","B2")]-coeff[("Engine","B2")])
  i1=sum(da[j]*eta^(j-1) for j in 1:7);i2=sum(db[j]*eta^(j-1) for j in 1:7)
  c1=1/(1+mm*(8eta-2eta^2)/(1-eta)^4+(1-mm)*(20eta-27eta^2+12eta^3-2eta^4)/((1-eta)*(2-eta))^2)
  a=zero(rho);b=zero(rho)
  for i in 1:n,j in 1:n
   e=z[i]*z[j]>0 ? 0. : sqrt(eps[i]*eps[j])*(1-k[i,j]);w=x[i]*x[j]*m[i]*m[j]*((sig[i]+sig[j])/2)^3
   a+=w*e/T;b+=w*(e/T)^2
  end
  result=-2pi*rn*i1*a-pi*rn*mm*c1*i2*b
  if n==1 && get(cs[1],"mu",0)>0
   coupling=cs[1]["mu"]^2/(T*m[1]*sig[1]^3)
   chain=min(m[1],2);legs=[1,(chain-1)/chain,(chain-1)*(chain-2)/chain^2]
   function polar(lib)
    si=(1e-21/299792458)^2/(4pi*8.8541878128e-12*1.380649e-23*1e-30)
    reduced=coupling*(lib=="Engine" || get(requests[rec["id"]],"dipole_convention","")=="matched" ? si : 1e-19/1.380649e-23)
    j2=sum(sum(legs.*(dip[(lib,"A")][:,j]+eps[1]/T*dip[(lib,"B")][:,j]))*eta^(j-1) for j in 1:5)
    j3=sum(sum(legs.*dip[(lib,"C")][:,j])*eta^(j-1) for j in 1:5)
    a2=-pi*rn*sig[1]^3*reduced^2*j2;a3=-4pi^2/3*rn^2*sig[1]^6*reduced^3*j3
    return a2/(1-a3/a2)
   end
   result+=polar("FeOs")-polar("Engine")
  end
  return result
 end
 for ph in get(rec,"phases",[])
  u=Float64.(vcat(ph["rho"],ph["x"]));g=ForwardDiff.gradient(delta,u);a=delta(u);dz=u[1]*g[1];dm=a .+ dz .+ g[2:end] .- sum(ph["x"].*g[2:end])
  push!(rows,Dict("id"=>q["id"],"phase"=>ph["phase"],"delta_a"=>a,"delta_Z"=>dz,"delta_mu"=>dm,"delta_lnphi"=>ph["Z"]>0 ? dm .- log((ph["Z"]+dz)/ph["Z"]) : nothing))
 end
end
open(io->JSON.print(io,rows),ARGS[4],"w")
