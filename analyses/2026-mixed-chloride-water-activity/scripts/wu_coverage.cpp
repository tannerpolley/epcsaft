#include "surface/parameter_builder.hpp"
#include "properties/helmholtz_model.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
using namespace epcsaft;
using Row = std::map<std::string,std::string>;
using Rows = std::vector<Row>;
const std::filesystem::path leaf="analyses/2026-mixed-chloride-water-activity", input=leaf/"inputs/issue-902-p12-r1", output=leaf/"results/issue-902-p12-r1";
std::vector<std::string> split(const std::string& s) {
    std::vector<std::string> v; std::string word; bool quoted=false;
    for(std::size_t i=0;i<s.size();++i) { char c=s[i]; if(c=='"') { if(quoted && i+1<s.size() && s[i+1]=='"') { word+=c; ++i; } else quoted=!quoted; } else if(c==',' && !quoted) { v.push_back(word); word.clear(); } else if(c!='\r') word+=c; }
    v.push_back(word); return v;
}
Rows read(const std::string& file) {
    auto path=input/file;
    if(file=="acid-replay.csv") path="analyses/2026-hydrochloric-acid-chloride-basis/results/packet-readback-r4/numerical/training.csv";
    else if(file=="frozen-realizations.csv") path="analyses/2026-lithium-chloride-refit/inputs/frozen-realizations.csv";
    else if(file=="yao1992-licl-mgcl2-table1.csv") path="data/packets/yao-1992-licl-mgcl2-observations/1/observations/"+file;
    else if(!std::filesystem::exists(path)) path="analyses/2026-hydrochloric-acid-chloride-basis/inputs/"+file;
    std::ifstream f(path); if(!f) throw std::runtime_error("missing "+file);
    std::string s; do { std::getline(f,s); } while(s.rfind('#',0)==0);
    auto names=split(s); Rows rows;
    while(std::getline(f,s)) { if(s.empty()) continue; auto values=split(s); if(values.size()!=names.size()) throw std::runtime_error("CSV width "+file); Row row; for(std::size_t i=0;i<names.size();++i) row[names[i]]=values[i]; rows.push_back(row); } return rows;
}
double num(const Row& r,const std::string& k) { auto it=r.find(k); if(it==r.end() || it->second.empty() || it->second=="-") return NAN; try { return std::stod(it->second.front()=='('?it->second.substr(1):it->second); } catch(...) { throw std::runtime_error("nonnumeric "+k+"="+it->second); } }
std::string str(const Row& r,const std::string& k) { auto it=r.find(k); return it==r.end()?"":it->second; }
std::string quote(const std::string& s) { std::string t="\""; for(char c:s) t+=c=='"'?"\"\"":std::string(1,c); return t+'"'; }
std::ofstream table(const std::string& name,const std::string& header) { std::ofstream f(output/name); f<<std::setprecision(17)<<header<<'\n'; return f; }
struct Sample {
    std::string set,id,source,system,line,role,phase="not reported",reference="printed", exclusion;
    std::array<double,6> m{}; double aw=NAN,phi=NAN,gamma=NAN,uaw=NAN,uphi=NAN; int salt=3;
    double I() const { return m[0]+m[1]+m[2]+3*m[3]+3*m[4]+m[5]; }
    double q() const { return 2*(m[0]+m[1]+m[2]+m[5])+3*(m[3]+m[4]); }
};
bool same_molality(double a,double b) {
    return std::abs(a-b)<=1e-9*std::max(std::abs(a),std::abs(b));
}
int band(double I) { return I<=3?0:I<=9.9?1:I<=12.8?2:I<=15?3:4; }
int main() try {
    ParameterDocument d;
    for(const auto& r:read("parameters.csv")) {
        const auto f=str(r,"field"),v=str(r,"value");
        if(f=="component_ids") d.component_ids.push_back(v);
        else if(f=="site_component") d.site_component.push_back(std::stoi(v));
        else if(f=="site_multiplicity") d.site_multiplicity.push_back(std::stoi(v));
        else if(f=="base_choice") d.base_choice=v; else if(f=="association_choice") d.association_choice=v;
        else if(f=="electrolyte_choice") d.electrolyte_choice=v; else if(f=="ionic_dispersion_choice") d.ionic_dispersion_choice=v;
        else if(f=="permittivity_choice") d.permittivity_choice=v; else if(f=="polar_choice") d.polar_choice=v;
        else { ParameterEntry e{str(r,"identity"),f,{},num(r,"value"),str(r,"unit")}; std::stringstream s(str(r,"indices")); std::size_t i; while(s>>i) e.index.push_back(i); d.entries.push_back(e); }
    }
    auto p=build_parameters(d); HelmholtzModel model(p); const auto n=p.component_count();
    const auto index=[&](const std::string& id) { auto it=std::find(d.component_ids.begin(),d.component_ids.end(),id); if(it==d.component_ids.end()) throw std::runtime_error(id); return std::size_t(it-d.component_ids.begin()); };
    const auto water=index("water"),cl=index("chloride-anion");
    const std::array<std::size_t,6> ions={index("lithium-cation"),index("sodium-cation"),index("potassium-cation"),index("magnesium-cation"),index("calcium-cation"),index("hydronium-cation")};
    const double mw=p.molar_mass_kg_per_mol[water];
    if(p.k_ij[ions[5]*n+water]!=-2.4255203694263883 || p.k_ij[ions[5]*n+cl]!=-2.25) throw std::runtime_error("H endpoint changed");
    auto coordinates=table("parameter-readback.csv","identity,field,indices,value,unit");
    for(const auto& e:d.entries) { coordinates<<quote(e.identity)<<','<<e.field<<','; for(auto i:e.index) coordinates<<i<<' '; coordinates<<','<<e.magnitude<<','<<e.unit<<'\n'; }
    std::vector<Sample> sources; const auto yao_roles=read("yao-roles.csv");
    const auto load=[&](const std::string& file) {
        std::size_t row=0;
        for(const auto& r:read(file)) {
            Sample s; s.set=str(r,"set"); s.source=file; s.id=file+":"+std::to_string(++row); s.phase=str(r,"phase"); s.line=str(r,"y_I_printed"); s.role="external measured comparison; already inspected";
            const auto salt=str(r,"salt"); s.salt=salt=="CaCl2"?4:3; s.system=s.salt==4?"Ca":"Mg";
            if(file=="el-guendouzi2001-mg-ca.csv" || file=="azougen2010-binary-Ca.csv") { s.m[s.salt]=num(r,"m"); s.aw=num(r,"aw_obs"); s.phi=num(r,"phi_obs"); }
            else if(file=="azougen2010-mixtures.csv") { s.m[0]=num(r,"m_LiCl"); s.m[s.salt]=num(r,"m_divalent"); s.system="Li-"+s.system; s.aw=num(r,"aw_obs"); s.phi=num(r,"phi_obs"); }
            else if(file=="el-guendouzi2005-table1.csv") {
                s.m[3]=num(r,"m_MgCl2"); s.m[str(r,"system")=="MgCl2-NaCl"?1:4]=num(r,"m_other"); s.system=str(r,"system"); s.aw=num(r,"aw_obs"); s.phi=num(r,"phi_obs");
                if(s.system=="MgCl2-NaCl") s.exclusion="report-only printed-property conflicts";
                const int page=int(num(r,"printed_page")),j=int(num(r,"row"));
                if(s.system=="MgCl2-CaCl2" && ((s.line=="0.33" && page==116 && j==11)||(s.line=="0.67" && page==117 && (j==8||j==9||j==11)))) s.exclusion="contradictory printed composition";
            } else if(file=="goldberg-high-m.csv" || file=="rard1981-tableIX.csv") {
                s.m[s.salt]=num(r,"m"); s.aw=num(r,file=="goldberg-high-m.csv"?"aw_evaluated":"aw_calc"); s.phi=num(r,file=="goldberg-high-m.csv"?"phi_evaluated":"phi_calc"); s.gamma=num(r,file=="goldberg-high-m.csv"?"gamma_evaluated":"gamma_calc"); s.role="dependent evaluated comparison";
            } else if(file.rfind("rard-",0)==0) {
                s.m[3]=num(r,"m_MgCl2_token");
                s.phi=num(r,"phi");
                s.aw=num(r,"aw");
                s.uphi=num(r,"propagated_reference_phi_scale");
                s.uaw=num(r,"propagated_reference_aw_linear_scale");
                s.reference=str(r,"reference_equation")+" "+str(r,"reference_column");
                s.line=str(r,"reference_column");
                s.role=file=="rard-reconstructed.csv"?"external source-referenced primary":"correlated reference alternative";
            } else if(file.rfind("robinson-",0)==0) {
                s.m[4]=str(r,"mCa_token")=="-"?0:num(r,"mCa_token"); s.m[3]=str(r,"mMg_token")=="-"?0:num(r,"mMg_token"); if(file=="robinson-reference-reconstructed.csv") { s.m[4]=num(r,"mCa_reference_token"); s.m[3]=0; }
                s.system=s.m[4]>0 && s.m[3]>0?"MgCl2-CaCl2":s.m[3]>0?"Mg":"Ca"; s.line=str(r,"series"); s.role="comparison only; absolute_target_ready=False";
                for(const auto& ref:{"anchor","goldberg1b","goldberg2b"}) { auto a=s;
                a.reference=ref;
                a.id+="/"+std::string(ref);
                a.phi=num(r,std::string("phi_")+ref+(std::string(ref)=="anchor"?"_estimate":"_diagnostic"));
                a.aw=num(r,std::string("aw_")+ref+(std::string(ref)=="anchor"?"_estimate":"_diagnostic"));
                if(!std::isfinite(a.aw)) a.aw=std::exp(-.0180153*a.q()*a.phi);
                sources.push_back(a);
                } continue;
            } else if(file.rfind("yao",0)==0) {
                s.m[0]=num(r,"m_LiCl"); s.m[3]=num(r,"m_MgCl2"); s.aw=num(r,"a_w_set"); s.phi=num(r,"phi"); s.system="Li-Mg"; s.line=str(r,"set"); s.role="calibration; high-I same-source extension"; if(str(r,"flag")=="illegible") s.phi=NAN;
            } else if(file=="hamer-wu-hcl.csv" || file=="messnaoui2020-hcl.csv") {
                s.salt=5; s.system="HCl"; s.m[5]=num(r,"molality"); s.phi=num(r,file.rfind("hamer",0)==0?"phi":"phi_adjusted"); s.aw=num(r,"aw_adjusted"); s.gamma=num(r,"gamma"); s.role=file.rfind("hamer",0)==0?"calibration curve; extension above 7m":"calibration-selected adjusted-water guard";
                if(file.rfind("messnaoui",0)==0) { auto a=s; a.id+="/unadjusted"; a.reference="unadjusted sensitivity"; a.aw=num(r,"aw_unadjusted"); a.phi=num(r,"phi_unadjusted"); sources.push_back(a); }
            } else {
                s.salt=5; s.system=file.rfind("jiang",0)==0?"H-Li":file=="khoo1977-hcl-ca.csv"?"H-Ca":"H-Mg";
                const double I=num(r,"ionic_strength"),y=num(r,"y_background"); const int bg=s.system=="H-Li"?0:s.system=="H-Ca"?4:3;
                s.m[5]=bg==0?num(r,"m_HCl"):I*(1-y); s.m[bg]=bg==0?num(r,"m_LiCl"):I*y/3;
                s.gamma=std::pow(10.,num(r,bg==0?"log10_gamma":"negative_log10_gamma")*(bg==0?1:-1)); s.line=std::to_string(I); if(bg==0 && I==5) s.role="calibration";
            }
            if(file.rfind("yao",0)==0) { s.role += "; "+str(yao_roles.at(row-1),"role"); if(str(yao_roles.at(row-1),"role")=="excluded") s.exclusion="inherited illegible phi exclusion"; }
            if(std::count_if(s.m.begin(),s.m.end(),[](double m){return m>0;})==1) s.role+="; pure endpoint";
            if(s.line.empty()) s.line="binary";
            sources.push_back(s);
        }
    };
    for(const auto& f:{"el-guendouzi2001-mg-ca.csv","el-guendouzi2005-table1.csv","azougen2010-mixtures.csv","azougen2010-binary-Ca.csv","goldberg-high-m.csv","rard1981-tableIX.csv","rard-reconstructed.csv","rard-reference-paths.csv","robinson-reconstructed.csv","robinson-reference-reconstructed.csv","yao1992-licl-mgcl2-table1.csv","hamer-wu-hcl.csv","messnaoui2020-hcl.csv","jiang1996-hcl-li.csv","roy1980-hcl-mg.csv","khoo1977-hcl-mg.csv","khoo1977-hcl-ca.csv"}) load(f);
    const auto nominal_sources=sources;
    for(const auto& scenario:read("frozen-realizations.csv")) if(str(scenario,"source")=="Yao") {
        const auto target=str(scenario,"row"),field=str(scenario,"field"); const double sign=num(scenario,"sign");
        for(const auto& original:nominal_sources) if(original.source.rfind("yao",0)==0 && original.exclusion.empty() && (target=="all" || original.id==original.source+":"+target)) {
            auto a=original; a.id+="/"+str(scenario,"id"); a.reference="printed rounding "+str(scenario,"id");
            if((field=="m_LiCl"||field=="all") && a.m[0]>0) a.m[0]+=sign*.00005;
            if((field=="m_MgCl2"||field=="all") && a.m[3]>0) a.m[3]+=sign*.00005;
            if((field=="phi"||field=="all") && std::isfinite(a.phi)) a.phi+=sign*.00005;
            sources.push_back(a);
        }
    }
    auto inventory=table("source-inventory.csv","id,source,system,line,role,phase,reference,exclusion,I,q,mLi,mNa,mK,mMg,mCa,mH,aw,phi,gamma,evaluate");
    std::map<std::string,double> upper;
    for(const auto& s:sources) if(s.I()>17) { auto key=s.source+s.system+s.line+s.reference; if(upper.find(key)==upper.end()||s.I()<upper[key]) upper[key]=s.I(); }
    std::vector<bool> selected;
    for(const auto& s:sources) {
        bool use=s.exclusion.empty() && (s.salt==5 || s.I()<=17 || s.I()==upper[s.source+s.system+s.line+s.reference]); selected.push_back(use);
        inventory<<quote(s.id)<<','<<s.source<<','<<s.system<<','<<s.line<<','<<quote(s.role)<<','<<quote(s.phase)<<','<<quote(s.reference)<<','<<quote(s.exclusion)<<','<<s.I()<<','<<s.q(); for(auto m:s.m) inventory<<','<<m; inventory<<','<<s.aw<<','<<s.phi<<','<<s.gamma<<','<<use<<'\n';
    }
    auto pred=table("predictions.csv","id,pressure,b,I,q,mLi,mNa,mK,mMg,mCa,mH,mCl,fMg,yLiI,yNaI,yKI,yMgI,yCaI,yHI,aw,phi,lnGammaLi,lnGammaNa,lnGammaK,lnGammaMg,lnGammaCa,lnGammaH,lnALi,lnANa,lnAK,lnAMg,lnACa,lnAH");
    auto errors=table("errors.csv","id,pressure,observable,prediction,reference,signed_error,relative_percent,source_scale,availability,dependence");
    auto diag=table("numerical.csv","id,pressure,owner,density,slope,residual,stop_limit,identity_error,cause,passes");
    std::map<std::string,std::array<double,6>> values;
    struct Error { std::size_t row; std::string obs; double e,u; }; std::vector<Error> residuals;
    bool complete=true;
    const auto close=[&](const State& state,const std::string& id,double pressure,const std::string& owner,double identity,const std::string& cause) {
        const auto c=state.density_convergence; const double stop=c && c->criterion?c->criterion->value:NAN;
        const bool ok=c && c->residual<=stop && state.molar_density>0 && state.pressure_density_derivative>0 && std::isfinite(identity) && std::abs(identity)<=1e-10;
        diag<<quote(id)<<','<<pressure<<','<<owner<<','<<state.molar_density<<','<<state.pressure_density_derivative<<','<<(c?c->residual:NAN)<<','<<stop<<','<<identity<<','<<quote(cause)<<','<<ok<<'\n'; complete=complete&&ok;
    };
    std::vector<Sample> wu; std::vector<double> bases;
    for(const auto& name:{"feed","scrub","regeneration","strip"}) for(double b:{.8,1.,1.2}) {
        Sample s; s.id=std::string(name)+"/"+std::to_string(b); s.system=name;
        if(s.system=="feed") s.m={1.36/6.94,5.75/22.98976928,6.24/39.0983,52.1/24.305,84.4/40.078,0};
        if(s.system=="scrub") s.m={.01,0,0,4.5,0,0}; if(s.system=="regeneration") s.m[3]=4; if(s.system=="strip") { s.m[5]=6; s.salt=5; }
        for(auto& m:s.m) { m/=b; }
        wu.push_back(s); bases.push_back(b);
    }
    auto checks=table("checks.csv","check,id,value,reference,error,tolerance,passes");
    const auto check=[&](const std::string& name,const std::string& id,double value,double expected,double tolerance) {
        const bool ok=std::isfinite(value) && std::abs(value-expected)<=tolerance; complete=complete&&ok;
        checks<<name<<','<<quote(id)<<','<<value<<','<<expected<<','<<value-expected<<','<<tolerance<<','<<ok<<'\n';
    };
    for(const auto& s:sources) if(s.source=="rard-reconstructed.csv" && std::abs(s.m[3]-4.9698)<1e-12) {
        check("Rard phi",s.id,s.phi,3.05894495999748,3.1e-10);
        check("Rard water",s.id,std::exp(-.0180153*s.q()*s.phi),.439716798363251,4.4e-11);
    }
    for(std::size_t i=0;i<wu.size();++i) { const auto& s=wu[i]; const double b=bases[i];
        if(s.system=="scrub") { check("scrub Li",s.id,s.m[0]*b,.01,1e-12); check("scrub Mg",s.id,s.m[3]*b,4.5,4.5e-10); }
        if(s.system=="strip") check("strip H",s.id,s.m[5]*b,6,6e-10);
        if(s.system=="feed") for(std::size_t j=0;j<5;++j) check("feed g/L",s.id,s.m[j]*b*std::array<double,5>{6.94,22.98976928,39.0983,24.305,40.078}[j],std::array<double,5>{1.36,5.75,6.24,52.1,84.4}[j],1e-8);
    }
    for(double pressure:{100000.,101325.}) {
        std::vector<double> solvent(n,0); solvent[water]=1; const auto wr=solvent_reference(model,298.15,pressure,solvent);
        std::vector<InfiniteDilutionReference> refs;
        for(std::size_t i=0;i<6;++i) refs.push_back(infinite_dilution_reference(model,298.15,pressure,solvent,{ions[i],cl},{1,p.charge_number[ions[i]]}));
        if(!wr.pure_solvent_state) throw std::runtime_error("missing solvent reference: "+wr.solvent_cause);
        if(wr.pure_solvent_state) close(*wr.pure_solvent_state,"water",pressure,"solvent_reference",0,wr.solvent_cause);
        for(std::size_t i=0;i<6;++i) if(refs[i].solvent_state) close(*refs[i].solvent_state,"reference/"+std::to_string(i),pressure,"infinite_dilution",0,refs[i].solvent_cause);
        const auto evaluate=[&](const Sample& s,double b) {
            std::vector<double> m(n,0); for(std::size_t i=0;i<6;++i) { m[ions[i]]=s.m[i]; m[cl]+=s.m[i]*p.charge_number[ions[i]]; }
            double charge=0,I=0; for(std::size_t i=0;i<n;++i) { charge+=m[i]*p.charge_number[i]; I+=.5*m[i]*p.charge_number[i]*p.charge_number[i]; }
            auto x=m; x[water]=1/mw; const double total=std::accumulate(x.begin(),x.end(),0.); for(auto& v:x) v/=total;
            if(std::abs(charge)>1e-12+1e-10*s.q() || std::abs(I-s.I())>1e-12+1e-10*I || std::abs(std::accumulate(x.begin(),x.end(),0.)-1)>1e-12) throw std::runtime_error("composition identity");
            const auto a=solvent_activity(model,298.15,pressure,x,wr); const double aw=a.solvent_activity.value_or(NAN),phi=aw>0?-std::log(aw)/(mw*s.q()):NAN;
            std::array<double,6> logs; logs.fill(NAN);
            for(std::size_t i=0;i<6;++i) if(s.m[i]>0) {
                const auto g=electrolyte_activity(model,298.15,pressure,solvent,{ions[i],cl},{1,p.charge_number[ions[i]]},s.m[i],refs[i],std::nullopt,m); logs[i]=g.mean_ionic_log_molality_scale;
                double direct=NAN; if(g.state.log_fugacity_coefficient) { const auto& f=*g.state.log_fugacity_coefficient; const double z=p.charge_number[ions[i]]; direct=(require_component_property(f[ions[i]])-refs[i].log_fugacity[0]+z*(require_component_property(f[cl])-refs[i].log_fugacity[1]))/(1+z)-std::log(1+mw*s.q()); }
                close(g.state,s.id,pressure,"activity/"+std::to_string(i),logs[i]-direct,"");
            }
            close(a.state,s.id,pressure,"water",0,a.cause); values[s.id+"/"+std::to_string(int(pressure))]={aw,phi,logs[s.salt],0,0,0};
            pred<<quote(s.id)<<','<<pressure<<','<<b<<','<<I<<','<<s.q(); for(auto v:s.m) pred<<','<<v; pred<<','<<m[cl]<<','<<s.m[3]/(s.m[3]+s.m[4]);
            for(std::size_t i=0;i<6;++i) { pred<<','<<s.m[i]*(i==3||i==4?3:1)/I; }
            pred<<','<<aw<<','<<phi;
            for(auto v:logs) pred<<','<<v;
            for(std::size_t i=0;i<6;++i) { const double z=p.charge_number[ions[i]]; pred<<','<<(s.m[i]>0?(1+z)*logs[i]+std::log(s.m[i])+z*std::log(m[cl]):NAN); } pred<<'\n';
            return std::array<double,3>{aw,phi,logs[s.salt]};
        };
        for(std::size_t i=0;i<wu.size();++i) evaluate(wu[i],bases[i]);
        for(std::size_t i=0;i<sources.size();++i) if(selected[i] && (pressure==100000 || sources[i].salt==5)) {
            const auto& s=sources[i]; const auto v=evaluate(s,NAN); const std::array<double,3> observed={s.aw,s.phi,s.gamma>0?std::log(s.gamma):NAN};
            for(int j=0;j<4;++j) {
                const int k=j==3?2:j; double ref=observed[k]; if(k==0 && !std::isfinite(ref) && std::isfinite(s.phi)) ref=std::exp(-.0180153*s.q()*s.phi);
                const double factor=j==3?1+p.charge_number[ions[s.salt]]:1, e=factor*(v[k]-ref),u=k==0?s.uaw:k==1?s.uphi:s.source=="roy1980-hcl-mg.csv"?factor*.0006*std::log(10.):NAN;
                double product_offset=0; if(j==3) { double mcl=0; for(std::size_t c=0;c<6;++c) mcl+=s.m[c]*p.charge_number[ions[c]]; product_offset=std::log(s.m[s.salt])+p.charge_number[ions[s.salt]]*std::log(mcl); } const std::string obs=j==0?"aw":j==1?"phi":j==2?"lnGamma":"lnA";
                errors<<quote(s.id)<<','<<pressure<<','<<obs<<','<<factor*v[k]+product_offset<<','<<factor*ref+product_offset<<','<<e<<','<<(k==2?100*std::expm1(e):NAN)<<','<<u<<','<<(!std::isfinite(ref)?"no source observable":!std::isfinite(v[k])?"numerical failure":"available")<<','<<(k==0 && !std::isfinite(s.aw)?"derived from phi":"source representation; aw and phi dependent")<<'\n';
                if(pressure==100000) residuals.push_back({i,obs,e,u});
            }
        }
    }
    for(const auto& row:read("acid-replay.csv")) {
        const double m=num(row,"molality"); const Sample* sample=nullptr;
        for(const auto& s:sources) if(s.source=="hamer-wu-hcl.csv" && s.m[5]==m) sample=&s;
        if(!sample) throw std::runtime_error("missing acid replay row");
        const auto v=values.at(sample->id+"/100000");
        check("acid replay phi",sample->id,v[1],num(row,"phi"),std::max(1e-10*std::abs(num(row,"phi")),1e-12/(mw*sample->q())));
        check("acid replay gamma",sample->id,std::exp(v[2]),num(row,"gamma"),std::max(1e-12,1e-10*num(row,"gamma")));
    }
    auto guards=table("acid-guards.csv","quantity,bucket,count,mean,worst,mean_limit,worst_limit,passes");
    for(int quantity=0;quantity<3;++quantity) for(bool high:{false,true}) {
        double sum=0,worst=0; int count=0;
        for(const auto& s:sources) if((quantity<2?s.source=="hamer-wu-hcl.csv" && s.m[5]<=7:s.source=="messnaoui2020-hcl.csv" && s.reference=="printed") && (s.m[5]>3)==high) {
            const auto v=values.at(s.id+"/100000"); const double error=quantity==0?v[1]-s.phi:quantity==1?100*std::expm1(v[2]-std::log(s.gamma)):v[0]-s.aw;
            sum+=quantity==1?std::abs(error):error*error; worst=std::max(worst,std::abs(error)); ++count;
        }
        const double mean=quantity==1?sum/count:std::sqrt(sum/count),ml=quantity==0?.011:quantity==1?2:.003,wl=quantity==0?.03:quantity==1?5:.006;
        const bool ok=count==(quantity<2?(high?7:23):(high?6:8)) && mean<=ml && worst<=wl; complete=complete&&ok;
        guards<<(quantity==0?"phi":quantity==1?"gamma":"aw")<<','<<(high?"high":"low")<<','<<count<<','<<mean<<','<<worst<<','<<ml<<','<<wl<<','<<ok<<'\n';
    }
    auto equality=table("robinson-isopiestic.csv","mixture,reference,reference_form,predicted_aw_difference,observed_difference,meaning");
    for(const auto& s:sources) if(s.source=="robinson-reconstructed.csv" && s.system=="MgCl2-CaCl2") for(const auto& ref:sources) if(ref.source=="robinson-reference-reconstructed.csv" && ref.set==s.set && ref.line==s.line && ref.reference==s.reference) {
        equality<<quote(s.id)<<','<<quote(ref.id)<<','<<ref.reference<<','<<values.at(s.id+"/100000")[0]-values.at(ref.id+"/100000")[0]<<",0,relational comparison; not independent absolute datum\n";
    }
    auto pressure_table=table("pressure-sensitivity.csv","id,observable,value_100000,value_101325,difference");
    for(const auto& [key,v]:values) if(key.size()>7 && key.substr(key.size()-7)=="/100000") {
        const auto id=key.substr(0,key.size()-7); auto other=values.find(id+"/101325"); if(other==values.end()) continue;
        for(int j=0;j<3;++j) pressure_table<<quote(id)<<','<<(j==0?"aw":j==1?"phi":"lnGamma")<<','<<v[j]<<','<<other->second[j]<<','<<other->second[j]-v[j]<<'\n';
    }
    struct Stats { int n=0,missing=0; double sum=0,sq=0,max=0,stated=0; std::string worst; bool allu=true; };
    std::map<std::string,Stats> stats;
    for(const auto& e:residuals) { const auto& s=sources[e.row];
    const auto key=s.source+","+s.system+","+s.line+","+quote(s.role)+","+quote(s.phase)+","+quote(s.reference)+","+std::to_string(band(s.I()))+","+e.obs;
    auto& a=stats[key];
    if(!std::isfinite(e.e)) { ++a.missing;
    continue;
    } ++a.n;
    a.sum+=e.e;
    a.sq+=e.e*e.e;
    if(std::abs(e.e)>=a.max) { a.max=std::abs(e.e);
    a.worst=s.id;
    } a.allu=a.allu&&std::isfinite(e.u);
    if(std::isfinite(e.u)) a.stated=std::max(a.stated,std::abs(e.e)+e.u);
    }
    auto summary=table("source-statistics.csv","source,system,line,role,phase,reference,I_band,observable,count,unavailable,signed_mean,rms,B_obs,worst_row,B_stated,uncertainty");
    for(const auto& [key,a]:stats) summary<<key<<','<<a.n<<','<<a.missing<<','<<(a.n?a.sum/a.n:NAN)<<','<<(a.n?std::sqrt(a.sq/a.n):NAN)<<','<<(a.n?a.max:NAN)<<','<<quote(a.worst)<<','<<(a.allu&&a.n?a.stated:NAN)<<','<<(a.allu?"stated reference scale only":"total source uncertainty unavailable")<<'\n';
    auto neighbours=table("nearest.csv","wu,source_row,system,line,reference,delta_I,delta_yLiI,delta_yMgI,delta_mLi,delta_mMg,delta_mCa,delta_mNa,delta_mK,delta_mH,missing_species,side");
    auto attached=table("coverage.csv","wu,observable,source,system,reference,coverage,attachment,count,unavailable,B_obs,B_stated,cause,role,contributing_rows,worst_row,signed_error,phase");
    auto state_summary=table("state-summary.csv","state,b,I,observable,prediction,signed_error,B_obs,coverage,evidence,meaning");
    for(const auto& w:wu) {
        std::map<std::string,double> union_bands;
        std::map<std::string,std::string> contributors,matched_errors,union_worst;
        std::map<std::string,std::vector<std::size_t>> groups;
        for(std::size_t i=0;i<sources.size();++i) if(selected[i] && sources[i].reference.rfind("printed rounding",0)!=0) { const auto& s=sources[i]; if((w.system=="strip") != (s.salt==5)) continue; groups[s.source+","+s.system+","+quote(s.reference)].push_back(i); }
        for(const auto& [group,rows]:groups) {
            std::map<std::string,std::vector<std::size_t>> lines; for(auto i:rows) lines[sources[i].line].push_back(i);
            const auto fraction=[&](const Sample& s) { return s.system=="MgCl2-CaCl2"?s.m[3]/(s.m[3]+s.m[4]):s.system=="Li-Mg"||s.system=="Li-Ca"?s.m[0]/s.I():0.; };
            const double target=sources[rows[0]].system=="MgCl2-CaCl2"?w.m[3]/(w.m[3]+w.m[4]):sources[rows[0]].system=="Li-Mg"||sources[rows[0]].system=="Li-Ca"?w.m[0]/w.I():0.;
            double lower=-INFINITY,upper_fraction=INFINITY;
            for(auto i:rows) { double f=fraction(sources[i]); if(f<=target) lower=std::max(lower,f); if(f>=target) upper_fraction=std::min(upper_fraction,f); }
            std::vector<std::size_t> enclosing;
            for(auto i:rows) if(fraction(sources[i])==lower || fraction(sources[i])==upper_fraction) enclosing.push_back(i);
            std::vector<std::size_t> envelope_rows; bool bracket=false;
            for(const auto& [line,ids]:lines) {
                if(std::any_of(ids.begin(),ids.end(),[&](auto i){return std::find(enclosing.begin(),enclosing.end(),i)!=enclosing.end();})) envelope_rows.insert(envelope_rows.end(),ids.begin(),ids.end());
                std::size_t lo=sources.size(),hi=lo; for(auto i:ids) { if(sources[i].I()<=w.I() && (lo==sources.size()||sources[i].I()>sources[lo].I())) lo=i; if(sources[i].I()>=w.I() && (hi==sources.size()||sources[i].I()<sources[hi].I())) hi=i; }
                const std::array<std::size_t,2> endpoints={lo,hi};
                if(lo<sources.size() && hi<sources.size() && std::all_of(endpoints.begin(),endpoints.end(),[&](auto i){return std::equal(w.m.begin(),w.m.end(),sources[i].m.begin(),[&](double a,double b){return same_molality(a/w.I(),b/sources[i].I());});})) bracket=true;
                for(auto i:std::vector<std::size_t>{lo,hi==lo?sources.size():hi}) if(i<sources.size()) { const auto& s=sources[i];
                neighbours<<quote(w.id)<<','<<quote(s.id)<<','<<s.system<<','<<s.line<<','<<quote(s.reference)<<','<<s.I()-w.I()<<','<<s.m[0]/s.I()-w.m[0]/w.I()<<','<<3*s.m[3]/s.I()-3*w.m[3]/w.I()<<','<<s.m[0]-w.m[0]<<','<<s.m[3]-w.m[3]<<','<<s.m[4]-w.m[4]<<','<<s.m[1]-w.m[1]<<','<<s.m[2]-w.m[2]<<','<<s.m[5]-w.m[5]<<','<<quote(w.system=="feed"?"simultaneous Li Na K Mg Ca; acidity Fe water unknown":w.system=="scrub"?"trace Li composition":w.system=="regeneration"?"carbonate acid contacted inventory":"outlet residual acid metals")<<','<<(lo<sources.size()&&hi<sources.size()?"I bracket only; composition must also bracket":"one-sided")<<'\n';
                }
            }
            for(const auto& obs:{"aw","phi","lnGamma","lnA"}) {
                Stats a; std::string row_ids,matched;
                for(const auto& e:residuals) if(e.obs==obs && band(sources[e.row].I())>=(w.system=="strip"?1:2) && band(sources[e.row].I())<=(w.system=="strip"?1:w.system=="regeneration"?3:4) && std::find(envelope_rows.begin(),envelope_rows.end(),e.row)!=envelope_rows.end()) { if(!std::isfinite(e.e)) { ++a.missing;
                continue;
                } ++a.n;
                if(std::abs(e.e)>=a.max) { a.max=std::abs(e.e); a.worst=sources[e.row].id; }
                row_ids+=sources[e.row].id+" [line "+sources[e.row].line+"]; ";
                if(std::equal(w.m.begin(),w.m.end(),sources[e.row].m.begin(),same_molality)) { std::ostringstream error; error<<std::setprecision(17)<<e.e; matched+=sources[e.row].id+"="+error.str()+"; "; }
                a.allu=a.allu&&std::isfinite(e.u);
                if(std::isfinite(e.u)) a.stated=std::max(a.stated,std::abs(e.e)+e.u);
                }
                const bool water=std::string(obs)=="aw"||std::string(obs)=="phi";
                const auto& first=sources[rows[0]];
                const auto role=first.reference=="unadjusted sensitivity"?"unadjusted correction sensitivity":first.source=="rard1981-tableIX.csv"?"smoothed dependent evaluated comparison":first.role;
                const bool relevant=water?(w.system=="strip"?first.system=="HCl":w.system=="regeneration"?first.system=="Mg":w.system=="scrub"?first.system=="Mg"||first.system=="Li-Mg":first.system=="Mg"||first.system=="Ca"||first.system=="MgCl2-CaCl2"):w.system=="strip"?first.source=="hamer-wu-hcl.csv":first.source=="goldberg-high-m.csv" && first.system=="Mg";
                std::ostringstream scale; scale<<std::setprecision(17)<<a.max;
                if(relevant && a.n && !a.missing) { if(a.max>=union_bands[obs]) { union_bands[obs]=a.max; union_worst[obs]=a.worst; } contributors[obs]+=first.source+"/"+first.system+"/"+first.reference+" ["+role+"; B_obs="+scale.str()+"; maximum at "+a.worst+"]; "; matched_errors[obs]+=matched; }
                const bool binary=w.system=="regeneration" && sources[rows[0]].system=="Mg",acid=w.system=="strip" && sources[rows[0]].system=="HCl";
                const bool available=a.n&&!a.missing,calibration=acid && ((first.source=="hamer-wu-hcl.csv" && w.I()<=7) || (first.source=="messnaoui2020-hcl.csv" && first.reference=="printed" && !matched.empty()));
                const auto label=calibration?"inside training":binary && bracket && water && (first.source=="rard-reconstructed.csv" || first.source=="el-guendouzi2001-mg-ca.csv")?"inside independent evidence; reference dependent":"extrapolation";
                const auto attachment=!available?"NA: no finite evidence-based band":!matched.empty()?"matching complete composition":bracket?"bracketing-evidence summary; between-row error unverified":"borrowed comparison scale; extrapolation error unbounded by available evidence";
                attached<<quote(w.id)<<','<<obs<<','<<group<<','<<quote(label)<<','<<quote(attachment)<<','<<a.n<<','<<a.missing<<','<<(available?a.max:NAN)<<','<<(available&&a.allu?a.stated:NAN)<<','<<quote(a.n?"source-line summary; total source uncertainty unavailable":"NA: no source activity observable")<<','<<quote(role)<<','<<quote(row_ids)<<','<<quote(a.worst)<<','<<quote(matched.empty()?"NA: no exact source match":matched)<<",phase unresolved\n";
            }
        }
        const double b=w.system=="strip"?6/w.m[5]:w.system=="regeneration"?4/w.m[3]:w.system=="scrub"?4.5/w.m[3]:(52.1/24.305)/w.m[3];
        const auto v=values.at(w.id+"/100000");
        double chloride=0; for(std::size_t i=0;i<6;++i) chloride+=w.m[i]*p.charge_number[ions[i]];
        const double product_log=w.system=="strip"?2*v[2]+2*std::log(w.m[5]):3*v[2]+std::log(w.m[3])+2*std::log(chloride);
        for(int j=0;j<4;++j) {
            const std::string obs=j==0?"aw":j==1?"phi":j==2?"lnGamma":"lnA";
            state_summary<<w.system<<','<<b<<','<<w.I()<<','<<(j<2?obs:obs+(w.system=="strip"?"HCl":"MgCl2"))<<','<<(j<2?v[j]:j==3?product_log:v[2])<<','<<quote(matched_errors[obs].empty()?"NA: no exact source match":matched_errors[obs])<<','<<(union_bands.count(obs)?union_bands[obs]:NAN)<<','<<(w.system=="strip"&&w.I()<=7?"inside training":w.system=="regeneration"&&j<2?"inside independent evidence; source dependent":"extrapolation")<<','<<quote("union maximum at "+union_worst[obs]+"; contributors: "+contributors[obs])<<",union comparison scale; Wu error unknown unless matched; not qualification; phase unresolved\n";
        }
        for(const auto& obs:{"Li_product","Na_product","K_product","actual_acidity","contacted_or_outlet_inventory"}) if(w.system!="strip" || std::string(obs)=="contacted_or_outlet_inventory") attached<<quote(w.id)<<','<<obs<<",NA,NA,NA,extrapolation,NA,0,0,nan,nan,NA: no finite evidence-based band,NA,NA,NA,NA,phase unresolved\n";
    }
    auto inventory_wu=table("wu-inventory.csv","state,b,molarity_basis,MgCO3_dose_g_L,dissolved_carbonate,acidity,Fe,water_mass_basis");
    for(std::size_t i=0;i<wu.size();++i) inventory_wu<<wu[i].system<<','<<bases[i]<<','<<quote(wu[i].system=="feed"?"Li Na K Mg Ca g/L = 1.36 5.75 6.24 52.1 84.4":wu[i].system=="scrub"?"MgCl2 4.5 M plus LiCl 0.01 M":wu[i].system=="regeneration"?"MgCl2 4 M":"HCl 6 M")<<','<<(wu[i].system=="regeneration"?15:0)<<",unreported,"<<(wu[i].system=="strip"?"analytical HCl input":"unreported; feed projection H=0 only")<<",unreported; Fe-free projection,assumed sensitivity; kg water/L\n";
    auto spread=table("water-mass-spread.csv","state,observable,min_prediction,max_prediction,meaning");
    for(const auto& state:{"feed","scrub","regeneration","strip"}) for(int j=0;j<2;++j) {
        double lo=INFINITY,hi=-INFINITY; for(const auto& w:wu) if(w.system==state) { const double v=values.at(w.id+"/100000")[j]; lo=std::min(lo,v); hi=std::max(hi,v); }
        spread<<state<<','<<(j==0?"aw":"phi")<<','<<lo<<','<<hi<<",discrete predictions; not measurement uncertainty\n";
    }
    std::ofstream maps(output/"native-mappings.txt"); std::ifstream proc("/proc/self/maps"); std::string line; while(std::getline(proc,line)) if(line.find("libepcsaft_native.so")!=std::string::npos) maps<<line<<'\n';
    return complete?0:1;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
