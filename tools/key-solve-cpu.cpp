#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
constexpr int cases=15, samples=20;
struct Observation {double log_ratio=0,clock_radius=0;};
struct Interval {double center=0,half=0;};
void need(bool b,const char* message){if(!b)throw std::runtime_error(message);}
Interval interval(const std::vector<Observation>& rows) {
    need(rows.size()==samples,"expected 20 balanced CPU observations");
    Interval out;double quantization=0;std::vector<double> sorted;
    for(auto r:rows){sorted.push_back(r.log_ratio);quantization=std::max(quantization,r.clock_radius);}
    std::sort(sorted.begin(),sorted.end());out.center=.5*(sorted[samples/2-1]+sorted[samples/2]);
    double mass=std::ldexp(1.,-samples),tail=0;int k=0;
    for(int i=0;i<samples/2;++i) {
        tail+=mass;if(2*tail<=.05/(3*cases))k=i+1;
        mass*=double(samples-i)/(i+1);
    }
    need(k>0,"insufficient CPU observations for simultaneous interval");
    out.half=std::max(out.center-sorted[k-1],sorted[samples-k]-out.center)+quantization;
    return out;
}
struct Judgment {std::string state;double lower,upper,noise;};
Judgment judge(const Interval& pre,const Interval& ab,const Interval& post) {
    const double noise=std::max(std::abs(pre.center)+pre.half,std::abs(post.center)+post.half);
    Judgment out{"inconclusive",ab.center-ab.half-noise,ab.center+ab.half+noise,noise};
    const bool control=std::abs(pre.center)<=pre.half && std::abs(post.center)<=post.half;
    if(!control)return out;
    if(out.upper<=std::log1p(.01))out.state="pass";
    else if(out.lower>std::log1p(.01))out.state="fail";
    return out;
}
int main(int argc,char**argv) {
    std::cout<<std::setprecision(17);
    if(argc==2 && std::string(argv[1])=="--selfcheck") {
        Interval control{0,.0001};
        need(judge(control,{0,.0001},control).state=="pass","zero overhead did not pass");
        need(judge(control,{std::log1p(.02),.0001},control).state=="fail","2 percent overhead did not fail");
        need(judge(control,{std::log1p(.01),.002},control).state=="inconclusive","uncertainty did not remain inconclusive");
        need(judge({.01,.0001},{0,.0001},control).state=="inconclusive","drifted control accepted");
        std::vector<Observation> rows;for(int i=0;i<samples;++i)rows.push_back({(i%2?.001:-.001),0});
        auto bounded=interval(rows);need(bounded.half>0 && bounded.center==0,"CPU interval calculation failed");
        std::cout<<"zero=pass,2pct=fail,threshold_uncertain=inconclusive,control_drift=inconclusive\n"
            <<"alternating_0.001_log_ratio_half_width="<<bounded.half<<'\n';return 0;
    }
    need(argc==2,"usage: r2-cpu sample-file | --selfcheck");
    std::ifstream in(argv[1]);need(bool(in),"missing samples");
    using Phases=std::array<std::vector<Observation>,3>;
    std::map<std::string,Phases> all;std::map<std::string,std::array<std::array<bool,samples>,3>> seen;
    std::map<std::string,int> batches;
    std::string phase,name;int block,k;double tick,a,b,c,d;
    std::string line;
    while(std::getline(in,line)) {
        if(line.empty())continue;
        std::istringstream row(line);std::string extra;
        need(bool(row>>phase>>name>>block>>k>>tick>>a>>b>>c>>d),"malformed CPU sample row");
        need(!(row>>extra),"extra CPU sample fields");
        const int p=phase=="aa_pre"?0:phase=="ab"?1:phase=="aa_post"?2:-1;
        need(p>=0 && block>=0 && block<samples && k>0 && tick>0 && std::isfinite(tick),"invalid CPU sample metadata");
        need(!seen[name][p][block],"duplicate CPU block");seen[name][p][block]=true;
        if(!batches.count(name))batches[name]=k;need(batches[name]==k,"batch size changed");
        const double delta=2*k*tick;double radius=0;
        for(double t:{a,b,c,d}){need(std::isfinite(t)&&t>delta,"invalid CPU duration");radius+=-.5*std::log1p(-delta/t);}
        all[name][p].push_back({.5*(std::log(b/a)+std::log(d/c)),radius});
    }
    need(!in.bad(),"CPU sample read failed");need(all.size()==cases,"expected every one of 15 CPU cases");
    for(const auto& base:{"case-6/figure-2-point-7/","case-7/flash-5918.5286-Pa/"})
        for(const auto* mode:{"cold","accepted"})need(all.count(std::string(base)+mode),"missing cold/accepted case");
    for(int i=1;i<=8;++i)need(all.count("khudaida/point-"+std::to_string(i)),"missing sequence point");
    for(const auto* name:{"mea-nominally-unloaded/353.15-K","mea-f1/residual-jacobian","mea-cp/hilliard-353.15-K-0.358"})
        need(all.count(name),"missing unloaded/F1/Cp case");
    std::cout<<"case,median_log_ratio,noise_log,lower_percent,upper_percent,judgment\n";
    bool fail=false,inconclusive=false;
    for(const auto& [name,phases]:all) {
        auto pre=interval(phases[0]),ab=interval(phases[1]),post=interval(phases[2]);auto j=judge(pre,ab,post);
        fail|=j.state=="fail";inconclusive|=j.state=="inconclusive";
        std::cout<<name<<','<<ab.center<<','<<j.noise<<','<<100*std::expm1(j.lower)<<','<<100*std::expm1(j.upper)<<','<<j.state<<'\n';
    }
    return fail?1:inconclusive?2:0;
}
