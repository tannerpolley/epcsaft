#include <Eigen/Dense>
#include "core/linear_error.hpp"
#include <iostream>
using namespace epcsaft::core;
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
constexpr double epsilon = std::numeric_limits<double>::epsilon();
struct Operand {
 double absolute=0, chart=0; bool available=true, resolved=true;
 operator epcsaft::core::NumericalEstimate() const { return {absolute,chart,available,resolved,""}; }
};
#include <quadmath.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>
using Csv = std::map<std::string,std::string>;
std::vector<std::string> split(const std::string& s) { std::vector<std::string> v;std::stringstream in(s);std::string x;while(std::getline(in,x,',')){if(!x.empty()&&x.back()=='\r')x.pop_back();v.push_back(x);}return v; }
std::vector<Csv> read(const std::string& p) {
 std::ifstream in(p); if(!in)throw std::runtime_error(p);std::string line;std::getline(in,line);auto keys=split(line);std::vector<Csv> result;
 while(std::getline(in,line)){auto v=split(line);Csv row;for(std::size_t k=0;k<keys.size();++k)row[keys[k]]=v.at(k);result.push_back(row);}return result;
}
double number(const Csv& r,const std::string& key) { return std::stod(r.at(key)); }
void require(bool b,const std::string& message) {if(!b)throw std::runtime_error(message);}
double ratio(double difference,double radius) {return radius>0?difference/radius:difference==0?0:INFINITY;}
void manufactured(const std::string& out) {
 std::ofstream f(out+"/manufactured.csv");f<<std::setprecision(17)<<"case,row_scale,column_scale,rank,correction,floor_radius,known_error\n";
 const Matrix j=Matrix::Identity(1,1),g=j;const Vector zero=Vector::Zero(1),ones=Vector::Ones(1);
 const Vector eta=canonical_row_floor(j,ones,Vector::Constant(1,2),zero);
 for(int r:{0,-20,20})for(int c:{0,-15,15}) {
  auto e=linear_error(j,zero,g,eta,Vector::Constant(1,std::ldexp(1.,r)),Vector::Constant(1,std::ldexp(1.,c)));
  require(e.radius(0)==2*epsilon,"counterexample scaling");require(e.correction(0)==0,"counterexample correction");
  f<<"y_minus_one,"<<std::ldexp(1.,r)<<','<<std::ldexp(1.,c)<<','<<e.rank<<','<<e.correction(0)<<','<<e.radius(0)<<','<<epsilon<<'\n';
 }
 for(double small:{1.,1e-8,1e-12,0.}) {
  Matrix a=Matrix::Identity(2,2);a(1,1)=small;
  auto e=linear_error(a,Vector::Zero(2),Matrix::Identity(2,2),Vector::Constant(2,1e-15),Vector::Ones(2),Vector::Ones(2));
  const double known=small?1e-15/small:INFINITY;
  if(small)require(e.radius.maxCoeff()>=known&&e.correction.norm()==0,"independent row perturbation");
  else require(e.rank==1,"singularity");
  f<<"diagonal_perturbation,1,1,"<<e.rank<<','<<e.correction.norm()<<','<<(small?e.radius.maxCoeff():INFINITY)<<','<<known<<'\n';
 }
 Matrix a=Matrix::Identity(2,2);a(0,1)=1048576.;
 Vector independent(2);independent<<0,std::ldexp(1.,-50);
 const Vector eta_n=independent.cwiseAbs();auto e=linear_error(a,Vector::Zero(2),Matrix::Identity(2,2),eta_n,Vector::Ones(2),Vector::Ones(2));
 const Vector known=a.fullPivLu().solve(independent);
 require((e.radius.array()>=known.cwiseAbs().array()).all()&&e.correction.norm()==0,"nonnormal floor");
 f<<"nonnormal,1,1,"<<e.rank<<",0,"<<e.radius.maxCoeff()<<','<<known.cwiseAbs().maxCoeff()<<'\n';
 for(int k:{-1,1}) {
  Vector r(2),c(2);r<<std::ldexp(1.,k),std::ldexp(1.,-k);c<<std::ldexp(1.,-k),std::ldexp(1.,k);
  auto scaled=linear_error(a,Vector::Zero(2),Matrix::Identity(2,2),eta_n,r,c);
  require((scaled.radius-e.radius).cwiseAbs().maxCoeff()==0,"nonnormal scaling");
  f<<"nonnormal_scaled,"<<r(0)<<','<<c(0)<<','<<scaled.rank<<",0,"<<scaled.radius.maxCoeff()<<','<<known.cwiseAbs().maxCoeff()<<'\n';
 }
}
void boundaries(const std::string& out) {
 std::ofstream f(out+"/comparison-boundaries.csv"); f<<std::setprecision(17)<<"case,judgment,reason,difference,allowance,chart_difference,radius_sum\n";
 const auto run=[&](const char* name,double a,double b,double allowance,Operand ea,Operand eb,const char* expected,Quantity q=Quantity::Linear,ComparisonEligibility eligible={},ComparisonMode mode=ComparisonMode::EstimatedError) {
  auto r=compare_equilibrium_quantity(q,a,b,allowance,ea,eb,eligible,mode);
  f<<name<<','<<r.judgment<<','<<r.cause<<','<<r.difference<<','<<allowance<<','<<r.chart_difference<<','<<r.radius_sum<<'\n';require(r.judgment==expected,name);
 };
 Operand small{0.75,0.75},large{2,2},missing{2,2,false,true},unresolved{2,2,true,false},nonfinite{INFINITY,INFINITY};
 run("two_small_radii",0,1.5,1,small,small,"failed");
 run("ordinary_boundary",0,1,1,{},missing,"ordinary");
 run("both_radii_equal_allowance",0,1.5,1,{1,1},{1,1},"failed");
 run("one_radius_large",0,2,1,{2,2},{0,0},"estimated_error");
 run("radius_sum_boundary",0,4,1,large,large,"estimated_error");
 run("no_added_allowance",0,4.5,1,large,large,"failed");
 run("no_margin",0,std::nextafter(4.,INFINITY),1,large,large,"failed");
 run("missing_a",0,2,1,missing,large,"failed");
 run("missing_b",0,2,1,large,missing,"failed");
 run("unresolved_a",0,2,1,unresolved,large,"failed");
 run("unresolved_b",0,2,1,large,unresolved,"failed");
 run("nonfinite_estimate",0,2,1,nonfinite,large,"failed");
 run("negative_estimate",0,2,1,{-1,-1},large,"failed");
 run("density_log_pass",1,2,0.01,{1,0.35},{1,0.35},"estimated_error",Quantity::Positive);
 run("packing_log_fail",0.1,0.2,0.001,{1,0.3},{1,0.3},"failed",Quantity::Positive);
 run("log_packing",std::log(0.1),std::log(0.2),0.001,{0.35,0.35},{0.35,0.35},"estimated_error",Quantity::LogPositive);
 run("positive_domain",0,2,0.01,large,large,"failed",Quantity::Positive);
 run("strict_pressure",0,2,1,large,large,"failed",Quantity::Linear,{},ComparisonMode::Strict);
 run("strict_key_fixture",1,1.001,1e-10,large,large,"failed",Quantity::Positive,{},ComparisonMode::Strict);
 run("pressure_own_estimate",100,101,1e-8,{1e-4,1e-4},{1e-4,1e-4},"failed");
 run("nonfinite_operand",NAN,2,1,large,large,"failed");
 run("negative_allowance",0,2,-1,large,large,"failed");
 run("composition_own_estimate",0.4,0.400001,1e-10,{1e-9,1e-9},{1e-9,1e-9},"failed");
 for(int k=0;k<8;++k) {
  ComparisonEligibility eligible;bool* fields[]={&eligible.declaration,&eligible.branch,&eligible.phase_matching,&eligible.inventory,&eligible.residual_a,&eligible.residual_b,&eligible.admission_a,&eligible.admission_b};*fields[k]=false;
  const std::string name="prerequisite_"+std::to_string(k);run(name.c_str(),0,2,1,large,large,"failed",Quantity::Linear,eligible);
 }
 // A compensated dot must retain the unit term lost by an ordinary double sum.
 Matrix cancellation(1,3);cancellation<<1e16,1,-1e16;
 require(compensated_defect(cancellation,Vector::Ones(3),Vector::Zero(1))(0)==-1,"compensated defect");
}
int run(int argc,char** argv) {
 const std::string out=argv[1];manufactured(out);boundaries(out);
 auto rows=read(out+"/estimates.csv");
 using Key=std::tuple<std::string,std::string,std::string>;std::map<Key,std::map<std::pair<std::string,int>,Csv>> endpoints;
 for(auto r:rows)endpoints[{r.at("set"),r.at("row"),r.at("arm")}][{r.at("quantity"),std::stoi(r.at("index"))}]=r;
 std::ofstream domains(out+"/domain-resolution.csv"); domains<<std::setprecision(17)<<"set,row,arm,phase_radius_over_gap,phase_resolved,packing_domain,composition_domain,amount_domain,rank\n";
 std::map<Key,bool> packing_resolved;
 for(const auto& r:read(out+"/packing.csv")) {
  const Key key={r.at("set"),r.at("row"),r.at("arm")};if(!packing_resolved.count(key))packing_resolved[key]=true;
  const double eta=number(r,"packing"),e=number(r,"log_radius");
  packing_resolved[key]=packing_resolved[key]&&eta*std::exp(-e)>=number(r,"lower")&&eta*std::exp(e)<=number(r,"upper");
 }
 std::map<Key,bool> resolved;
 for(const auto& [key,q]:endpoints) {
  const auto& l=q.at({"density",0}); const auto& v=q.at({"density",1});
  const double el=number(l,"estimate_abs"),ev=number(v,"estimate_abs"),gap=number(l,"value")-number(v,"value");
  bool comp=true,amount=true,packing=packing_resolved.at(key);const bool phase=gap>el+ev;
  for(const auto& [name,r]:q) {
   const double x=number(r,"value"),e=number(r,"estimate_abs");
   if(name.first=="composition")comp=comp&&(x-e>=0&&x+e<=1);
   if(name.first=="amount")amount=amount&&(x-e>=0);
  }

  resolved[key]=phase&&comp&&amount&&packing&&number(l,"rank")==number(l,"n");
  domains<<std::get<0>(key)<<','<<std::get<1>(key)<<','<<std::get<2>(key)<<','<<(el+ev)/gap<<','<<phase<<','<<packing<<','<<comp<<','<<amount<<','<<number(l,"rank")<<'\n';
 }
 std::map<Key,std::array<double,3>> criteria;
 for(const auto& r:read(out+"/original-criteria.csv"))criteria[{r.at("set"),r.at("row"),r.at("arm")} ]={number(r,"success"),number(r,"residual"),number(r,"criterion")};
 struct Group {int count=0,uncovered=0,ordinary=0,conditional=0,failed=0;double maximum=0,density_max=0,largest_density_difference=0,largest_pair_radius=0;};
 std::map<std::pair<std::string,double>,Group> groups;
 std::ofstream comparisons(out+"/comparisons.csv");comparisons<<std::setprecision(17)<<"set,T,a,b,quantity,index,difference,chart_difference,radius_a,radius_b,ratio,judgment,reason\n";
 for(auto ia=endpoints.begin();ia!=endpoints.end();++ia)for(auto ib=std::next(ia);ib!=endpoints.end();++ib) {
  const auto& [sa,ra,aa]=ia->first;const auto& [sb,rb,ab]=ib->first;
  if(sa!=sb||(sa!="854"&&ra!=rb))continue;
  const auto& qa=ia->second;const auto& qb=ib->second;
  if(number(qa.begin()->second,"T")!=number(qb.begin()->second,"T"))continue;
  auto& group=groups[{sa,number(qa.begin()->second,"T")}];
  for(const auto& [name,a]:qa) {
   const auto& b=qb.at(name);const double av=number(a,"value"),bv=number(b,"value"),ea=number(a,"estimate_abs"),eb=number(b,"estimate_abs");
   const bool positive=name.first=="density";
   const double eca=positive?number(a,"estimate_chart"):ea,ecb=positive?number(b,"estimate_chart"):eb;
   ComparisonEligibility eligibility;
   const auto ca=criteria.at(ia->first),cb=criteria.at(ib->first);eligibility.admission_a=bool(ca[0]);eligibility.admission_b=bool(cb[0]);eligibility.residual_a=ca[0]&&ca[1]<=ca[2];eligibility.residual_b=cb[0]&&cb[1]<=cb[2];
   const double scale=std::max(std::abs(av),std::abs(bv));
   const double allowance=sa=="854"?1e-10*std::max(positive||name.first=="pressure"?0.:1.,scale):std::max(1e-12,1e-10*scale);
   const bool strict=sa=="854"&&name.first=="pressure";
   const auto result=compare_equilibrium_quantity(positive?Quantity::Positive:Quantity::Linear,av,bv,allowance,{ea,eca,true,resolved.at(ia->first)},{eb,ecb,true,resolved.at(ib->first)},eligibility,strict?ComparisonMode::Strict:ComparisonMode::EstimatedError);
   const double difference=positive?log_difference(av,bv):std::abs(av-bv),r=ratio(difference,eca+ecb);
   ++group.count;group.uncovered+=r>1;group.maximum=std::max(group.maximum,r);
   group.ordinary+=result.judgment=="ordinary";group.conditional+=result.judgment=="estimated_error";group.failed+=result.judgment=="failed";
   if(positive){group.density_max=std::max(group.density_max,r);if(std::abs(av-bv)>group.largest_density_difference){group.largest_density_difference=std::abs(av-bv);group.largest_pair_radius=eca+ecb;}}
   comparisons<<sa<<','<<number(a,"T")<<','<<ra<<'/'<<aa<<','<<rb<<'/'<<ab<<','<<name.first<<','<<name.second<<','<<std::abs(av-bv)<<','<<difference<<','<<eca<<','<<ecb<<','<<r<<','<<result.judgment<<','<<result.cause<<'\n';
  }
 }
 std::ofstream summary(out+"/coverage.csv");summary<<std::setprecision(17)<<"set,T,comparisons,max_ratio,max_density_ratio,uncovered,ordinary,conditional,failed,largest_density_difference,log_radius_sum_at_largest\n";
 for(const auto& [key,g]:groups)summary<<key.first<<','<<key.second<<','<<g.count<<','<<g.maximum<<','<<g.density_max<<','<<g.uncovered<<','<<g.ordinary<<','<<g.conditional<<','<<g.failed<<','<<g.largest_density_difference<<','<<g.largest_pair_radius<<'\n';
 std::map<double,std::array<__float128,2>> references;
 for(const auto& r:read(out+"/binary128-roots.csv"))if(r.at("method")=="binary128")references[number(r,"T_K")]={strtoflt128(r.at("rho_L").c_str(),nullptr),strtoflt128(r.at("rho_V").c_str(),nullptr)};
 std::ofstream absolute(out+"/absolute-reference.csv");absolute<<std::setprecision(17)<<"row,arm,T,phase,absolute_error,absolute_radius,ratio,log_error,log_radius,log_ratio\n";
 for(const auto& [key,q]:endpoints)if(std::get<0>(key)=="854")for(int phase=0;phase<2;++phase) {
  const auto& r=q.at({"density",phase});const double t=number(r,"T"),value=number(r,"value"),radius=number(r,"estimate_abs");
  const __float128 ref=references.at(t)[phase];const double error=double(fabsq(__float128(value)-ref));
  const double log_error=double(fabsq(logq(__float128(value)/ref))),log_radius=number(r,"estimate_chart");
  absolute<<std::get<1>(key)<<','<<std::get<2>(key)<<','<<t<<','<<phase<<','<<error<<','<<radius<<','<<error/radius<<','<<log_error<<','<<log_radius<<','<<log_error/log_radius<<'\n';
 }
 std::cout<<"endpoints="<<endpoints.size()<<", scalar_rows="<<rows.size()<<", groups="<<groups.size()<<'\n';
 return 0;
}

int main(int argc,char** argv) {try {run(argc,argv);return 0;} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
