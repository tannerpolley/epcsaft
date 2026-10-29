#include "equilibrium/gibbs_minimizer.hpp"
#include "properties_support.hpp"
#include "properties/helmholtz_tape.hpp"
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace epcsaft;
using namespace epcsaft::equilibrium;
int main() {
 const double t=330,v=6e-5;
 const std::vector<double> x{.3,.7};const auto p=properties_support::associating_binary();
 HelmholtzModel model(p,{},properties_support::p3_binary_record());
 std::cout<<std::setprecision(17)<<"case,expected_inner,produced_inner,covered\n";
 bool passed=true;
 const auto check=[&](const std::string& label,double expected,double produced) {
  const bool ok=produced==expected;passed=passed&&ok;
  std::cout<<label<<','<<expected<<','<<produced<<','<<ok<<'\n';
 };
 for(int kind=0;kind<4;++kind) {
  EquilibriumProblem problem;problem.feed.values=x;problem.pressure=FreeValue::fixed(1e5);
  problem.temperature=kind==0?FreeValue::fixed(t):FreeValue::variable(t);
  PhaseSpec phase;phase.name="L";phase.kind=PhaseKind::Liquid;phase.support={0,1};problem.phases={phase};
  if(kind==1)problem.enthalpy_j=0;if(kind==2)problem.entropy_j_per_k=0;
  if(kind==3){problem.pressure=FreeValue::variable(1e5);problem.critical=CriticalDeclaration{1/v,{}};}
  const auto c=compile_problem(model,problem);auto point=c.generic_point();
  const auto offsets=c.phase_offsets();point[offsets[0]]=std::log(x[0]);point[offsets[0]+1]=std::log(x[1]);
  point[offsets[0]+2]=1;point[offsets[0]+3]=std::log(model.packing_fraction_coefficient(t,x)/v);
  auto at=c.prepare(point);at.derivative[0].site_fractions[0]+=1e-7;
  const auto u=c.row_uncertainty(at);if(!u.cause.empty())throw std::runtime_error(u.cause);
  std::vector<double> base{at.temperature,at.volume[0],x[0],x[1]},g(4);
  double energy=0;
  for(std::size_t j=0;j<4;++j){std::vector<double> unit(4,0);unit[j]=1;auto e=model.residual_coefficient_errors({base,unit},at.derivative[0].site_fractions);energy=e[0];g[j]=e[1];}
  for(std::size_t i=0;i<c.rows.size();++i) {
   const auto& row=c.rows[i];
   if(row.kind==RowKind::Chemical)check("chemical_"+std::to_string(kind)+"_"+std::to_string(row.index),g[2+row.index],u.inner[i]);
   if(row.kind==RowKind::Pressure)check("pressure_"+std::to_string(kind),std::abs(at.volume[0])*g[1],u.inner[i]);
   if(row.kind==RowKind::Enthalpy)check("enthalpy",at.amount[0]/c.feed_total*(at.temperature*g[0]+std::abs(at.volume[0])*g[1]),u.inner[i]);
   if(row.kind==RowKind::Entropy)check("entropy",at.amount[0]/c.feed_total*(at.temperature*g[0]+energy),u.inner[i]);
   if(row.kind==RowKind::CriticalMode||row.kind==RowKind::CriticalCubic) {
    std::vector<double> w(4,0);for(std::size_t j=0;j<2;++j)w[2+j]=std::sqrt(c.n_feed[j]/c.feed_total)*at.point[c.mode_offset()+j];
    const auto error=[&](const std::vector<double>& dir,std::size_t order) {std::vector<std::vector<double>> path(order+1,std::vector<double>(4,0));path[0]=base;path[1]=dir;return model.residual_coefficient_errors(path,at.derivative[0].site_fractions)[order];};
    if(row.kind==RowKind::CriticalCubic)check("critical_cubic",6*error(w,3),u.inner[i]);
    else {std::vector<double> unit(4,0),sum=w;unit[2+row.index]=1;sum[2+row.index]+=1;check("critical_mode_"+std::to_string(row.index),std::sqrt(c.n_feed[row.index]/c.feed_total)*(error(sum,2)+error(unit,2)+error(w,2)),u.inner[i]);}
   }
  }
 }
 return passed?0:1;
}
