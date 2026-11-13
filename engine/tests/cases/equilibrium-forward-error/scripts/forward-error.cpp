#include "equilibrium/gibbs_minimizer.hpp"
#include "properties_support.hpp"
#include "core/trace.hpp"
#include <Eigen/Dense>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace epcsaft;
using namespace epcsaft::equilibrium;
using Dense = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;


int main(int argc,char**argv) {
 std::cout<<std::setprecision(17);
 std::ofstream coordinates(argv[3]); coordinates<<std::setprecision(17)<<"set,row,arm,index,value\n";
 std::ofstream floors(argv[1]); floors<<std::setprecision(17)<<"set,row,arm,row_index,kind,term_scale,inner,floor\n";
 std::ofstream packing(argv[2]); packing<<std::setprecision(17)<<"set,row,arm,phase,packing,log_radius,lower,upper\n";
 std::cout<<"set,row,arm,T,P,residual,n,rank,condition,inverse_inf,quantity,index,value,estimate_abs,estimate_chart\n";
 std::string line;
 while(std::getline(std::cin,line)) {
  if(line.empty())continue;
  std::istringstream in(line);std::string set,arm;int row,components;double t,p;
  in>>set>>row>>arm>>components>>t>>p;
  std::vector<double>x(2*components),rho(2);
  for(auto&v:x)in>>v;for(auto&v:rho)in>>v;
  auto par=components==2?properties_support::neutral_binary():properties_support::parameters(1);
  if(components==1){par.segment_count={2.002};par.sigma_constant={3.6184};par.dispersion_energy_over_k={208.11};par.molar_mass_kg_per_mol={0.044097};}
  HelmholtzModel model(par);EquilibriumProblem problem;
  PhaseSpec liquid;liquid.name="L";liquid.kind=PhaseKind::Liquid;
  for(int i=0;i<components;++i)liquid.support.push_back(i);
  PhaseSpec vapor=liquid;vapor.name="V";vapor.kind=PhaseKind::Vapor;vapor.amount=AmountSpec::pinned(0);
  problem.phases={liquid,vapor};problem.feed.values=std::vector<double>(components,1./components);
  problem.temperature=FreeValue::fixed(t);problem.pressure=FreeValue::variable(p);
  const auto c=compile_problem(model,problem);
  auto point=c.generic_point();point.resize(c.primal_size);
  const auto offsets=c.phase_offsets();
  for(int a=0;a<2;++a){
   for(int i=0;i<components;++i) point[offsets[a]+i]=std::log(x[a*components+i]);
   const auto eta=offsets[a]+components+(a==0?1:0);
   point[eta]=std::log(model.packing_fraction_coefficient(t,{x.begin()+a*components,x.begin()+(a+1)*components})*rho[a]);
   if(a==0)point[offsets[a]+components]=1.;
  }
  point[c.primal_size-1]=p;point=add_multipliers(c,point);
  EquilibriumResult result;
  result.success=true; result.temperature=t; result.pressure=p; result.amounts={1,0};
  result.molar_densities=rho; result.mole_fractions=x; result.coordinates=point;
  result.declaration=std::make_shared<const EquilibriumProblem>(c.problem);result.model_identity=&model;
  for(std::size_t j=0;j<point.size();++j)coordinates<<set<<','<<row<<','<<arm<<','<<j<<','<<point[j]<<'\n';
  const auto at=c.prepare(point);result.packing_fractions=at.eta;
  const auto e=equilibrium_numerical_estimates(c,result);
  if(!e.cause.empty())throw std::runtime_error(e.cause);
  for(std::size_t i=0;i<c.rows.size();++i)
   floors<<set<<','<<row<<','<<arm<<','<<i<<','<<int(c.rows[i].kind)<<','<<e.row_terms[i]<<','<<e.row_inner[i]<<','<<e.row_floor[i]<<'\n';
  for(int phase=0;phase<2;++phase)
   packing<<set<<','<<row<<','<<arm<<','<<phase<<','<<at.eta[phase]<<','<<e.packing[phase].chart<<','<<policy::PACKING_FRACTION_FLOOR.value<<','<<policy::PACKING_FRACTION_CEILING.value<<'\n';
  for(int k=0;k<1+2+2*components+2;++k) {
   std::string q=k==0?"pressure":k<3?"density":k<3+2*components?"composition":"amount";
   int index=k==0?0:k<3?k-1:k<3+2*components?k-3:k-3-2*components;
   double value=k==0?p:k<3?rho[index]:k<3+2*components?x[index]:result.amounts[index];
   const auto& error=k==0?e.pressure:k<3?e.densities[index]:k<3+2*components?e.composition[index]:e.amounts[index];
   std::cout<<set<<','<<row<<','<<arm<<','<<t<<','<<p<<','<<maximum_abs(e.residuals)<<','<<e.dimension<<','<<e.rank<<','<<e.condition<<','<<e.inverse_norm<<','<<q<<','<<index<<','<<value<<','<<error.absolute<<','<<error.chart<<'\n';
  }
 }
}
