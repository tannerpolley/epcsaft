#include "properties_support.hpp"
#include "properties/helmholtz_tape.hpp"
#include <iomanip>
#include <iostream>
int main() {
 using namespace epcsaft;
 const auto p=properties_support::associating_binary();HelmholtzModel m(p);
 const double t=330,v=6e-5;const std::vector<double> composition{0.3,0.7};
 auto d=m.evaluate(t,v,composition,DerivativeOrder::Hessian);
 std::cout<<std::setprecision(17)<<"site_perturbation,coordinate,prescribed_site_floor,reported_inner,covered\n";
 bool passed=true;
 for(double perturb:{0.,1e-7}) {
  auto sites=d.site_fractions;sites[0]+=perturb;
  std::vector<double> point{t,v,0.3,0.7},input=point;input.insert(input.end(),sites.begin(),sites.end());
  auto tape=detail::HelmholtzTape::record(p,{},point,sites,{});
  const auto values=tape->function.Forward(0,input);
  double b0=sites.size()*std::numeric_limits<double>::epsilon();
  for(std::size_t a=0;a<sites.size();++a)b0=std::max(b0,std::abs(values[1+a]));
  const auto sigma=segment_sigma(t,p),diameter=segment_diameter(t,p,sigma);
  const auto lu=properties::site_factorization(sites,composition,p,NUMBER_DENSITY_PER_ANGSTROM3/v,
      association_delta(t,1/v,composition,p,sigma,diameter));
  for(std::size_t j=0;j<4;++j) {
   std::vector<double> unit(input.size(),0);unit[j]=1;
   tape->function.Forward(0,input);const auto first=tape->function.Forward(1,unit);
   Eigen::VectorXd rhs(sites.size());for(std::size_t a=0;a<sites.size();++a)rhs[a]=-first[1+a];
   const Eigen::VectorXd tangent=lu.solve(rhs);
   double norm=0;for(std::size_t a=0;a<sites.size();++a)
    norm+=std::abs(tangent[a]*composition[p.site_component[a]]*p.site_multiplicity[a]/sites[a]);
   std::vector<double> direction(4,0);direction[j]=1;
   const auto e=m.residual_coefficient_errors({point,direction},sites);
   const bool covered=e[1]>=norm*b0;passed=passed&&covered;
   std::cout<<perturb<<','<<j<<','<<norm*b0<<','<<e[1]<<','<<covered<<'\n';
  }
 }
 return passed?0:1;
}
