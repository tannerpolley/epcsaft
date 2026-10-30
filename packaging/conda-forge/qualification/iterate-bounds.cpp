#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_support.hpp"
#include "properties_support.hpp"
#include <iostream>
#include <iomanip>
using namespace epcsaft;using namespace epcsaft::equilibrium;
int main(){HelmholtzModel model(properties_support::neutral_binary());auto problem=equilibrium_support::bubble_at(200);problem.feed.values={0.3,0.7};const auto c=compile_problem(model,problem);GibbsMinimizer nlp(c,c.generic_point());Ipopt::Index n,m,j,h;Ipopt::TNLP::IndexStyleEnum style;if(!nlp.get_nlp_info(n,m,j,h,style))return 1;std::vector<double> lo(n),hi(n);if(!nlp.physical_bounds(lo,hi))return 1;std::cout<<std::setprecision(17)<<"coordinate,lower,upper\n";for(std::size_t k=0;k<lo.size();++k){std::cout<<k<<','<<lo[k]<<','<<hi[k]<<'\n';if(std::isnan(lo[k])||std::isnan(hi[k])||lo[k]>hi[k])return 1;}core::Capture capture(true);core::Scope scope("windows_iterate_bounds");auto r=run_driver(model,c,boundary_start(c),structural_result(c,scope),scope);std::cout<<"ipopt_diagnostics,"<<r.ipopt_diagnostics.has_value()<<"\n";for(std::size_t k=0;k<r.coordinates.size();++k)std::cout<<"returned_coordinate,"<<k<<","<<r.coordinates[k]<<"\n";std::cout<<"iterations,"<<r.iterations<<",evaluations,"<<r.evaluations<<",pressure,"<<*r.pressure<<'\n';return r.success&&r.requested_tolerance_met.value_or(false)&&r.ipopt_diagnostics.has_value()?0:1;}
