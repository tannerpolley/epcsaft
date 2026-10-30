#include "equilibrium/gibbs_minimizer.hpp"
#include "one_component_two_density_inputs.hpp"
#include <iostream>
#include <iomanip>
using namespace epcsaft;using namespace epcsaft::equilibrium;
int main(){HelmholtzModel m(caloric_specification_checks::propane_parameters());std::cout<<std::setprecision(20)<<"T,rhoL,rhoV,radiusL,radiusV,success,bits,corrections\n";for(double T:{375.13,375.139}){core::Capture capture(true);auto r=solve_equilibrium(m,pure_saturation_checks::saturation(T));if(!r.success||!r.numerical_estimates)throw std::runtime_error(r.message);std::cout<<T<<','<<r.molar_densities[0]<<','<<r.molar_densities[1]<<','<<r.numerical_estimates->densities[0].chart<<','<<r.numerical_estimates->densities[1].chart<<','<<r.success<<','<<(r.linear_solve_diagnostics?r.linear_solve_diagnostics->precision_bits:53)<<','<<(r.linear_solve_diagnostics?r.linear_solve_diagnostics->refinement_iterations:0)<<'\n';}}
