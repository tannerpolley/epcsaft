#include "equilibrium/gibbs_minimizer.hpp"
#include "one_component_two_density_inputs.hpp"
#include "qualification_clock.hpp"
#include <iostream>
#include <iomanip>
using namespace epcsaft; using namespace epcsaft::equilibrium;
int main(int argc, char** argv) {
 auto ethane=properties_support::parameters(1);ethane.segment_count={1.6069};ethane.sigma_constant={3.5206};ethane.dispersion_energy_over_k={191.42};ethane.molar_mass_kg_per_mol={0.03007};
 std::cout<<std::setprecision(20)<<"case,region,repetition,calls,elapsed,process_cpu,value,coefficient,error,iterations,evaluations,hessians,bits,corrections,site_factorizations,site_solves,pressure_factorizations,pressure_solves,thermal_evaluations,recordings,directional_paths,forward_sweeps,fixed0,pressure0,pressure1,pressure_floor0,pressure_floor1,temperature_coefficient\n";
 for(const auto& [name,parameters,T]:std::vector<std::tuple<std::string,Parameters,double>>{
  {"propane",caloric_specification_checks::propane_parameters(),300},
  {"propane-low-pressure",caloric_specification_checks::propane_parameters(),110},
  {"propane-packing-conversion",caloric_specification_checks::propane_parameters(),216.4},
  {"ethane-source",ethane,280},{"acetone",pure_saturation_checks::acetone(),300},{"methanol",pure_saturation_checks::methanol(),300}}) {
  if(argc>1 && name!=argv[1])continue;
  QualificationClock construction; HelmholtzModel model(parameters);
  auto compiled=compile_problem(model,pure_saturation_checks::saturation(T));
  std::cout<<name<<",construction,0,1,"<<construction.seconds()<<','<<construction.cpu_seconds()<<",nan,nan,nan,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable\n";
  EquilibriumResult endpoint;
  for(int repetition=-3;repetition<=7;++repetition)for(int control=0;control<2;++control) {
   QualificationClock clock;std::size_t calls=0;
   do {endpoint=solve_saturation(compiled);++calls;}while(repetition>0 && clock.seconds()<1.0);
   const double elapsed=clock.seconds(),cpu=clock.cpu_seconds();
   if(!endpoint.success)throw std::runtime_error(endpoint.message);
   std::cout<<name<<(control?"_B":"_A")<<",complete,"<<repetition<<','<<calls<<','<<elapsed<<','<<cpu<<','<<*endpoint.pressure<<','<<endpoint.molar_densities[0]<<','<<endpoint.molar_densities[1]<<','<<endpoint.iterations<<','<<endpoint.evaluations<<','<<endpoint.hessian_calls<<','<<(endpoint.linear_solve_diagnostics?endpoint.linear_solve_diagnostics->precision_bits:53)<<','<<(endpoint.linear_solve_diagnostics?endpoint.linear_solve_diagnostics->refinement_iterations:0)<<",unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable\n";
  }
  const std::vector<double> logs{std::log(endpoint.molar_densities[0]),std::log(endpoint.molar_densities[1])};
  const auto prepared=one_component_two_density_point(compiled,logs);
  // This enclosing interval includes the shared Base=double reduction's wider accumulation.
  for(int repetition=-3;repetition<=7;++repetition)for(int control=0;control<2;++control) {
   QualificationClock clock;std::size_t calls=0;OneComponentTwoDensityPoint value;
   do {value=one_component_two_density_point(compiled,logs,nullptr,&prepared.phase,std::log(model.packing_fraction_coefficient(T,{1.0})));++calls;}while(repetition>0 && clock.seconds()<1.0);
   const double elapsed=clock.seconds(),cpu=clock.cpu_seconds();
   std::cout<<name<<(control?"_B":"_A")<<",two_density_enclosing,"<<repetition<<','<<calls<<','<<elapsed<<','<<cpu<<','<<value.rows[0]<<','<<value.rows[1]<<",nan,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable,unavailable\n";
  }
  // Same pure-state insertion and thermal direction, paired with optional error capture.
  PressureActionDirection direction;direction.temperature=T;
  const double volume=name=="ethane-source" ? 1/model.state_at_pressure(T,3e6,{1.0},DensityBranch::Liquid).molar_density : 1/endpoint.molar_densities[0];
  for(int repetition=-3;repetition<=7;++repetition)for(int control=0;control<2;++control)for(int arm=0;arm<2;++arm) {
   const bool capture=(arm+(repetition%2+2)%2)%2;ReferenceTaylorError error;ReferenceTaylor value;
   QualificationClock clock;std::size_t calls=0;
   do{value=model.reference_taylor(T,volume,{1.0},{1.0},direction,1,false,capture?&error:nullptr);++calls;}while(repetition>0 && clock.seconds()<1.0);
   const double elapsed=clock.seconds(),cpu=clock.cpu_seconds();
   if(value.status!=PropertyActionStatus::Available)throw std::runtime_error("reference Taylor unavailable");
   std::cout<<name<<(control?"_B":"_A")<<(capture?",reference_error,":",reference_control,")<<repetition<<','<<calls<<','<<elapsed<<','<<cpu<<','<<value.contraction[0]<<','<<value.contraction[1]<<','<<(capture?error.fixed[1]:NAN)<<",unavailable,unavailable,unavailable,unavailable,unavailable";
   const auto& d=value.diagnostics;
   std::cout<<','<<d.site_factorizations<<','<<d.site_solves<<','<<d.pressure_factorizations<<','<<d.pressure_solves<<','<<d.thermal_evaluations<<','<<d.recordings<<','<<d.directional_paths<<','<<d.forward_sweeps;
   if(capture)std::cout<<','<<error.fixed[0]<<','<<error.pressure[0]<<','<<error.pressure[1]<<','<<error.pressure_floor[0]<<','<<error.pressure_floor[1]<<','<<error.temperature_coefficient;
   else std::cout<<",unavailable,unavailable,unavailable,unavailable,unavailable,unavailable";
   std::cout<<'\n';
  }
 }
}
