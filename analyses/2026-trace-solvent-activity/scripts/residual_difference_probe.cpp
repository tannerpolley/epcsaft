#include "surface/parameter_builder.hpp"
#include "properties/helmholtz_model.hpp"
#include "core/physical_constants.hpp"
#include <iostream>
#include <iomanip>
int main(){ epcsaft::ParameterDocument d;
d.component_ids={"water","sodium-cation","chloride-anion"};
d.site_component={0,0};
d.site_multiplicity={1,1};
d.base_choice="pc-saft";
d.association_choice="general-site";
d.electrolyte_choice="fully-dissociated-debye-huckel";
d.ionic_dispersion_choice="exclude-same-sign-ion-pairs";
d.permittivity_choice="solvent-only";
d.polar_choice="none";
d.entries.push_back({"component/water/molar_mass","molar_mass",{0},18.01528,"gram / mole"});
d.entries.push_back({"component/water/charge_number","charge_number",{0},0.0,"elementary-charge"});
d.entries.push_back({"component/water/dispersion_energy_over_k","dispersion_energy_over_k",{0},353.9449,"kelvin"});
d.entries.push_back({"component/water/segment_count","segment_count",{0},1.204659,"dimensionless"});
d.entries.push_back({"component/sodium-cation/charge_number","charge_number",{1},1.0,"elementary-charge"});
d.entries.push_back({"component/sodium-cation/debye_huckel_diameter","debye_huckel_diameter",{1},2.122736,"angstrom"});
d.entries.push_back({"component/sodium-cation/dispersion_energy_over_k","dispersion_energy_over_k",{1},646.0504,"kelvin"});
d.entries.push_back({"component/sodium-cation/packing_diameter","packing_diameter",{1},2.122736,"angstrom"});
d.entries.push_back({"component/sodium-cation/segment_count","segment_count",{1},1.0,"dimensionless"});
d.entries.push_back({"component/sodium-cation/segment_diameter","sigma_constant",{1},2.4122,"angstrom"});
d.entries.push_back({"component/chloride-anion/charge_number","charge_number",{2},-1.0,"elementary-charge"});
d.entries.push_back({"component/chloride-anion/debye_huckel_diameter","debye_huckel_diameter",{2},2.6906,"angstrom"});
d.entries.push_back({"component/chloride-anion/dispersion_energy_over_k","dispersion_energy_over_k",{2},47.2878,"kelvin"});
d.entries.push_back({"component/chloride-anion/packing_diameter","packing_diameter",{2},2.6906,"angstrom"});
d.entries.push_back({"component/chloride-anion/segment_count","segment_count",{2},1.0,"dimensionless"});
d.entries.push_back({"component/chloride-anion/segment_diameter","sigma_constant",{2},3.0575,"angstrom"});
d.entries.push_back({"component/water/relative_permittivity/constant/constant","relative_permittivity",{0},78.09,"dimensionless"});
d.entries.push_back({"component/water/segment_diameter/constant-plus-sum-of-exponentials/constant","sigma_constant",{0},2.7927,"angstrom"});
d.entries.push_back({"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/amplitude","sigma_term_amplitude",{0,0},10.11,"angstrom"});
d.entries.push_back({"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/exponent_coefficient","sigma_term_exponent",{0,0},-0.01775,"1 / kelvin"});
d.entries.push_back({"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/amplitude","sigma_term_amplitude",{1,0},-1.417,"angstrom"});
d.entries.push_back({"component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/exponent_coefficient","sigma_term_exponent",{1,0},-0.01146,"1 / kelvin"});
d.entries.push_back({"pair/chloride-anion/sodium-cation/k_ij","k_ij",{2,1},1.0,"dimensionless"});
d.entries.push_back({"pair/chloride-anion/water/k_ij","k_ij",{2,0},0.0,"dimensionless"});
d.entries.push_back({"pair/sodium-cation/water/k_ij","k_ij",{1,0},0.0,"dimensionless"});
d.entries.push_back({"association/water/a/water/b/energy_over_k","association_energy_over_k",{0,1},2425.6714,"kelvin"});
d.entries.push_back({"association/water/a/water/b/volume","association_volume",{0,1},0.0450989,"dimensionless"});

auto model=epcsaft::HelmholtzModel(epcsaft::build_parameters(d));
const double T=298.15,P=100000,MW=.01801528,rt=epcsaft::GAS_CONSTANT_J_PER_MOL_K*T;
auto pure=model.state_at_pressure(T,P,{1,0,0},epcsaft::DensityBranch::Liquid,std::nullopt);
std::cout<<std::setprecision(17);
for(double m:{1e-3,1e-6,1e-8,1e-10,1e-12}){
 auto x=epcsaft::formula_composition({1,0,0},model.parameters(),{1,2},{1,1},m);
 auto f=model.state_at_pressure(T,P,x,epcsaft::DensityBranch::Liquid,std::nullopt);
 double s=2*m*MW;
 double dv=(s*pure.molar_density-(f.molar_density-pure.molar_density))/(f.molar_density*pure.molar_density);
 auto integ=[&](bool four){
  double sum=0;
  std::vector<double> nodes=four?std::vector<double>{.06943184420297371,.33000947820757187,.6699905217924281,.9305681557970262}:std::vector<double>{.21132486540518713,.7886751345948129};
  std::vector<double> weights=four?std::vector<double>{.17392742256872692,.32607257743127305,.32607257743127305,.17392742256872692}:std::vector<double>{.5,.5};
  for(size_t i=0;i<nodes.size();++i){
   double u=nodes[i],t=u*u;
   std::vector<double> amounts{1,t*m*MW,t*m*MW};
   auto h=model.evaluate(T,1/pure.molar_density+t*dv,amounts,epcsaft::DerivativeOrder::Hessian);
   size_t n=h.gradient.size();
   double v=h.hessian[2*n+1]*dv;
   for(size_t j=1;j<3;++j)v+=h.hessian[2*n+2+j]*m*MW;
   sum+=weights[i]*2*u*v;
  }
  return sum;
 };
 auto corr=[&](const epcsaft::State& st,const std::vector<double>& amounts){
  auto h=model.evaluate(T,1/st.molar_density,amounts,epcsaft::DerivativeOrder::Hessian);
  double rho=st.molar_density;
  double a=1/rho-h.hessian[2*h.gradient.size()+1]/(rho*rho);
  double correction=a*(P-st.pressure)/st.pressure_density_derivative;
  double floor=a*rt*rho*st.density_convergence->floor.value_or(0)/st.pressure_density_derivative;
  return std::array<double,3>{correction,floor,h.gradient_floor};
 };
 auto c=corr(f,x),c0=corr(pure,{1,0,0});
 double logx=-std::log1p(s), lr=std::log1p((f.molar_density-pure.molar_density)/pure.molar_density);
 double mu=integ(true),err=std::abs(mu-integ(false));
 std::cout<<m<<","<<-(logx+mu+lr)/s<<","<<-(logx+mu+lr+c[0]-c0[0])/s<<","<<err/s<<","<<(std::abs(c[1])+std::abs(c0[1]))/s<<","<<c[2]<<","<<c0[2]<<"\n";
}
}
