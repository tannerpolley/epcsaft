#include "surface/parameter_builder.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "regression/observations.hpp"
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <numeric>
namespace nb = nanobind;
using namespace nb::literals;
using namespace epcsaft;
using namespace epcsaft::equilibrium;
using namespace epcsaft::regression;
template <class... A> nb::dict data(A&&... a) {
    return nb::cast<nb::dict>(nb::module_::import_("builtins").attr("dict")(std::forward<A>(a)...));
}
nb::dict run(nb::dict a, nb::list tuples, nb::dict row) {
    const bool mixture=nb::cast<std::string>(row["table_label"])=="S6";
    auto number = [](nb::handle h) { return std::stod(nb::cast<std::string>(h)); };
    const double T=number(row["T_K"]), measured=number(row["density_kg_m3"]);
    ParameterDocument d;
    d.base_choice="pc-saft"; d.association_choice="general-site";
    d.electrolyte_choice=d.permittivity_choice=d.polar_choice="none";
    d.component_ids = mixture ? std::vector<std::string>{"DES","TOPO"} : std::vector<std::string>{"DES"};
    d.site_component = mixture ? std::vector<int>{0,0,1} : std::vector<int>{0,0};
    d.site_multiplicity.assign(d.site_component.size(),1);
    auto entry = [&](std::string field, std::vector<std::size_t> ix, double v, std::string unit) {
        d.entries.push_back({field+std::to_string(d.entries.size()),field,ix,v,unit});
    };
    for (std::size_t i=0;i<d.component_ids.size();++i) {
        auto r=nb::cast<nb::dict>(tuples[i]);
        entry("molar_mass",{i},nb::cast<double>(a["molar_mass_kg_mol"][i ? "topo" : "des"]),"kg / mol");
        entry("charge_number",{i},0,"elementary-charge");
        entry("segment_count",{i},number(r["m"]),"dimensionless");
        entry("sigma_constant",{i},number(r["sigma_A"]),"angstrom");
        entry("dispersion_energy_over_k",{i},number(r["epsilon_over_k_K"]),"kelvin");
        for (auto field : {"association_energy_over_k","association_volume"})
            entry(field,{i ? 2u : 0u,i ? 2u : 1u},number(r[std::string(field)=="association_volume" ? "kappa_assoc" : "epsilon_assoc_over_k_K"]),std::string(field)=="association_volume" ? "dimensionless" : "kelvin");
    }
    if (mixture) {
        entry("k_ij",{0,1},nb::cast<double>(a["DES_TOPO_dispersion_k_ij"]),"dimensionless");
        entry("no_cross_association",{0,1},0,"dimensionless");
    }
    HelmholtzModel model(build_parameters(d),{});
    std::vector<double> w = mixture ? std::vector<double>{nb::cast<double>(a["DES_mass_fraction_S6"]),
        nb::cast<double>(a["TOPO_mass_fraction_S6"])} : std::vector<double>{1};
    auto x=resolve_feed(model,{FeedKind::MassFractions,w});
    const double sum=std::accumulate(x.begin(),x.end(),0.0); for (auto& v:x) v/=sum;
    const auto& masses=model.parameters().molar_mass_kg_per_mol;
    const double mean=std::inner_product(x.begin(),x.end(),masses.begin(),0.0), P=nb::cast<double>(a["pressure_Pa"]);
    std::vector<double> inverse; for (std::size_t i=0;i<x.size();++i) inverse.push_back(x[i]*masses[i]/mean-w[i]);
    Observation o{ObservationKind::Density,T,P,x}; o.branch=DensityBranch::Liquid; o.mass_basis=true;
    const auto p=predict(model,o);
    nb::dict r=data("admitted"_a=false,"status"_a=observation_status_name(p.status),"cause"_a=p.message,
        "x"_a=x,"mass_fraction_defects"_a=inverse,"mean_mass_kg_mol"_a=mean,"P_Pa"_a=P,
        "numerical_limits"_a=std::vector<double>{1e-8,1e-12},"observed_kg_m3"_a=measured);
    if (!p.state || !p.value) return r;
    const auto& s=*p.state; const auto direct=model.state(T,s.molar_density,x);
    const double pressure=std::abs(direct.pressure/P-1), identity=std::abs(*p.value/(s.molar_density*mean)-1);
    const bool ok=p.status==ObservationStatus::Available && std::isfinite(*p.value) && *p.value>0 && s.molar_density>0 &&
        s.density_convergence && s.density_convergence->criterion && s.density_convergence->residual<=s.density_convergence->criterion->value &&
        (s.density_branch==DensityBranch::Liquid || (s.density_branch==DensityBranch::Unspecified && s.stable_root_count==1)) && s.stable_root_count.value_or(0)>0 && direct.pressure_density_derivative>0 &&
        std::isfinite(direct.pressure_density_derivative) && pressure<=1e-8 && identity<=1e-12 && maximum_abs(inverse)<=1e-12 &&
        std::abs(std::accumulate(x.begin(),x.end(),0.0)-1)<=1e-12;
    r.attr("update")(data("admitted"_a=ok,"density_kg_m3"_a=*p.value,"rho_mol_m3"_a=s.molar_density,
        "branch"_a=int(s.density_branch),"stable_roots"_a=s.stable_root_count.value_or(0),"pressure_Pa"_a=direct.pressure,
        "pressure_relative_defect"_a=pressure,"pressure_slope"_a=direct.pressure_density_derivative,"mass_density_defect"_a=identity,
        "convergence"_a=core::describe(*s.density_convergence),"minus_measured_kg_m3"_a=*p.value-measured,
        "relative_error_percent"_a=100*(*p.value/measured-1)));
    PressureActionDirection direction{1,0,std::vector<double>(x.size())};
    const auto action=model.pressure_actions(T,P,x,DensityBranch::Liquid,{}, {direction},s.molar_density);
    r["action_cause"]=action.density_message; r["action_status"]=int(action.density_status);
    if (action.density_status==PropertyActionStatus::Available && !action.density_action.empty()) {
        const double derivative=action.density_action[0]*mean;
        r.attr("update")(data("density_T_kg_m3_K"_a=derivative,"action_rho_mol_m3"_a=*action.molar_density,
            "local_uncertainty_kg_m3"_a=std::hypot(number(row["density_uncertainty_kg_m3"]),number(row["temperature_uncertainty_K"])*derivative)));
    }
    return r;
}
NB_MODULE(density,m) { m.def("run",&run); }
