#include "surface/parameter_builder.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "contributions/helmholtz.hpp"
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <numeric>
namespace nb = nanobind;
using namespace nb::literals;
using namespace epcsaft;
using namespace epcsaft::equilibrium;
template <class... A> nb::dict data(A&&... a) {
    return nb::cast<nb::dict>(nb::module_::import_("builtins").attr("dict")(std::forward<A>(a)...));
}
nb::dict run(nb::dict input, int solvent, double w, double measured, double published) {
    auto a = nb::cast<nb::dict>(input["assumptions"]), raw = nb::cast<nb::dict>(input["parameters"]);
    const double T = nb::cast<double>(a["temperature_K"]), P = nb::cast<double>(a["pressure_Pa"]);
    std::vector<int> indices{0, solvent, 3, 4};
    auto values = [&](const char* key) {
        std::vector<double> v; for (auto x : raw[key]) v.push_back(x.is_none() ? 0 : nb::cast<double>(x)); return v;
    };
    auto matrix = [&](const char* key) { return nb::cast<std::vector<std::vector<double>>>(raw[key]); };
    const auto mass = values("molar_mass_kg_mol"), energy = values("association_energy_over_k_k"),
               volume = values("association_volume");
    const auto kij = matrix("k_ij"), khb = matrix("k_hb_ij");
    ParameterDocument d;
    d.component_ids = {"water", "solvent", "Li+", "Cl-"};
    d.base_choice = "pc-saft"; d.association_choice = "none"; d.electrolyte_choice = "born";
    d.ionic_dispersion_choice = "all-ion-pairs"; d.permittivity_choice = "component-permittivity-mixing";
    d.polar_choice = "none";
    auto entry = [&](std::string field, std::vector<std::size_t> index, double v, std::string unit) {
        d.entries.push_back({field + std::to_string(d.entries.size()), field, index, v, unit});
    };
    for (std::size_t i = 0; i < 4; ++i) {
        const int j = indices[i];
        auto sigma = nb::cast<nb::dict>(raw["segment_diameter_angstrom"][j]);
        auto eps = nb::cast<nb::dict>(raw["relative_permittivity"][j]);
        const double s = nb::cast<double>(sigma[i == 0 ? "constant" : "value"]);
        entry("molar_mass", {i}, mass[j], "kg / mol");
        entry("solvation_factor", {i}, 0, "dimensionless");
        entry("segment_count", {i}, values("segment_count")[j], "dimensionless");
        entry("sigma_constant", {i}, s, "angstrom");
        entry("dispersion_energy_over_k", {i}, values("dispersion_energy_over_k_k")[j], "kelvin");
        entry("charge_number", {i}, i < 2 ? 0 : i == 2 ? 1 : -1, "elementary-charge");
        entry("relative_permittivity", {i}, i ? nb::cast<double>(eps["value"]) :
              nb::cast<double>(eps["a"]) * std::log(T) + nb::cast<double>(eps["b"]), "dimensionless");
        if (i >= 2) for (auto field : {"packing_diameter", "debye_huckel_diameter", "born_diameter"})
            entry(field, {i}, s * (std::string(field) == "packing_diameter" ? 1 : .88), "angstrom");
        for (std::size_t k = i + 1; k < 4; ++k) entry("k_ij", {i, k}, kij[j][indices[k]], "dimensionless");
    }
    for (auto field : {"ionic_region_relative_permittivity", "born_c_shell", "born_c_dielectric"})
        entry(field, {}, std::string(field) == "ionic_region_relative_permittivity" ? 8 : 0, "dimensionless");
    auto p = build_parameters(d);
    for (auto term : raw["segment_diameter_angstrom"][0]["terms"]) {
        p.sigma_term_component.push_back(0);
        p.sigma_term_amplitude.push_back(nb::cast<double>(term["amplitude"]));
        p.sigma_term_exponent.push_back(nb::cast<double>(term["exponent_coefficient_per_k"]));
    }
    p.site_component = {0,0,1,1,2,2}; p.site_multiplicity.assign(6,1);
    p.association_active.assign(36,0); p.association_energy_over_k.assign(36,0); p.association_volume.assign(36,0);
    const auto sigma = segment_sigma(T,p);
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        const double e = (energy[indices[i]] + energy[indices[j]]) / 2 * (1 - khb[indices[i]][indices[j]]);
        if (e <= 0) continue;
        const double v = std::sqrt(volume[indices[i]] * volume[indices[j]]) *
                         std::pow(std::sqrt(sigma[i]*sigma[j]) / ((sigma[i]+sigma[j])/2),3);
        for (int slot : {(2*i)*6+2*j+1,(2*j+1)*6+2*i}) {
            p.association_active[slot]=1; p.association_energy_over_k[slot]=e; p.association_volume[slot]=v;
        }
    }
    HelmholtzModel model(p,{});
    std::vector<double> physical{(1-w)*nb::cast<double>(a["aqueous_solution_mass_kg"])/mass[0],
        nb::cast<double>(a["pure_solvent_mass_kg"])/mass[solvent],
        w*nb::cast<double>(a["aqueous_solution_mass_kg"])/(mass[3]+mass[4])};
    physical.push_back(physical[2]);
    auto feed = resolve_feed(model,{FeedKind::Amounts,physical});
    const double total = std::accumulate(feed.begin(),feed.end(),0.0);
    for (auto& n : feed) n /= total;
    EquilibriumProblem problem{{},FreeValue::fixed(T),FreeValue::fixed(P),{FeedKind::Amounts,feed}};
    problem.phases = {{"aqueous",PhaseKind::Liquid,{0,1,2,3},true,{}, {.98,.0001,.00995,.00995}},
                      {"organic",PhaseKind::Liquid,{0,1,2,3},true,{}, {.05,.9498,.0001,.0001}}};
    problem.solver_options.set("EQUILIBRIUM_ROW_RESIDUAL",1e-10);
    const auto cell = compile_problem(model,problem);
    const auto s = solve_equilibrium(model,problem);
    nb::dict r = data("physical_feed_mol"_a=physical,"feed_mol"_a=feed,"solved"_a=s,
        "T_K"_a=T,"P_Pa"_a=P,"source_P_Pa"_a=100000,"topology"_a=int(s.topology_event),
        "numerical_limits"_a=std::vector<double>{1e-12,1e-10,1e-10,1e-8,1e-8},
        "convergence"_a=core::describe(s.convergence),"separation_limit"_a=policy::PHASE_COINCIDENCE_DISTANCE.value);
    bool ok = s.has_state() && s.requested_tolerance_met == true && s.topology_event == TopologyEvent::None;
    r["admitted"] = false;
    if (!s.has_state()) return r;
    const auto raw_rows = cell.evaluate(s.coordinates).values;
    nb::list phases; std::vector<double> balance=feed, water_mass, observable;
    double chemical=0, distance=0;
    for (std::size_t i=0;i<raw_rows.size();++i)
        if (cell.rows[i].kind==RowKind::Chemical) chemical=std::max(chemical,std::abs(raw_rows[i]));
    for (std::size_t phase=0;phase<2;++phase) {
        std::vector<double> x(s.mole_fractions.begin()+phase*4,s.mole_fractions.begin()+(phase+1)*4);
        const auto state=model.state(T,s.molar_densities[phase],x);
        double mean=0; for (int i=0;i<4;++i) {
            mean+=x[i]*mass[indices[i]]; balance[i]-=s.amounts[phase]*x[i];
            distance=std::max(distance,std::abs(s.mole_fractions[i]-s.mole_fractions[4+i]));
        }
        const double norm=std::accumulate(x.begin(),x.end(),0.0)-1, charge=x[2]-x[3], pressure=std::abs(state.pressure/P-1);
        water_mass.push_back(x[0]*mass[0]/mean);
        observable.push_back(solvent==2 ? x[1]*nb::cast<double>(a["emim_cation_molar_mass_kg_mol"])/mean : x[1]);
        phases.append(data("x"_a=x,"amount_mol"_a=s.amounts[phase],"rho_mol_m3"_a=s.molar_densities[phase],
            "water_mass_fraction"_a=water_mass.back(),"normalization_defect"_a=norm,"charge_mol_per_mol"_a=charge,
            "pressure_Pa"_a=state.pressure,"pressure_relative_defect"_a=pressure,"pressure_slope"_a=state.pressure_density_derivative));
        ok = ok && std::isfinite(s.amounts[phase]) && s.amounts[phase]>0 && *std::min_element(x.begin(),x.end())>=0 &&
             std::abs(norm)<=1e-12 && std::abs(charge)<=1e-10 && pressure<=1e-8 && std::isfinite(state.pressure_density_derivative) && state.pressure_density_derivative>0;
    }
    const int aq=water_mass[0]>water_mass[1] ? 0 : 1;
    const double value=observable[solvent==2 ? aq : 1-aq];
    r.attr("update")(data("minus_measured"_a=value-measured,"minus_published"_a=value-published,
        "absolute_measured_error"_a=std::abs(value-measured),"absolute_published_error"_a=std::abs(value-published)));
    ok=ok && maximum_abs(balance)<=1e-10 && chemical<=1e-8 && std::isfinite(chemical) && distance>policy::PHASE_COINCIDENCE_DISTANCE.value;
    r.attr("update")(data("admitted"_a=ok,"phases"_a=phases,"balance_defects_mol"_a=balance,"raw_rows"_a=raw_rows,
        "chemical_max_over_RT"_a=chemical,"phase_distance"_a=distance,"aqueous_phase"_a=aq,"candidate_value"_a=observable[solvent==2 ? aq : 1-aq]));
    return r;
}
NB_MODULE(subsystems,m) { m.def("run",&run); }
