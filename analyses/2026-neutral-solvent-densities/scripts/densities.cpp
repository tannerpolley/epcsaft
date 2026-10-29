// Fixed pure-solvent comparisons; thermodynamics stays with shared Engine owners.
#include "equilibrium/problem.hpp"
#include "surface/parameter_builder.hpp"
#include "regression/observations.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace epcsaft;
using namespace epcsaft::equilibrium;

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    const std::string input = argv[1], output = argv[2];
    std::ifstream parameters(input + "/parameters.tsv");
    std::string id, cas;
    double mass, m, sigma, eps, ae, av, dipole, count;
    std::map<std::string, Parameters> models;
    std::ofstream readback(output + "/readback.tsv"), failures(output + "/failures.tsv");
    readback << std::setprecision(17) << "component\tCAS\tmass_kg_mol\tm\tsigma_A\tepsilon_K\tassociation_K\tassociation_volume\tdipole_D\tdipole_count\tsites\tpass\n";
    failures << "row\tattempt\tcause\n";
    while (parameters >> id >> cas >> mass >> m >> sigma >> eps >> ae >> av >> dipole >> count) {
        ParameterDocument d;
        d.component_ids = {cas};
        d.base_choice = "pc-saft"; d.electrolyte_choice = "none"; d.permittivity_choice = "none";
        d.association_choice = ae > 0 ? "general-site" : "none";
        d.polar_choice = count > 0 ? "dipole-dipole" : "none";
        d.entries = {{id,"charge_number",{0},0,"elementary-charge"}, {id,"molar_mass",{0},mass,"kg / mol"}, {id,"segment_count",{0},m,"dimensionless"},
            {id,"sigma_constant",{0},sigma,"angstrom"}, {id,"dispersion_energy_over_k",{0},eps,"kelvin"}};
        if (ae > 0) {
            d.site_component = {0,0}; d.site_multiplicity = {1,1};
            d.entries.push_back({id,"association_energy_over_k",{0,1},ae,"kelvin"});
            d.entries.push_back({id,"association_volume",{0,1},av,"dimensionless"});
        }
        if (count > 0) {
            d.entries.push_back({id,"dipole_moment",{0},dipole,"debye"});
            d.entries.push_back({id,"dipole_count",{0},count,"dimensionless"});
        }
        auto p = build_parameters(d);
        const double energy = ae > 0 ? p.association_energy_over_k[1] : 0;
        const double volume = ae > 0 ? p.association_volume[1] : 0;
        const bool match = p.molar_mass_kg_per_mol[0] == mass && p.segment_count[0] == m
            && p.sigma_constant[0] == sigma && p.dispersion_energy_over_k[0] == eps
            && energy == ae && volume == av && p.charge_number[0] == 0
            && (ae == 0 || (p.association_energy_over_k[0] == 0 && p.association_energy_over_k[2] == ae
                && p.association_energy_over_k[3] == 0 && p.association_volume[2] == av))
            && (count > 0 || p.dipole_count.empty())
            && p.site_component == d.site_component && p.site_multiplicity == d.site_multiplicity
            && (count == 0 || (p.dipole_moment_debye[0] == dipole && p.dipole_count[0] == count));
        readback << id << '\t' << cas << '\t' << p.molar_mass_kg_per_mol[0] << '\t' << p.segment_count[0]
            << '\t' << p.sigma_constant[0] << '\t' << p.dispersion_energy_over_k[0] << '\t' << energy
            << '\t' << volume << '\t' << dipole << '\t' << count << '\t' << p.site_count() << '\t' << match << '\n';
        if (!match) failures << id << "\treadback\tparameter mismatch\n";
        models.emplace(id, p);
    }
    std::ifstream densities(input + "/densities.tsv");
    std::string row, component, variant;
    double t, pressure, observed, uncertainty;
    std::ofstream density(output + "/densities.tsv"), summary(output + "/screens.tsv");
    density << std::setprecision(17) << "row\tcomponent\tvariant\tT_K\tP_Pa\tpredicted_kg_m3\terror_kg_m3\terror_percent\tdelta_central_kg_m3\tpressure_relative\tdp_drho\tbranch_mismatch\tnumerical_pass\tengineering_pass\tsource_interval_pass\tinterpretation\n";
    summary << std::setprecision(17) << "source\tmax_error_percent\tengineering_pass\n";
    std::map<std::string, double> maxima, central, closure_max;
    while (densities >> row >> component >> variant >> t >> pressure >> observed >> uncertainty) {
        try {
            const HelmholtzModel model(models.at(component));
            regression::Observation observation;
            observation.temperature=t; observation.pressure=pressure; observation.composition={1};
            observation.branch=DensityBranch::Liquid; observation.mass_basis=true; observation.observed=observed;
            const auto prediction=regression::predict(model,observation);
            if (!prediction.state || !prediction.value) throw std::runtime_error(prediction.message);
            const auto& state=*prediction.state;
            const double rho=*prediction.value;
            const double error = observed > 0 ? rho-observed : NAN, percent = 100*error/observed;
            const double closure = std::abs(state.pressure-pressure)/std::max(pressure,1.);
            const bool numerical = std::isfinite(rho) && rho > 0 && closure <= 1e-8
                && state.pressure_density_derivative > 0 && !state.density_branch_mismatch();
            closure_max["density_pressure_relative"]=std::max(closure_max["density_pressure_relative"],closure);
            if (variant == "central") central[row] = rho;
            density << row << '\t' << component << '\t' << variant << '\t' << t << '\t' << pressure << '\t' << rho
                << '\t' << error << '\t' << percent << '\t' << rho-central.at(row) << '\t' << closure
                << '\t' << state.pressure_density_derivative << '\t' << state.density_branch_mismatch() << '\t' << numerical
                << '\t' << (observed > 0 ? (std::abs(percent)<=1 ? "1" : "0") : "NA") << '\t' << (observed > 0 && uncertainty > 0 ? (std::abs(error)<=uncertainty ? "1" : "0") : "NA") << '\t'
                << (component=="dodecane" ? "technical sulfonated-kerosene carrier represented by n-dodecane (assumption); dependent on the n-dodecane stand-in" : "fixed-input comparison") << '\n';
            if (observed>0 && variant=="central") maxima[component] = std::max(maxima[component], std::abs(percent));
            if (!numerical) failures << row << '\t' << variant << "\tdensity closure failed\n";
        } catch (const std::exception& e) { failures << row << '\t' << variant << '\t' << e.what() << '\n'; maxima[component]=INFINITY; }
    }
    for (const auto& [name, value] : maxima) summary << name << '\t' << value << '\t' << (value<=1) << '\n';
    std::ifstream vapors(input + "/vapors.tsv");
    std::string source;
    double source_t, source_p;
    std::ofstream vapor(output + "/vapors.tsv"), diagnostics(output + "/residuals.tsv");
    vapor << std::setprecision(17) << "row\tsource\tT_K\tsource_P_Pa\tstart_factor\tsuccess\tP_Pa\terror_Pa\terror_percent\trhoL\trhoV\tslopeL\tslopeV\tpressure_relative\tchemical_RT\tmaterial_scaled\tstart_relative\tnumerical_pass\tengineering_pass\n";
    diagnostics << std::setprecision(17) << "row\tstart_factor\tkind\tindex\tvalue\n";
    maxima.clear();
    const HelmholtzModel model(models.at("MIBK"));
    while (vapors >> row >> source >> source_t >> source_p) {
        t = source == "Fuge" ? source_t+273.15 : source_t;
        pressure = source_p * (source == "Fuge" ? 133.322387415 : 1000.);
        double first = NAN;
        for (double factor : {.5,2.}) {
            try {
                EquilibriumProblem problem;
                problem.temperature = FreeValue::fixed(t); problem.pressure = FreeValue::variable(factor*pressure);
                problem.feed = {FeedKind::Amounts,{1},{},1};
                PhaseSpec liquid, gas;
                liquid.name="liquid"; liquid.support={0}; liquid.composition_guess={1};
                gas=liquid; gas.name="vapor"; gas.kind=PhaseKind::Vapor; gas.amount=AmountSpec::pinned(0);
                problem.phases={liquid,gas};
                const auto result = solve_equilibrium(model, problem);
                for (const auto& [name, values] : std::vector<std::pair<std::string,std::vector<double>>>{
                    {"raw",result.residuals},{"classified",result.classified_residuals},{"amounts",result.amounts},
                    {"fractions",result.mole_fractions},{"densities",result.molar_densities}})
                    for (std::size_t k=0;k<values.size();++k) diagnostics << row << '\t' << factor << '\t' << name << '\t' << k << '\t' << values[k] << '\n';
                if (!result.pressure || result.molar_densities.size()!=2) throw std::runtime_error(result.message);
                const double p=*result.pressure, rl=result.molar_densities[0], rv=result.molar_densities[1];
                const auto l=model.state(t,rl,{1}), v=model.state(t,rv,{1});
                const double chem=std::abs(std::log(rl/rv)+*l.residual_chemical_potential_over_rt[0]-*v.residual_chemical_potential_over_rt[0]);
                const double material=std::abs(result.amounts.at(0)*result.mole_fractions.at(0)+result.amounts.at(1)*result.mole_fractions.at(1)-1);
                const double pc=std::max(std::abs(l.pressure-p),std::abs(v.pressure-p))/p;
                if (factor==.5) first=p;
                const double agreement=std::abs(p/first-1), percent=100*(p-pressure)/pressure;
                const bool numerical=result.success && rl>rv && rv>0 && l.pressure_density_derivative>0
                    && v.pressure_density_derivative>0 && chem<=1e-8 && material<=1e-10 && pc<=1e-8 && agreement<=1e-7;
                vapor << row << '\t' << source << '\t' << t << '\t' << pressure << '\t' << factor << '\t' << result.success
                    << '\t' << p << '\t' << p-pressure << '\t' << percent << '\t' << rl << '\t' << rv << '\t' << l.pressure_density_derivative
                    << '\t' << v.pressure_density_derivative << '\t' << pc << '\t' << chem << '\t' << material << '\t' << agreement
                    << '\t' << numerical << '\t' << (std::abs(percent)<=5) << '\n';
                maxima[source]=std::max(maxima[source],std::abs(percent));
                for (const auto& [name,value] : std::vector<std::pair<std::string,double>>{{"vapor_pressure_relative",pc},
                    {"chemical_RT",chem},{"material_scaled",material},{"start_relative",agreement}})
                    closure_max[name]=std::isfinite(value) ? std::max(closure_max[name],value) : INFINITY;
                if (!numerical) failures << row << '\t' << factor << '\t' << result.message << "; coexistence closure failed\n";
            } catch (const std::exception& e) { failures << row << '\t' << factor << '\t' << e.what() << '\n'; maxima[source]=INFINITY; }
        }
    }
    std::ofstream closures(output+"/closure-maxima.tsv");
    closures << std::setprecision(17) << "quantity\tmaximum\n";
    for (const auto& [name,value] : closure_max) closures << name << '\t' << value << '\n';
    for (const auto& [name,value] : maxima) summary << name << '\t' << value << '\t' << (value<=5) << '\n';
}
