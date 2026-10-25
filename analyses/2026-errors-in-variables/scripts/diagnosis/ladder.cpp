// Native reference-ladder probes of the frozen water/cation/anion fixture.
// T in K, P in Pa, lengths in angstrom, salt molality in mol/kg.
// Values and input/parameter actions use the Engine; no fit is run.
#include "binary128_support_probe.hpp"
#include "properties_support.hpp"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <tuple>
using namespace epcsaft;
namespace fs = std::filesystem;
fs::path destination;
Parameters parameters(double length) {
  auto p = properties_support::aqueous_electrolyte();
  p.sigma_constant[1] = length;
  p.packing_diameter_angstrom[1] = p.debye_huckel_diameter_angstrom[1] =
      0.88 * length;
  p.born_diameter_angstrom[1] = 3.5;
  return p;
}
std::vector<ActiveParameter> active() {
  return {{detail::ParameterFamily::PackingDiameter, {1}},
          {detail::ParameterFamily::SigmaConstant, {1}},
          {detail::ParameterFamily::DebyeHuckelDiameter, {1}}};
}
void ladder(std::string name, double t, double p, double length) {
  HelmholtzModel model(parameters(length), active());
  std::ofstream out(destination / (name + "-ladder.csv"));
  out << std::setprecision(17);
  out << "m,rho,lnphi_plus,lnphi_minus,delta_plus,delta_minus,T_plus,T_minus,P_"
         "plus,P_minus,parameter_plus,parameter_minus,pressure_residual,"
         "iterations,evaluations,microseconds\n";
  std::optional<double> anchor;
  std::array<double, 2> previous{};
  for (double m : {1e-6, 1e-8, 1e-10, 4e-12, 2.25e-12, 1.21e-12, 1.0201e-12,
                   1.002001e-12, 1e-12, 1e-14, 1e-16}) {
    const auto x =
        formula_composition({1, 0, 0}, model.parameters(), {1, 2}, {1, 1}, m);
    const auto start = std::chrono::steady_clock::now();
    const auto state =
        model.state_at_pressure(t, p, x, DensityBranch::Liquid, anchor);
    anchor = state.molar_density;
    const auto a =
        model.pressure_actions(t, p, x, DensityBranch::Liquid,
                               {PropertyObservable::LogFugacityCoefficient},
                               {{1, 0, {0, 0, 0}, {}}}, anchor);
    const auto b =
        model.pressure_actions(t, p, x, DensityBranch::Liquid,
                               {PropertyObservable::LogFugacityCoefficient},
                               {{0, 1, {0, 0, 0}, {}}}, anchor);
    std::array<double, 2> value{}, T{}, P{}, theta{};
    for (std::size_t i = 0; i < 2; ++i) {
      value[i] =
          require_component_property((*state.log_fugacity_coefficient)[i + 1]);
      theta[i] =
          0.88 * require_component_property(
                     state.log_fugacity_parameter_derivative[(i + 1) * 3]) +
          require_component_property(
              state.log_fugacity_parameter_derivative[(i + 1) * 3 + 1]) +
          0.88 * require_component_property(
                     state.log_fugacity_parameter_derivative[(i + 1) * 3 + 2]);
    }
    for (std::size_t k = 0; k < a.results.size(); ++k) {
      const auto &r = a.results[k];
      if (!r.component || *r.component == 0)
        continue;
      if (!r.action)
        throw std::runtime_error("unavailable ladder action");
      T[*r.component - 1] = *r.action;
    }
    for (const auto &r : b.results)
      if (r.component && *r.component != 0) {
        if (!r.action)
          throw std::runtime_error("unavailable ladder pressure action");
        P[*r.component - 1] = *r.action;
      }
    out << m << ',' << state.molar_density << ',' << value[0] << ',' << value[1]
        << ',' << value[0] - previous[0] << ',' << value[1] - previous[1] << ','
        << T[0] << ',' << T[1] << ',' << P[0] << ',' << P[1] << ',' << theta[0]
        << ',' << theta[1] << ',' << a.diagnostics.pressure_residual << ','
        << state.density_convergence->iterations << ','
        << state.density_convergence->evaluations << ','
        << std::chrono::duration<double, std::micro>(
               std::chrono::steady_clock::now() - start)
               .count()
        << '\n';
    previous = value;
  }
  // Identical terminal composition, both original and refined predecessor
  // seeds.
  const auto x =
      formula_composition({1, 0, 0}, model.parameters(), {1, 2}, {1, 1}, 1e-12);
  const auto cold = model.state_at_pressure(t, p, x, DensityBranch::Liquid);
  const auto ref =
      infinite_dilution_reference(model, t, p, {1, 0, 0}, {1, 2}, {1, 1});
  std::ofstream identity(destination / (name + "-terminal.csv"));
  identity << std::setprecision(17);
  identity << "ion,analytic_limit,finite,cold_minus_limit\n";
  for (std::size_t i = 0; i < 2; ++i)
    identity
        << i << ',' << ref.log_fugacity[i] << ','
        << require_component_property((*cold.log_fugacity_coefficient)[i + 1])
        << ','
        << require_component_property((*cold.log_fugacity_coefficient)[i + 1]) -
               ref.log_fugacity[i]
        << '\n';
}
void wide() {
  const auto p = parameters(3.2);
  HelmholtzModel model(p);
  const auto ref = infinite_dilution_reference(model, 298.15, 1e5, {1, 0, 0},
                                               {1, 2}, {1, 1});
  std::ofstream out(destination / "truth-wide.csv");
  out << std::setprecision(17);
  std::ofstream tail(destination / "wide-reference-tail.csv");
  tail << std::setprecision(17);
  tail << "m_ref,gamma_wide,mean_lnphi_difference_from_terminal\n";
  out << "state,m,gamma_double,gamma_wide,gamma_difference,osmotic_difference,"
         "wide_pressure_residual,wide_pressure_limit\n";
  for (int g = 0; g < 16; ++g) {
    const double m = .1 * std::pow(45., g / 15.);
    const auto value = electrolyte_activity(model, 298.15, 1e5, {1, 0, 0},
                                            {1, 2}, {1, 1}, m, ref);
    equilibrium::EquilibriumProblem problem;
    problem.temperature = equilibrium::FreeValue::fixed(298.15);
    problem.pressure = equilibrium::FreeValue::fixed(1e5);
    ThermochemistryRecord records;
    binary128_support::WideIonic oracle(
        p, records, problem,
        formula_composition({1, 0, 0}, p, {1, 2}, {1, 1}, m)[1]);
    const auto w = oracle.activity_values(m, value.state.molar_density);
    if (g == 0)
      for (double mr :
           {1e-10, 1.21e-12, 1.0201e-12, 1.002001e-12, 1e-12, 1e-14, 1e-16}) {
        const auto v = oracle.activity_values(m, value.state.molar_density, mr);
        tail << mr << ',' << double(v[0]) << ','
             << double(CppAD::log(w[0] / v[0])) << '\n';
      }
    out << g << ',' << m << ',' << value.mean_ionic_molality_scale << ','
        << double(w[0]) << ','
        << double(binary128_support::Wide(value.mean_ionic_molality_scale) -
                  w[0])
        << ','
        << double(binary128_support::Wide(
                      std::get<double>(*value.osmotic_coefficient)) -
                  w[1])
        << ',' << oracle.pressure_residual() << ',' << oracle.pressure_limit()
        << '\n';
  }
}
using Table = std::vector<std::vector<std::string>>;
Table read(std::string p) {
  std::ifstream s(p);
  std::string line;
  Table rows;
  std::getline(s, line);
  while (std::getline(s, line)) {
    std::istringstream r(line);
    std::string v;
    std::vector<std::string> row;
    while (std::getline(r, v, ',')) {
      v.erase(std::remove(v.begin(), v.end(), '"'), v.end());
      row.push_back(v);
    }
    rows.push_back(row);
  }
  return rows;
}
void allowances() {
  std::ofstream summary(destination / "allowance-predictions.csv");
  summary << std::setprecision(17)
          << "option,ratio,reference_eiv,reference_fixed,epsilon_conditional\n";
  const auto components = read("analyses/2026-errors-in-variables/results/"
                               "ionic/numerical-components.csv");
  const auto decision =
      read("analyses/2026-errors-in-variables/results/ionic/decision.csv");
  const double e = std::stod(components.at(1).at(2)),
               f = std::stod(components.at(0).at(2));
  const double epsilon = std::stod(decision.at(0).at(2)),
               rest = epsilon - 2 * (e + f);
  const auto ladder = read((destination / "truth-ladder.csv").string());
  std::map<double, std::array<double, 2>> phi;
  for (const auto &r : ladder)
    phi[std::stod(r[0])] = {std::stod(r[2]), std::stod(r[3])};
  auto difference = [&](double a, double b) {
    return (std::abs(phi.at(a)[0] - phi.at(b)[0]) + std::abs(phi.at(a)[1] - phi.at(b)[1])) /
           2;
  };
  const double base = difference(1e-10, 1e-12);
  for (auto [name, a, b] : std::vector<std::tuple<std::string, double, double>>{
           {"unchanged", 1e-10, 1e-12},
           {"fine_1.21e-12", 1.21e-12, 1e-12},
           {"fine_1.0201e-12", 1.0201e-12, 1e-12},
           {"fine_1.002001e-12", 1.002001e-12, 1e-12},
           {"witness_or_new_terminal_1e-14", 1e-12, 1e-14},
           {"new_terminal_1e-16", 1e-14, 1e-16}}) {
    double ratio = difference(a, b) / base;
    summary << name << ',' << ratio << ',' << e * ratio << ',' << f * ratio
            << ',' << rest + 2 * (e + f) * ratio << '\n';
  }
}
int main(int argc, char **argv) {
  if (argc != 2)
    throw std::runtime_error("expected output directory");
  destination = argv[1];
  fs::create_directories(destination);
  ladder("truth", 298.15, 1e5, 3.2);
  ladder("endpoint-0", 298.37981917157344, 97200.481720674681,
         3.245846989133768);
  wide();
  allowances();
}
