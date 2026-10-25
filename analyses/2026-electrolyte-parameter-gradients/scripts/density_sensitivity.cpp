// Row-3 diagnosis: differentiate the common-pressure Euler expression at fixed
// T, x and reference.
#define main electrolyte_reference_main
#include "electrolyte_reference.cpp"
#undef main
int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  std::ifstream in(argv[1]);
  auto p = epcsaft::build_parameters(read_parameter_document(in));
  std::size_t count, n;
  in >> count;
  std::string rest;
  std::getline(in, rest);
  for (std::size_t i = 0; i < count; ++i)
    std::getline(in, rest);
  in >> n;
  std::vector<epcsaft::ActiveParameter> active;
  for (std::size_t i = 0; i < n; ++i) {
    int family;
    std::size_t width;
    in >> family >> width;
    epcsaft::ActiveParameter a{static_cast<epcsaft::ParameterFamily>(family),
                               {}};
    for (std::size_t j = 0; j < width; ++j) {
      std::size_t index;
      in >> index;
      a.indices.push_back(index);
    }
    active.push_back(a);
  }
  if (!in)
    throw std::runtime_error("incomplete sensitivity inputs");
  epcsaft::HelmholtzModel model(p, active, std::nullopt,
                                epcsaft::policy::SolverOptions{}.resolved({}));
  std::vector<double> water(p.component_count(), 0.0);
  water[0] = 1.0;
  auto ref = epcsaft::solvent_reference(model, 298.15, 1e5, water);
  auto x = epcsaft::formula_composition(water, p, {2, 5}, {1, 1}, 0.1);
  double D = p.molar_mass_kg_per_mol[0] * 2 * 0.1;
  const double rho0 = 55437.903503625,
               dr = std::nextafter(rho0, INFINITY) - rho0;
  auto e =
      model.evaluate(298.15, 1 / rho0, x, epcsaft::DerivativeOrder::Hessian);
  const std::size_t width = e.gradient.size();
  double dL =
      (-e.gradient[1] - 1e5 / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * 298.15)) /
      (rho0 * rho0);
  for (std::size_t i = 0; i < x.size(); ++i)
    if (i != 0 && x[i] > 0)
      dL += x[i] * e.hessian[(i + 2) * width + 1] / (rho0 * rho0);
  dL = dL / x[0] + 1 / rho0;
  std::cerr << std::setprecision(17) << "dphi_drho=" << -dL / D
            << " ulp_rho=" << dr << " propagation=" << -dL / D * dr << '\n';
  // Off-root finite differences verify the smooth Euler expression, not
  // pressure-closed states.
  const auto phi_at = [&](double rho) {
    auto state = model.state(298.15, rho, x);
    const auto &pure = *ref.pure_solvent_state;
    long double c = static_cast<long double>(state.residual_helmholtz) -
                    pure.residual_helmholtz;
    c += 1e5 / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * 298.15) *
         ((static_cast<long double>(pure.molar_density) - rho) /
          (rho * static_cast<long double>(pure.molar_density)));
    for (std::size_t i = 0; i < x.size(); ++i)
      if (i != 0 && x[i] > 0)
        c -= x[i] * (static_cast<long double>(
                         *state.residual_chemical_potential_over_rt[i]) -
                     *pure.residual_chemical_potential_over_rt[0]);
    return -(-std::log1p(D) + static_cast<double>(c / x[0]) +
             std::log1p((rho - pure.molar_density) / pure.molar_density)) /
           D;
  };
  const double fd01 = (phi_at(rho0 + 0.01) - phi_at(rho0 - 0.01)) / 0.02;
  const double fd001 = (phi_at(rho0 + 0.001) - phi_at(rho0 - 0.001)) / 0.002;
  // The smaller step magnifies a 1e-12 value-level roundoff scale to 1e-9 in
  // the difference.
  if (std::max(std::abs(fd01 + dL / D), std::abs(fd001 + dL / D)) > 1e-9)
    throw std::runtime_error("density derivative finite difference disagrees");
  std::cerr << "step_mol_m3=0.01 finite_difference=" << fd01
            << " step_mol_m3=0.001 finite_difference=" << fd001 << '\n';
  std::cout
      << std::setprecision(17)
      << "density,phi,L,helmholtz,pure_helmholtz,mu_Li,mu_Cl,pure_mu,pressure,"
         "dphi_drho,density_ulp,smooth_change,fd_step_0_01,fd_step_0_001\n";
  for (double rho : {55437.90350362499, 55437.903503625}) {
    auto state = model.state_at_pressure(298.15, 1e5, x,
                                         epcsaft::DensityBranch::Liquid, rho);
    auto activity =
        epcsaft::solvent_activity(model, 298.15, 1e5, x, ref, state, D);
    if (!activity.osmotic_coefficient)
      throw std::runtime_error(activity.cause);
    double phi = std::get<double>(*activity.osmotic_coefficient);
    std::cout << rho << ',' << phi << ',' << -D * phi << ','
              << state.residual_helmholtz << ','
              << ref.pure_solvent_state->residual_helmholtz << ','
              << *state.residual_chemical_potential_over_rt[2] << ','
              << *state.residual_chemical_potential_over_rt[5] << ','
              << *ref.pure_solvent_state->residual_chemical_potential_over_rt[0]
              << ',' << state.pressure << ',' << -dL / D << ',' << dr << ','
              << -dL / D * dr << ',' << fd01 << ',' << fd001 << '\n';
  }
}
