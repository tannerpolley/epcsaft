#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_support.hpp"
#include <iomanip>
#include <iostream>
using namespace epcsaft;
using namespace epcsaft::equilibrium;
// Diagnostic only: inspect the internal corrector; this does not publish an admitted state.
int main(int, char** argv) {
    HelmholtzModel model(properties_support::neutral_binary());
    auto p=equilibrium_support::boundary();
    p.temperature=FreeValue::variable(std::stod(argv[1]));
    p.pressure=FreeValue::fixed(std::stod(argv[2]));
    ReactionSpec r; r.name="iso"; r.stoichiometry={-1.,1.};
    r.correlation=ReactionLogPolynomial{std::stod(argv[3]),std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]),std::stod(argv[7]),std::stod(argv[8]),std::stod(argv[9]),EOS_STANDARD_STATE_ID};
    p.reactions={r};
    // Reproduce the diagnosed self-ratio K=1 start with declared feed compositions.
    for (auto& phase : p.phases) phase.composition_guess={0.5,0.5};
    auto c=compile_problem(model,p);
    auto start=boundary_start(c);
    auto s=solve_equilibrium_from_coordinates(model,p,start);
    std::cout<<std::setprecision(17);
    const auto prepared=c.prepare(start);
    for (size_t phase=0;phase<2;++phase) {
        std::cout<<"SEED,"<<phase<<","<<prepared.temperature<<","<<prepared.pressure<<","<<1./prepared.volume[phase]<<","<<model.packing_fraction_coefficient(prepared.temperature,prepared.composition[phase])/prepared.volume[phase]<<","<<prepared.composition[phase][0]<<","<<prepared.composition[phase][1]<<"\n";
    }
    if (!s.success || !s.temperature || !s.pressure || s.molar_densities.size()!=2
        || s.packing_fractions.size()!=2 || s.mole_fractions.size()!=4) {
        std::cerr<<"State unavailable: "<<s.message<<"\n";
        return 2;
    }
    for (size_t phase=0;phase<2;++phase) {
        std::cout<<"FINAL,"<<phase<<","<<*s.temperature<<","<<*s.pressure<<","<<s.molar_densities[phase]<<","<<s.packing_fractions[phase]<<","<<s.mole_fractions[phase*2]<<","<<s.mole_fractions[phase*2+1]<<"\n";
    }
}
