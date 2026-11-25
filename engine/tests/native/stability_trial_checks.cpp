#include "properties_support.hpp"
#include "equilibrium_support.hpp"
#include "equilibrium/phase_stability.cpp"

namespace {
using namespace epcsaft;
using namespace epcsaft::equilibrium;

void recovery_support() {
    const HelmholtzModel model(properties_support::neutral_binary());
    PhasePoint methane, ethane;
    methane.point = ethane.point = {0, 1, std::log(0.1)};
    methane.composition = {{1, 0}};
    ethane.composition = {{0, 1}};
    methane.volume = {model.packing_fraction_coefficient(240, methane.composition[0]) / 0.1};
    ethane.volume = {model.packing_fraction_coefficient(240, ethane.composition[0]) / 0.1};
    std::map<std::vector<double>, std::optional<std::string>, ReferenceTaylorKeyLess> decisions;
    decisions.try_emplace(recovery_point_key(methane), "methane branch refused");
    CHECK(decisions.try_emplace(recovery_point_key(ethane)).second);
    CHECK(!decisions.try_emplace(recovery_point_key(methane)).second);
}

void split_point() {
    StabilityTrial a;
    a.phase.support = {0, 1};
    a.donor = 0;
    a.composition = {0.5, 0.5};
    a.quantities["volume_m3_per_mol"] = 0.02;
    CHECK(same_split_trial(a, a));
    auto b = a;
    b.composition[0] = std::nextafter(b.composition[0], 1.0);
    b.composition[1] = std::nextafter(b.composition[1], 0.0);
    CHECK(!same_split_trial(a, b));
    b = a;
    auto& volume = b.quantities.at("volume_m3_per_mol");
    volume = std::nextafter(volume, 1.0);
    CHECK(!same_split_trial(a, b));
}
}
int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.equilibrium.recovery-support", "recovery-support"}, recovery_support},
        {{"native.equilibrium.relaxed-split-point", "relaxed-split-point"}, split_point},
    });
}
