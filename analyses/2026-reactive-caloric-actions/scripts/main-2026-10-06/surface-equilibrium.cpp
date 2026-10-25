#include "equilibrium_e5_checks.cpp"
#include "equilibrium_reactive_caloric_checks.cpp"
using namespace equilibrium_checks;
#define main retained_surface_main
#include "surface_example.cpp"
#undef main
int main() {try {
    using namespace reactive_caloric_checks;
    for(bool fixed:{false,true}) {
        auto c=input(true,fixed);c.active.clear();c.problem.reactions[0].reference->source_basis=ReferenceSourceBasis::RawMoleFractionInfiniteDilution;
        const HelmholtzModel model(c.parameters,{},c.records);const auto compiled=compile_problem(model,c.problem);const auto solved=solve_compiled(compiled);
        const auto q=quantities(3);const auto v=values(model,compiled,solved,q);
        const auto exact=solved_state_actions(compiled,solved,{{q[index(q,"total-h")].observable,q[index(q,"total-cp")].observable},{zero(c)}});
        for(std::size_t k=0;k<2;++k) {
            const auto name=k?"total-cp":"total-h";auto id=std::string("reference-calorics/ionic/")+(fixed?"fixed/":"system/")+(k?"cp":"h");
            CHECK(compare(id,q[index(q,name)],*exact.results[k].value,v.values[index(q,name)],"S2 literal RawMoleFraction source; native independent central calorics",1e-10));
        }
        for(bool coupled:{false,true}) {
            e5::Fixture f;f.parameters=c.parameters;f.problem=c.problem;f.scaled=false;
            f.label=std::string(fixed?"fixed":"system")+(coupled?"-T+P":"-T");
            f.directions={e5::direction(1,coupled?1e5:0,{0,0,0})};
            f.observables={e5::obs(SolvedStateObservableKind::PhaseMoleFraction,0,1)};
            CHECK(e5::qualify(f,"native-row-S2",1));
        }
    }
    e5::Fixture f;f.label="surface-thermal";f.parameters=build_parameters(associating_binary());f.record=ideal_record();f.scaled=false;
    f.problem=equilibrium_support::flash(2);f.problem.phases={equilibrium_support::phase("L",PhaseKind::Liquid),equilibrium_support::phase("V",PhaseKind::Vapor)};
    f.problem.temperature=FreeValue::fixed(353.15);f.problem.pressure=FreeValue::fixed(1.8e4);f.problem.feed.values={0.25,0.75};
    f.directions={e5::direction(1,0,{0,0})};f.observables={e5::obs(SolvedStateObservableKind::PhaseComponentAmount,0,1),e5::obs(SolvedStateObservableKind::PhaseLogActivity,0,1),e5::obs(SolvedStateObservableKind::PhaseMoleFraction,0,1),e5::obs(SolvedStateObservableKind::PhaseProperty,0,0,PropertyObservable::TotalEnthalpy),e5::obs(SolvedStateObservableKind::PhaseMolarDensity,0)};
    CHECK(e5::qualify(f,"native-row-surface",1));
    return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
