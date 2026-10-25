#include "properties_checks_shared.hpp"
using namespace properties_checks;
#define main retained_surface_main
#include "surface_example.cpp"
#undef main
int main() {try {
    const HelmholtzModel model(build_parameters(associating_binary()), {}, ideal_record());
    const double T=353.15,P=1e5;const std::vector<double> x{0.25,0.75};
    const double rho=model.state_at_pressure(T,P,x,DensityBranch::Liquid).molar_density;
    const std::vector<PropertyObservable> os{PropertyObservable::ResidualEnthalpy,PropertyObservable::TotalIsobaricHeatCapacity};
    for(bool pressure:{false,true})for(std::size_t order:{1U,2U}) {
        const auto verify=[&](const auto& ds,double mechanical) {
            using D=typename std::decay_t<decltype(ds)>::value_type;
            PropertyActionBatch batch;
            if constexpr(std::is_same_v<D,DensityActionDirection>)batch=model.density_actions(T,mechanical,x,os,ds);
            else batch=model.pressure_actions(T,mechanical,x,DensityBranch::Liquid,os,ds);
            const auto ref=property_action_reference(model.parameters(),{},model.thermochemistry(),T,mechanical,x,os,ds,DensityBranch::Liquid);
            for(std::size_t k=0;k<os.size();++k) {
                const double tolerance=std::array<double,3>{5e-6,5e-5,5e-4}[order-1],allowance=tolerance*ref.scale[k];
                const double action=*batch.results[k].action;
                const auto id=std::string("action/")+(pressure?"pressure/":"density/")+std::to_string(2*(order-1)+k);
                CHECK(recording_support::compare_custom(id.c_str(),"caloric action","J/(mol K)",action,ref.value[k],
                    recording_support::fields({recording_support::field("T_K",T),recording_support::field("P_Pa",P),
                        recording_support::field("rho_mol_per_m3",rho),recording_support::field("order",order),
                        recording_support::field("allowance",allowance),recording_support::field("reference_uncertainty",ref.uncertainty[k])}),
                    "Numerical","unchanged P4 preceding-action ladder","Fresh native surface row against independent P4 reference",
                    "abs(action-reference)<=allowance && reference_uncertainty<=allowance",std::abs(action-ref.value[k])<=allowance && ref.uncertainty[k]<=allowance));
            }
        };
        if(pressure) {std::vector<PressureActionDirection> ds{{1,0,{0,0},{}}};if(order==2)ds.push_back({0,1e5,{0,0},{}});verify(ds,P);}
        else {std::vector<DensityActionDirection> ds{{1,0,{0,0},{}}};if(order==2)ds.push_back({0,1,{0,0},{}});verify(ds,rho);}
    }
    return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
