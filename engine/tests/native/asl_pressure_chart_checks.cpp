#include "../../native/bindings/asl_adapter.cpp"
#include "properties_checks_shared.hpp"

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.regression.asl-pressure-chart", "ASLPressureChart"}, [] {
            PropertyContext context;
            context.model = std::make_unique<HelmholtzModel>(properties_support::identical_ethane_pair());
            std::vector<double> point{300.0, 1e5, 1.0, 0.0};
            const std::vector<bool> live(point.size(), false);
            CHECK(std::isfinite(evaluate_property(context, 1, point, live, 0).values.at(0)));
            point.back() = -0.0;
            CHECK(std::signbit(evaluate_property(context, 1, point, live, 0).q.back()));
            CHECK(evaluate_property(context, 1, point, live, 0).hessians.empty());  // value-only states hold no derivative rows
            CHECK(evaluate_property(context, 1, point, std::vector<bool>(point.size(), true), 2).hessians.size() == 4);
        }},
        {{"native.regression.asl-pressure-center-counts", "ASLPressureCenterCounts"}, [] {
            HelmholtzModel model(properties_support::associating_binary());
            std::vector<std::size_t> factorizations;
            model.pressure_actions(300.0, 1e5, {0.5, 0.5}, epcsaft::DensityBranch::Liquid,
                {epcsaft::PropertyObservable::LogFugacityCoefficient}, std::nullopt, nullptr,
                [&](const HelmholtzModel::PressureAction& evaluate) {
                    for (int k = 0; k < 2; ++k)
                        factorizations.push_back(evaluate({{1.0, 0.0, {0.0, 0.0}, {}}}).diagnostics.site_factorizations);
                });
            CHECK((factorizations == std::vector<std::size_t>{1, 0}));  // one shared association setup, charged once
        }},
    });
}
