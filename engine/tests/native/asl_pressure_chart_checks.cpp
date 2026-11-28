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
        }},
    });
}
