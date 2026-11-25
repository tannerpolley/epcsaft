#include "properties_checks_shared.hpp"
#include "properties/density.cpp"

int main() {
    using namespace epcsaft;
    if constexpr (std::numeric_limits<WiderScalar>::digits > std::numeric_limits<double>::digits) {
        using Wide = WiderScalar;
        const HelmholtzModel model(properties_support::aqueous_electrolyte());
        const PressureSeries<Wide, true> a(model, Wide(298.15), {Wide(0.95), Wide(0.025), Wide(0.025)});
        const PressureSeries<Wide, true> b(model, Wide(330), {Wide(0.75), Wide(0.125), Wide(0.125)});
        Wide floor = 0;
        const auto evaluate = [&](const auto& series, Wide density) {
            return series.template residual<2>(density, floor, nullptr, nullptr, {Wide(1), Wide(0), Wide(0), Wide(0)});
        };
        const auto first = evaluate(a, Wide(54'000));
        (void)evaluate(b, Wide(43'000));
        const auto resumed = evaluate(a, Wide(54'000));
        CHECK(first.c == resumed.c);
    }
}
