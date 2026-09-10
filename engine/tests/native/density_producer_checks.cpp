#include "core/numerical_policy.hpp"
#include "recording_support.hpp"

#include <array>
#include <cstddef>

namespace epcsaft::testing {
std::array<bool, 3> density_producer_checks();
}

#define CHECK(condition) recording_support::check( \
    static_cast<bool>(condition), __FILE__, __LINE__, #condition)

namespace {

void producer_boundaries() {
    const auto checks = epcsaft::testing::density_producer_checks();
    const auto cap_iterations = static_cast<std::size_t>(
        epcsaft::policy::DENSITY_MAX_ITERATIONS.value
    ) + 1;
    recording_support::expect_group(
        recording_support::active_test_id(), "density-publication",
        "density-publication", 3, "C3 density producer boundaries", "check"
    );
    recording_support::compare(
        "density-publication", "refine-nonfinite", "1",
        checks[0] ? 1.0 : 0.0, 1.0, {},
        "Relationship", "density-refine",
        "refine rejects nonfinite root before publication", "absolute", 0.0, 0.0,
        checks[0], "native.diagnostics.c3-density-producers", "C3",
        "density-publication", "availability", "residual,floor,amplification",
        "root rejected"
    );
    recording_support::compare(
        "density-publication", "anchor-nonfinite", "1",
        checks[1] ? 1.0 : 0.0, 1.0, {},
        "Relationship", "density-anchor",
        "anchor rejects nonfinite root before publication", "absolute", 0.0, 0.0,
        checks[1], "native.diagnostics.c3-density-producers", "C3",
        "density-publication", "availability", "residual,floor,amplification",
        "root rejected"
    );
    recording_support::compare(
        "density-publication", "finite-cap", "1",
        checks[2] ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("iterations", cap_iterations)}),
        "Relationship", "density-root",
        "finite criterion root remains available at cap", "absolute", 0.0, 0.0,
        checks[2], "native.diagnostics.c3-density-producers", "C3",
        "density-publication", "availability", "residual,iterations",
        "cap root published"
    );
    CHECK(checks[0] && checks[1] && checks[2]);
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.diagnostics.c3-density-producers", "Density_producer_boundaries",
          "Numerical", "Properties", "Explicit", "Numerical", nullptr, "C3",
          "density-publication", "engine/tests/native/density_producer_checks.cpp",
          "density_producer_checks", "C3", "density-publication:3"}, producer_boundaries},
    });
}
