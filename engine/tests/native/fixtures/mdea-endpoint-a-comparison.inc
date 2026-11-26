// Source declaration comparison:157:problem on the same MDEA Endpoint A.
auto comparison = cap;
comparison.temperature = FreeValue::fixed(293.15);
comparison.pressure = FreeValue::fixed(500000.0);
comparison.feed.values = {1.168, 26.458028962081073, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
PhaseSpec vapor;
vapor.name = "V";
vapor.kind = PhaseKind::Vapor;
vapor.support = {0, 1, 6};
vapor.support_declared = true;
vapor.composition_guess = {0.9, 0.1, 1e-05};
