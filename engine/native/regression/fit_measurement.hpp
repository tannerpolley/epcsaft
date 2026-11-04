// Included once, by fit.cpp inside its anonymous namespace: the measurement-error chart.
// Statistical membership is declared; equal thermodynamic states never create a measurement.
struct MeasurementChart final {
    std::vector<std::size_t> membership, groups, inputs, rows;  // rows: compact matrix order, n + c for inputs
    std::vector<double> measured, sigma, lower, upper, factor;  // factor: row-major L of R = LL^T
    std::vector<std::vector<double>> composition;
    MeasurementChart() = default;

    explicit MeasurementChart(const FitRequest& request) {
        if (!request.uncertainty) return;
        membership.assign(request.observations.size(), request.observations.size());
        const auto& uncertainty = *request.uncertainty;
        require(!request.minimax, "uncertainty with minimax is unsupported");
        require(uncertainty.output_sigma.size() == membership.size(), "missing output sigma");
        for (std::size_t i = 0; i < membership.size(); ++i) {
            const double weight = request.weights.empty() ? 1.0 : request.weights.at(i);
            require(weight == 0.0 || weight == 1.0, "uncertainty requires unit positive relative weights");
            require(request.observations[i].form == ResidualForm::Difference, "physical-unit Gaussian errors require Difference, not LogRatio");
            require(weight == 0.0 || (std::isfinite(uncertainty.output_sigma[i]) && uncertainty.output_sigma[i] > 0.0), "positive finite output sigma");
        }
        const HelmholtzModel model(request.parameters, {}, request.thermochemistry, request.controls.solver_options);
        const auto& charges = request.parameters.charge_number;
        const std::size_t count = request.parameters.component_count();
        for (std::size_t g = 0; g < uncertainty.groups.size(); ++g) {
            const auto& group = uncertainty.groups[g];
            require(!group.rows.empty(), "empty measurement group");
            require(group.rows.front() < membership.size(), "measurement row index");
            const auto& anchor = request.observations.at(group.rows.front());
            const bool ionic = electrolyte_observation(anchor.kind), pressure = pressure_closed(anchor.kind);
            const bool law = anchor.kind == ObservationKind::LogEquilibriumConstant;
            bool objective = false;
            for (const auto row : group.rows) {
                require(row < membership.size() && membership[row] == membership.size(), "overlapping measurement groups");
                membership[row] = g;
                const auto& item = request.observations[row];
                require(!addition(item.kind) && item.kind != ObservationKind::SolvationGibbsEnergy
                    && item.kind != ObservationKind::TransferGibbsEnergy && item.kind != ObservationKind::BulkRelativePermittivity
                    && item.kind != ObservationKind::FeedEnthalpy && item.kind != ObservationKind::SystemEnthalpy
                    && !item.problem.critical && !item.problem.envelope, "uncertain inputs unsupported for this observation kind");
                auto a = anchor, b = item;
                if (law) { a.kind = b.kind = ObservationKind::PhaseDensity; }
                const bool same = pressure ? pressure_closed(item.kind) && item.temperature == anchor.temperature
                    && item.pressure == anchor.pressure && item.composition == anchor.composition && item.branch == anchor.branch
                    : same_state(a, b);
                require(same && (!law || (item.kind == anchor.kind && item.reaction == anchor.reaction)), "conflicting measurement state declarations");
                objective |= request.weights.empty() || request.weights[row] > 0.0;
            }
            std::vector<double> x;
            const bool feed = std::any_of(group.inputs.begin(), group.inputs.end(), [](const auto& input) { return input.kind == MeasuredInputKind::MoleFraction; });
            double total = 1.0;
            if (feed) {
                require(!ionic && !law, "composition uncertainty is unsupported for this observation");
                x = pressure ? anchor.composition : equilibrium::resolve_feed(model, anchor.problem.feed);
                total = std::accumulate(x.begin(), x.end(), 0.0);
                if (pressure) { require(std::abs(total - 1.0) <= policy::COMPOSITION_SUM_TOLERANCE.value, "measured composition sum"); total = 1.0; }
                require(pressure || (std::isfinite(group.feed_total_mol) && group.feed_total_mol > 0.0
                    && std::abs(total - group.feed_total_mol) <= policy::COMPOSITION_SUM_TOLERANCE.value * total), "fixed feed total disagrees with resolved amounts");
                for (auto& value : x) value /= total;
                require(std::abs(std::inner_product(x.begin(), x.end(), charges.begin(), 0.0))
                    <= policy::COMPOSITION_SUM_TOLERANCE.value, "measured feed is not charge neutral");
                require(group.dependent_species.size() == 1 || group.dependent_species.size() == 2, "feed chart needs one or two dependent species");
                for (std::size_t d = 0; d < group.dependent_species.size(); ++d) {
                    const auto species = group.dependent_species[d];
                    require(species < count && x[species] > 0.0, "dependent feed species is absent");
                    require(std::find(group.dependent_species.begin(), group.dependent_species.begin() + d, species)
                        == group.dependent_species.begin() + d, "duplicate dependent feed species");
                }
            } else require(group.dependent_species.empty(), "dependent species without feed coordinates");
            std::vector<std::vector<double>> directions;
            std::vector<double> values;
            for (std::size_t j = 0; j < group.inputs.size(); ++j) {
                const auto& input = group.inputs[j];
                using K = MeasuredInputKind;
                require(input.kind >= K::Temperature && input.kind <= K::Molality, "unknown measured input");
                require(std::isfinite(input.sigma) && input.sigma >= 0.0, "input sigma must be finite and nonnegative");
                require(std::isfinite(input.lower) && std::isfinite(input.upper) && input.lower <= input.upper, "finite physical input bounds");
                for (std::size_t k = 0; k < j; ++k)
                    require(input.kind != group.inputs[k].kind || (input.kind == K::MoleFraction && input.species != group.inputs[k].species), "duplicate measured input target");
                require(input.unit == (input.kind == K::Temperature ? "K" : input.kind == K::Pressure ? "Pa"
                    : input.kind == K::MoleFraction ? "mol/mol" : "mol/kg"), "unsupported measured input unit");
                require(!law || input.kind == K::Temperature, "source ln K admits only uncertain temperature");
                if (input.kind == K::Temperature || input.kind == K::Pressure) {
                    const auto& intensive = input.kind == K::Temperature ? anchor.problem.temperature : anchor.problem.pressure;
                    require(ionic || pressure || (!intensive.free && !intensive.absent), "free or absent intensive is not a measured input");
                    require(input.lower > 0.0, "positive temperature and pressure domain");
                    if (law) for (const auto& reaction : anchor.problem.reactions) if (reaction.name == anchor.reaction && reaction.correlation)
                        require(input.lower >= reaction.correlation->temperature_min && input.upper <= reaction.correlation->temperature_max, "temperature box leaves source-law interval");
                }
                require(input.kind != K::Molality || ionic, "molality uncertainty requires an electrolyte observation");
                require(input.kind != K::Molality || std::all_of(group.rows.begin(), group.rows.end(), [&](const auto row) {
                    return input.lower >= (request.observations[row].kind == ObservationKind::OsmoticCoefficient
                        ? policy::OSMOTIC_MOLALITY_FLOOR.value : 0.0);
                }), "molality box leaves the admitted domain");
                if (ionic) require(std::count(anchor.composition.begin(), anchor.composition.end(), 1.0) == 1
                    && (input.kind == K::Temperature || input.kind == K::Pressure || input.kind == K::Molality), "input uncertainty needs a one-salt declaration in one pure solvent");
                std::vector<double> dx(count, 0.0);
                double value = input.kind == K::Temperature ? (ionic || pressure ? anchor.temperature : anchor.problem.temperature.value)
                    : input.kind == K::Pressure ? (ionic || pressure ? anchor.pressure : anchor.problem.pressure.value) : anchor.molality;
                if (input.kind == K::MoleFraction) {
                    require(input.species < count && x[input.species] > 0.0
                        && std::find(group.dependent_species.begin(), group.dependent_species.end(), input.species) == group.dependent_species.end(), "independent feed species is absent or dependent");
                    value = x[input.species]; dx[input.species] = 1.0;
                    const auto d = group.dependent_species[0];
                    if (group.dependent_species.size() == 1) {
                        require(charges[d] == 0 && charges[input.species] == 0, "charged feed chart needs two distinct dependent charges");
                        dx[d] = -1.0;
                    } else {
                        const auto e = group.dependent_species[1];
                        require(charges[d] != charges[e], "singular normalization/charge block");
                        dx[d] = double(charges[e] - charges[input.species]) / double(charges[d] - charges[e]);
                        dx[e] = -1.0 - dx[d];
                    }
                }
                require(std::isfinite(value) && value >= input.lower && value <= input.upper
                    && (input.sigma == 0.0 || input.lower < input.upper), "measured input outside its physical box");
                directions.push_back(dx); values.push_back(value);
                if (objective && input.sigma > 0.0) {
                    groups.push_back(g); inputs.push_back(j); measured.push_back(value); sigma.push_back(input.sigma);
                    lower.push_back((input.lower - value) / input.sigma); upper.push_back((input.upper - value) / input.sigma);
                    for (auto& v : dx) v *= pressure ? 1.0 : total;
                    composition.push_back(std::move(dx));
                }
            }
            for (std::size_t i = 0; feed && i < count; ++i) {
                double lo = x[i], hi = x[i];
                for (std::size_t j = 0; j < directions.size(); ++j) {
                    const double a = directions[j][i] * (group.inputs[j].lower - values[j]);
                    const double b = directions[j][i] * (group.inputs[j].upper - values[j]);
                    lo += std::min(a, b); hi += std::max(a, b);
                }
                // The feed's unit-sum criterion admits roundoff at exact affine faces.
                require(lo >= -policy::COMPOSITION_SUM_TOLERANCE.value
                    && hi <= 1.0 + policy::COMPOSITION_SUM_TOLERANCE.value, "nonphysical feed box");
            }
        }
        factor_measurement_errors(request);
    }

    void factor_measurement_errors(const FitRequest& request) {
        const auto& uncertainty = *request.uncertainty;
        // One joint matrix over positive-weight outputs, then active inputs; refused rather than repaired (ALG-FIT).
        const auto& matrix = uncertainty.covariance ? uncertainty.covariance : uncertainty.correlation;
        if (!matrix) return;
        require(!uncertainty.covariance || !uncertainty.correlation, "declares both covariance and correlation");
        std::vector<double> deviations;
        for (std::size_t i = 0; i < membership.size(); ++i)
            if (request.weights.empty() || request.weights[i] > 0.0) { rows.push_back(i); deviations.push_back(uncertainty.output_sigma[i]); }
        for (std::size_t c = 0; c < groups.size(); ++c) { rows.push_back(membership.size() + c); deviations.push_back(sigma[c]); }
        const auto m = static_cast<Eigen::Index>(rows.size());
        require(matrix->size() == rows.size() * rows.size(), "measurement matrix dimension");
        require(finite(*matrix), "finite measurement matrix");
        Eigen::MatrixXd r = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(matrix->data(), m, m);
        const Eigen::Map<const Eigen::VectorXd> s(deviations.data(), m);
        if (uncertainty.covariance) r.array() /= (s * s.transpose()).array();  // R = D^-1 C D^-1
        const double allowance = 64.0 * static_cast<double>(m) * std::numeric_limits<double>::epsilon();
        for (Eigen::Index i = 0; i < m; ++i) {
            require(std::abs(r(i, i) - 1.0) <= allowance, "measurement matrix diagonal disagrees with the declared sigmas");
            for (Eigen::Index j = 0; j < i; ++j) require(std::abs(r(i, j) - r(j, i)) <= allowance, "measurement matrix is not symmetric");
        }
        const Eigen::LLT<Eigen::MatrixXd> llt(r);  // reads the lower triangle
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> l = llt.matrixL();
        core::Scope scope("fit_measurement");
        scope.probe("normalized_matrix", r.data(), r.size());
        scope.probe("llt_success", llt.info() == Eigen::Success ? 1.0 : 0.0);
        if (!(llt.info() == Eigen::Success && l.allFinite())) {
            core::Convergence record; record.residual = NAN; record.stop = core::Stop::infeasible;
            record.cause = "fit measurement matrix is not positive definite";
            scope.finish(record); throw core::DeclarationError(record.cause, record);
        }
        factor.assign(l.data(), l.data() + l.size());
    }

    void apply(const FitRequest& request, const std::vector<double>& point, std::vector<Observation>& observations) const {
        const auto width = request.coordinates.coordinates.size();
        for (std::size_t c = 0; c < groups.size(); ++c) {
            const auto& input = request.uncertainty->groups[groups[c]].inputs[inputs[c]];
            const double correction = sigma[c] * point[width + c], value = measured[c] + correction;
            for (const auto row : request.uncertainty->groups[groups[c]].rows) {
                auto& observation = observations[row];
                const bool pressure = pressure_closed(observation.kind), direct = pressure || electrolyte_observation(observation.kind);
                if (input.kind == MeasuredInputKind::Temperature) (direct ? observation.temperature : observation.problem.temperature.value) = value;
                if (input.kind == MeasuredInputKind::Pressure) (direct ? observation.pressure : observation.problem.pressure.value) = value;
                if (input.kind == MeasuredInputKind::Molality) observation.molality = value;
                if (input.kind == MeasuredInputKind::MoleFraction) {
                    if (!pressure && observation.problem.feed.kind != equilibrium::FeedKind::Amounts) {
                        const HelmholtzModel model(request.parameters, {}, request.thermochemistry, request.controls.solver_options);
                        observation.problem.feed.values = equilibrium::resolve_feed(model, observation.problem.feed);
                        observation.problem.feed.kind = equilibrium::FeedKind::Amounts;
                    }
                    auto& x = pressure ? observation.composition : observation.problem.feed.values;
                    for (std::size_t i = 0; i < x.size(); ++i) x[i] += composition[c][i] * correction;
                }
            }
        }
    }

    std::pair<std::vector<ObservationInputDirection>, std::vector<double>> directions(const FitRequest& request, std::size_t row) const {
        std::vector<ObservationInputDirection> result;
        std::vector<double> deviations;
        for (std::size_t c = 0; c < groups.size(); ++c) if (membership[row] == groups[c]) {
            const auto kind = request.uncertainty->groups[groups[c]].inputs[inputs[c]].kind;
            ObservationInputDirection d;
            d.state.feed_amounts_mol = composition[c];
            d.state.active_parameters.assign(request.coordinates.active_parameters().size() + request.holds.size(), 0.0);
            d.state.temperature_k = kind == MeasuredInputKind::Temperature ? 1.0 : 0.0;
            d.state.pressure_pa = kind == MeasuredInputKind::Pressure ? 1.0 : 0.0;
            d.molality = kind == MeasuredInputKind::Molality ? 1.0 : 0.0;
            result.push_back(std::move(d)); deviations.push_back(sigma[c]);
        }
        return {result, deviations};
    }
};
