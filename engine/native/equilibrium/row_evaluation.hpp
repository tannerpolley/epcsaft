// Included once, by request_compile.cpp in namespace epcsaft::equilibrium: CompiledSystem coordinate labels, row values, uncertainty, Jacobians and Hessians.
// The layout owners of generic_point name each column: phase compositions, free amounts and volume slots (log packing,
// or log molar density for an ideal-gas model), free solids, critical mode, loading, free T/P, then multipliers.
std::vector<std::string> CompiledSystem::coordinate_labels() const {
    core::Scope scope("coordinate_labels");
    const auto id = [](std::size_t component) { return "component " + std::to_string(component); };  // model order
    std::vector<std::string> labels(unknown_count());
    if (ideal) {
        for (std::size_t j = 0; j < primal_size; ++j) labels[j] = "ln(n[" + id(problem.phases[0].support[j]) + "]/N_feed)";
        return labels;
    }
    const auto shape = layouts(*this);
    for (std::size_t p = 0; p < shape.size(); ++p) {
        const auto& name = problem.phases[p].name;
        for (std::size_t i = 0; i < shape[p].x.size(); ++i)
            labels[shape[p].x[i]] = "ln x[" + name + "," + id(problem.phases[p].support[i]) + "]";
        if (shape[p].amount) labels[*shape[p].amount] = "n[" + name + "] / mol";
        labels[shape[p].volume] = model->parameters().ideal_gas ? "ln rho[" + name + "] / (mol/m3)" : "ln eta[" + name + "]";
    }
    for (std::size_t s = 0; s < solids.size(); ++s)
        if (solids[s].column) labels[*solids[s].column] = "n_solid[" + problem.solids[s].name + "] / mol";
    for (std::size_t k = 0; k < mode_count; ++k) labels[mode_offset() + k] = "u[" + id(k) + "]";
    if (problem.loading) labels[loading_column()] = "loading";
    std::size_t offset = primal_size - free_intensive_count;
    if (problem.temperature.free) labels[offset++] = "T / K";
    if (problem.pressure.free) labels[offset++] = "P / Pa";
    for (std::size_t k = 0; k < rank; ++k) labels[primal_size + k] = "lambda[moiety " + std::to_string(k) + "]";
    for (std::size_t r = 0; r < charge_row_count; ++r) labels[primal_size + rank + r] = "phi[charge row " + std::to_string(r) + "]";
    return labels;
}

std::vector<core::DifferenceCheck> CompiledSystem::jacobian_difference(const std::vector<double>& point, double relative_step,
    double tau, double kappa) const {
    core::Scope scope("jacobian_difference");
    core::admit_difference_parameters(relative_step, tau, kappa);
    if (point.size() != unknown_count()) throw std::invalid_argument("jacobian difference needs every solved-chart coordinate");
    const auto exact = jacobian(point);
    const auto labels = coordinate_labels();
    const std::size_t n = point.size();
    std::vector<core::DifferenceCheck> out;
    for (std::size_t j = 0; j < n; ++j) {
        auto plus = point, minus = point;
        plus[j] += relative_step * std::max(1.0, std::abs(point[j]));
        minus[j] -= relative_step * std::max(1.0, std::abs(point[j]));
        const auto up = evaluate(plus), down = evaluate(minus);
        for (std::size_t i = 0; i < rows.size(); ++i)
            out.push_back(core::difference_check(std::string(row_kind_name(rows[i].kind)) + "[" + std::to_string(rows[i].phase) + ","
                + std::to_string(rows[i].index) + "] / " + labels[j], exact[i * n + j], up.values[i], down.values[i], up.floor, down.floor,
                0.5 * (plus[j] - minus[j]), tau, kappa));
    }
    return out;
}

PhasePoint CompiledSystem::prepare(const std::vector<double>& point, bool state_only,
    std::optional<core::Convergence>* first) const { return prepare_values(point, state_only, first); }

template <typename Base>
PhasePointValues<Base> CompiledSystem::prepare_values(const std::vector<Base>& point, bool state_only,
    std::optional<core::Convergence>* first_numerical_failure) const {
    core::Scope scope("equilibrium_prepare");
    if (point.size() < primal_size) {
        throw std::invalid_argument("equilibrium point is too short");
    }
    if constexpr (std::is_same_v<Base, double>) if (ideal) {
        if (point.size() != primal_size) throw std::invalid_argument("ideal amount coordinate count");
        PhasePointValues<Base> result; result.point = point;
        result.temperature = problem.temperature.value; result.pressure = problem.pressure.value;
        std::vector<Base> amounts(n_feed.size(), 0.0);
        for (std::size_t j = 0; j < point.size(); ++j) {
            const Base amount = feed_total * wide_math::exp(point[j]);
            if (!wide_math::isfinite(amount) || amount <= 0.0) throw std::invalid_argument("invalid ideal amount");
            amounts[problem.phases[0].support[j]] = amount;
        }
        const Base total = sum(amounts);
        if (!wide_math::isfinite(total) || total <= 0.0) throw std::invalid_argument("invalid ideal total amount");
        for (Base& amount : amounts) amount /= total;
        result.composition.push_back(std::move(amounts)); result.amount.push_back(total);
        return result;
    }
    const std::size_t components = model->parameters().component_count();
    PhasePointValues<Base> result;
    result.point = point;
    result.shape = layouts(*this);
    std::size_t offset = primal_size - free_intensive_count;
    result.temperature = problem.temperature.free ? point[offset++] : problem.temperature.value;
    result.pressure = problem.pressure.free ? point[offset++] : problem.pressure.value;
    if constexpr (std::is_same_v<Base, double>) admit_calorics(*this, result.temperature);
    for (std::size_t phase = 0; phase < result.shape.size(); ++phase) {
        core::AssociationObservations phase_observations(phase);
        const auto& layout = result.shape[phase];
        std::vector<Base> values(components, 0.0);
        Base total = 0.0;
        for (std::size_t i = 0; i < layout.x.size(); ++i) {
            const Base value = wide_math::exp(point[layout.x[i]]);
            if (!wide_math::isfinite(value) || value <= 0.0) throw std::invalid_argument("invalid composition");
            values[problem.phases[phase].support[i]] = value;
            total += value;
        }
        if (!(total > 0.0)) throw std::invalid_argument("composition has zero total");
        result.packing.push_back(packing_terms(*model, result.temperature, values));
        result.eta.push_back(wide_math::exp(point[layout.volume]));
        if (!model->parameters().ideal_gas && log_packing_unbounded() && !admitted_packing_fraction(double(result.eta.back())))
            throw std::invalid_argument("phase " + problem.phases[phase].name + " packing fraction outside the fluid domain");
        result.volume.push_back(result.packing.back().coefficient / result.eta.back());
        result.amount.push_back(layout.amount.has_value()
            ? point[*layout.amount] : phase_amount(problem.phases[phase]));
        const auto order = !state_only ? DerivativeOrder::Hessian : problem.temperature.free
            ? DerivativeOrder::PhaseHessian : DerivativeOrder::IsothermalPhaseHessian;
        if constexpr (std::is_same_v<Base, double>) {
        if (state_only && one_component_two_density_eligible(*this)) {
            result.pure.push_back(model->pure_isothermal(result.temperature, result.volume.back(), total));
            result.derivative.push_back(result.pure.back().derivative);
        } else result.derivative.push_back(model->evaluate(result.temperature, result.volume.back(), values, order));
        } else result.derivative.push_back(model->isothermal(result.temperature, result.volume.back(), values));
        if (first_numerical_failure) core::retain_first(*first_numerical_failure, result.derivative.back().first_numerical_failure);
        scope.probe("gradient_floor", double(result.derivative.back().gradient_floor),
                    static_cast<std::uint32_t>(phase));
        result.composition.push_back(std::move(values));
    }
    result.standard = standard_term(
        *this, double(result.temperature), double(result.pressure), 2, problem.pressure.free, problem.temperature.free);
    return result;
}

template PhasePointValues<double> CompiledSystem::prepare_values(const std::vector<double>&, bool, std::optional<core::Convergence>*) const;
template PhasePointValues<WiderScalar> CompiledSystem::prepare_values(const std::vector<WiderScalar>&, bool, std::optional<core::Convergence>*) const;

PointEvaluation CompiledSystem::evaluate(const std::vector<double>& point) const {
    return evaluate_at(prepare(point));
}

PointEvaluation CompiledSystem::evaluate_at(const PhasePoint& at) const { return evaluate_at<double>(at); }

template <typename Base>
PointEvaluationValues<Base> CompiledSystem::evaluate_at(const PhasePointValues<Base>& prepared) const {
    core::Scope scope("equilibrium_evaluate");
    if constexpr (std::is_same_v<Base, double>) if (ideal) return {ideal_rows(*this, prepared).Forward(0, prepared.point), 0.0, 0.0};
    const auto& point = prepared.point;
    const auto balance = [&] { if constexpr (std::is_same_v<Base, double>) return matvec(B, feed_at(point)); else return this->balance; }();
    const double feed_total = [&] { if constexpr (std::is_same_v<Base, double>) return feed_total_at(prepared); else return this->feed_total; }();
    const auto& shape = prepared.shape;
    const std::size_t components = model->parameters().component_count();
    const auto& z = model->parameters().charge_number;
    const auto& x = prepared.composition;
    const auto& amount = prepared.amount;
    const auto& volume = prepared.volume;
    std::vector<Base> lambda(rank, 0.0);
    for (std::size_t row = 0; row < rank && primal_size + row < point.size(); ++row) {
        lambda[row] = point[primal_size + row];
    }
    std::vector<Base> phi(charge_row_count, 0.0);
    for (std::size_t row = 0; row < charge_row_count
         && primal_size + rank + row < point.size(); ++row) {
        phi[row] = point[primal_size + rank + row];
    }
    const Base rt = GAS_CONSTANT_J_PER_MOL_K * prepared.temperature;
    const auto& standard = prepared.standard.derivative[0];
    PointEvaluationValues<Base> result;
    for (const auto& derivative : prepared.derivative) result.floor = std::max(result.floor, derivative.gradient_floor);
    std::vector<std::optional<std::size_t>> phase_phi(shape.size());
    std::size_t phi_index = 0;
    for (const auto& item : rows) {
        if (item.kind == RowKind::Electro) phase_phi[item.phase] = phi_index++;
    }
    std::vector<Base> series;  // every condition row from one zeroth coefficient
    for (const auto& row : rows) {
        Base value = 0.0;
        if (row.kind == RowKind::Chemical) {
            const auto& derivative = prepared.derivative[row.phase];
            const std::size_t species = row.index;
            Base constraint = 0.0;
            for (std::size_t balance_row = 0; balance_row < rank; ++balance_row) {
                constraint += lambda[balance_row] * B[balance_row][species];
            }
            const Base electro = z[species] * (
                phase_phi[row.phase].has_value() ? phi[*phase_phi[row.phase]] : 0.0);
            value = wide_math::log(std::max(x[row.phase][species] / volume[row.phase],
                                      std::numeric_limits<Base>::min()))
                + derivative.gradient[2 + species] + standard[species]
                - constraint - electro;
        } else if (row.kind == RowKind::Pressure) {
            const auto& derivative = prepared.derivative[row.phase];
            value = volume[row.phase] * (
                (model->parameters().ideal_gas ? sum(x[row.phase]) : 1.0) / volume[row.phase] - derivative.gradient[1] - prepared.pressure / rt);
        } else if (row.kind == RowKind::Moiety) {
            for (std::size_t phase = 0; phase < shape.size(); ++phase) {
                for (std::size_t species = 0; species < components; ++species) {
                    value += amount[phase] * B[row.index][species] * x[phase][species];
                }
            }
            for (const auto& solid : solids)
                if (solid.column) value += solid.balance[row.index] * point[*solid.column];
            value = (value - balance[row.index]) / feed_total;
        } else if (row.kind == RowKind::Electro) {
            for (std::size_t species = 0; species < components; ++species) value += z[species] * x[row.phase][species];
        } else if (row.kind == RowKind::SolidAffinity) {
            if constexpr (std::is_same_v<Base, double>) { const auto& solid = solids[row.index];
            value = solid.present ? -solid_affinity(prepared, row.index) : point[*solid.column]; }
        } else if (condition_row(row.kind)) {
            if constexpr (std::is_same_v<Base, double>) { if (series.empty())
                series = row_series(prepared, {point}, std::vector<Base>(control_width()), 0).rows[0];
            value = series[result.values.size()]; }
        } else {
            for (const std::size_t species : problem.phases[row.phase].support) value += x[row.phase][species];
            value -= 1.0;
        }
        result.values.push_back(value);
    }
    if constexpr (std::is_same_v<Base, double>) scope.probe("row_values", result.values);
    for (std::size_t phase = 0; phase < shape.size(); ++phase) {
        const auto& derivative = prepared.derivative[phase];
        Base ideal = ideal_mechanical_helmholtz(x[phase], volume[phase], Base(1));
        result.objective += amount[phase] * (
            derivative.value + ideal + prepared.pressure * volume[phase] / rt
            #ifdef EPCSAFT_SOFTWARE_WIDER
            + std::inner_product(standard.begin(), standard.end(), x[phase].begin(), 0.0,
                [](double sum, Base product) { return static_cast<double>(sum + product); }, std::multiplies<Base>{})
#else
            + std::inner_product(standard.begin(), standard.end(), x[phase].begin(), 0.0)
#endif
        );
    }
    if constexpr (std::is_same_v<Base, double>) for (std::size_t s = 0; s < solids.size(); ++s)
        if (solids[s].column) result.objective += point[*solids[s].column] * solid_gibbs(prepared, s);
    return result;
}

template PointEvaluationValues<WiderScalar> CompiledSystem::evaluate_at(const PhasePointValues<WiderScalar>&) const;

RowUncertainty CompiledSystem::row_uncertainty(const PhasePoint& at) const {
    const auto balance = matvec(B, feed_at(at.point));
    const double feed_total = feed_total_at(at);
    RowUncertainty out;
    out.terms.assign(rows.size(), 0.0); out.inner = out.terms;
    if (ideal) {
        const auto& support = problem.phases[0].support;
        const double shift = *std::max_element(at.point.begin(), at.point.end());
        double total = 0.0;
        for (double y : at.point) total += std::exp(y - shift);
        const double log_total = shift + std::log(total);
        for (std::size_t r = 0; r < problem.reactions.size(); ++r) {
            out.terms[r] = std::abs(ideal_log_k[r]);
            for (std::size_t j = 0; j < support.size(); ++j) {
                const double activity = ideal->activity_basis == "molarity"
                    ? ideal->molarity_log_activity(at.point[j] + std::log(feed_total), support[j]) : at.point[j] - log_total;
                out.terms[r] += std::abs(problem.reactions[r].stoichiometry[support[j]] * activity);
            }
        }
        for (std::size_t k = 0; k < rank; ++k) {
            auto& terms = out.terms[problem.reactions.size() + k];
            terms = std::abs(balance[k] / feed_total);
            for (std::size_t j = 0; j < support.size(); ++j) terms += std::abs(B[k][support[j]] * std::exp(at.point[j]));
        }
        return out;
    }
    const auto components = model->parameters().component_count(), width = 2 + components + model->active_parameters().size();
    const double rt = GAS_CONSTANT_J_PER_MOL_K * at.temperature;
    std::vector<std::vector<double>> gradient_error(at.shape.size(), std::vector<double>(width, 0.0));
    std::vector<double> energy_error(at.shape.size(), 0.0);
    for (std::size_t p = 0; p < at.shape.size(); ++p) {
        std::vector<double> point{at.temperature, at.volume[p]};
        point.insert(point.end(), at.composition[p].begin(), at.composition[p].end());
        const auto active = model->active_parameter_values(); point.insert(point.end(), active.begin(), active.end());
        for (std::size_t j = 0; j < 2 + components; ++j) {
            if (j >= 2 && at.composition[p][j - 2] == 0.0) { gradient_error[p][j] = NAN; continue; }
            std::vector<double> unit(width, 0.0); unit[j] = 1.0;
            const auto error = model->residual_coefficient_errors({point, unit}, at.derivative[p].site_fractions);
            energy_error[p] = error[0]; gradient_error[p][j] = error[1];
        }
    }
    StandardTerm::NumericalError standard_error;
    if (!reference_diagnostics.empty()) standard_term(*this, at.temperature, at.pressure, 1, false, true, {}, &standard_error);
    if (!standard_error.cause.empty()) { out.cause = standard_error.cause; return out; }
    else if (!standard_error.absolute[0].empty()) for (auto& error : gradient_error)
        for (std::size_t j = 0; j < components; ++j) error[2 + j] += standard_error.absolute[0][j];
    out = row_uncertainty_at(at, gradient_error);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i]; const auto p = row.phase, species = row.index;
        double& terms = out.terms[i]; double& inner = out.inner[i];
        if (row.kind == RowKind::CriticalNormalization) {
            terms = 1;
            for (std::size_t j = 0; j < mode_count; ++j) terms += std::pow(at.point[mode_offset() + j], 2);
        } else if (row.kind == RowKind::AzeotropeCondition) {
            const auto local = problem.phases[0].support[species];
            terms = wide_math::abs(wide_math::log(at.composition[0][local])) + wide_math::abs(wide_math::log(at.composition[1][local]));
        } else if (row.kind == RowKind::SaturationPressure) {
            terms = (wide_math::abs(at.pressure) + problem.loading->pressure) / problem.loading->pressure;
        } else if (row.kind == RowKind::EnvelopeCondition) {
            const auto& e = *problem.envelope;
            terms = wide_math::abs(e.target);
            if (e.kind != EnvelopeCoordinate::LogK) terms += wide_math::abs(wide_math::log(e.kind == EnvelopeCoordinate::LogPressure ? at.pressure : at.temperature));
            else terms += wide_math::abs(wide_math::log(at.composition[0][e.component])) + wide_math::abs(wide_math::log(at.composition[1][e.component]));
        } else if (row.kind == RowKind::Enthalpy || row.kind == RowKind::Entropy) {
            const bool enthalpy = row.kind == RowKind::Enthalpy;
            const double denominator = feed_total * GAS_CONSTANT_J_PER_MOL_K * (enthalpy ? at.temperature : 1.0);
            terms = wide_math::abs(enthalpy ? *problem.enthalpy_j : *problem.entropy_j_per_k) / denominator;
            for (std::size_t a = 0; a < at.shape.size(); ++a) {
                const auto ideal = caloric_ideal(*this, a, at.temperature, at.volume[a], at.composition[a], at.temperature);
                const double weight = at.amount[a] / feed_total;
                terms += wide_math::abs(at.amount[a] * (enthalpy ? ideal.enthalpy : ideal.entropy)) / denominator;
                terms += weight * (at.temperature * wide_math::abs(at.derivative[a].gradient[0])
                    + (enthalpy ? wide_math::abs(at.volume[a] * at.derivative[a].gradient[1]) : wide_math::abs(at.derivative[a].value)));
                inner += weight * (at.temperature * gradient_error[a][0]
                    + (enthalpy ? wide_math::abs(at.volume[a]) * gradient_error[a][1] : energy_error[a]));
            }
        } else if (row.kind == RowKind::CriticalMode || row.kind == RowKind::CriticalCubic) {
            const auto& derivative = at.derivative[p];
            std::vector<double> point{at.temperature, at.volume[p]}, w(width, 0.0);
            point.insert(point.end(), at.composition[p].begin(), at.composition[p].end());
            const auto active = model->active_parameter_values(); point.insert(point.end(), active.begin(), active.end());
            for (std::size_t j = 0; j < components; ++j) w[2 + j] = wide_math::sqrt(n_feed[j] / feed_total) * at.point[mode_offset() + j];
            const auto coefficient_error = [&](const std::vector<double>& direction, std::size_t order) {
                std::vector<std::vector<double>> path(order + 1, std::vector<double>(width, 0.0));
                path[0] = point; path[1] = direction;
                return model->residual_coefficient_errors(path, derivative.site_fractions)[order];
            };
            if (row.kind == RowKind::CriticalMode) {
                const double d = wide_math::sqrt(n_feed[species] / feed_total);
                terms = wide_math::abs(d * w[2 + species] / at.composition[p][species]);
                for (std::size_t j = 0; j < components; ++j)
                    terms += wide_math::abs(d * derivative.hessian[(2 + species) * width + 2 + j] * w[2 + j]);
                // Polarization of E2 supplies this requested mixed contraction without a tensor.
                std::vector<double> unit(width, 0.0), sum = w; unit[2 + species] = 1.0; sum[2 + species] += 1.0;
                inner = wide_math::abs(d) * (coefficient_error(sum, 2) + coefficient_error(unit, 2) + coefficient_error(w, 2));
            } else {
                const auto nominal = model->residual_energy_action(at.temperature, at.volume[p], at.composition[p], {w, w, w});
                if (!nominal.action) { out.cause = "critical_contraction_unavailable"; continue; }
                terms = wide_math::abs(*nominal.action);
                for (std::size_t j = 0; j < components; ++j)
                    terms += wide_math::abs(std::pow(w[2 + j], 3) / std::pow(at.composition[p][j], 2));
                inner = 6 * coefficient_error(w, 3);
            }
        } else if (row.kind == RowKind::SolidAffinity) {
            const auto& solid = solids[species]; const auto& spec = problem.solids[species];
            if (!solid.present) { terms = wide_math::abs(at.point[*solid.column]); continue; }
            terms = wide_math::abs(solid.log_k) + wide_math::abs(spec.molar_volume_m3_per_mol * (at.pressure - *spec.reference->reference_pressure_pa) / rt);
            for (std::size_t j = 0; j < components; ++j) terms += wide_math::abs(spec.stoichiometry[j] * at.standard.derivative[0][j]);
            for (std::size_t j = 0; j < rank; ++j) terms += wide_math::abs(at.point[primal_size + j] * solid.balance[j]);
            auto reference = *solid.reference;
            reference_errors(*model, at.temperature, spec.stoichiometry, reference);
            if (!reference.numerical_error_cause.empty() || !standard_error.cause.empty()) {
                inner = NAN; out.cause = !reference.numerical_error_cause.empty() ? reference.numerical_error_cause : standard_error.cause;
            } else {
                const double source = spec.correlation ? evaluate_reaction_properties(*spec.correlation, at.temperature).ln_k : wide_math::log(*spec.equilibrium_constant);
                const double scale = spec.reference->source_basis == ReferenceSourceBasis::CommonMolalityInfiniteDilution
                    ? reaction_reference_log_scale(*model, *problem.source_reference, spec.stoichiometry) : 0.0;
                inner = reference.numerical_error[0] + std::numeric_limits<double>::epsilon()
                    * (wide_math::abs(source) + wide_math::abs(scale) + wide_math::abs(*reference.terminal_contraction));
                for (std::size_t j = 0; j < components; ++j) if (!standard_error.absolute[0].empty())
                    inner += wide_math::abs(spec.stoichiometry[j]) * standard_error.absolute[0][j];
            }

        }
    }
    return out;
}

template <typename Base>
RowUncertaintyValues<Base> CompiledSystem::row_uncertainty_at(const PhasePointValues<Base>& at,
    const std::vector<std::vector<Base>>& gradient_error) const {
    const auto balance = [&] { if constexpr (std::is_same_v<Base, double>) return matvec(B, feed_at(at.point)); else return this->balance; }();
    const double feed_total = [&] { if constexpr (std::is_same_v<Base, double>) return feed_total_at(at); else return this->feed_total; }();
    RowUncertaintyValues<Base> out;
    out.terms.assign(rows.size(), 0); out.inner = out.terms;
    const auto components = model->parameters().component_count();
    const Base rt = GAS_CONSTANT_J_PER_MOL_K * at.temperature;
    std::vector<std::optional<std::size_t>> phi(at.shape.size());
    std::size_t charge = 0;
    for (const auto& row : rows) if (row.kind == RowKind::Electro) phi[row.phase] = primal_size + rank + charge++;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i]; const auto p = row.phase, species = row.index;
        Base& terms = out.terms[i]; Base& inner = out.inner[i];
        if (row.kind == RowKind::Chemical) {
            terms = wide_math::abs(wide_math::log(at.composition[p][species] / at.volume[p]))
                + wide_math::abs(at.derivative[p].gradient[2 + species]) + wide_math::abs(at.standard.derivative[0][species]);
            inner = gradient_error[p][2 + species];
            for (std::size_t k = 0; k < rank; ++k) terms += wide_math::abs(at.point[primal_size + k] * B[k][species]);
            if (phi[p]) terms += wide_math::abs(model->parameters().charge_number[species] * at.point[*phi[p]]);
        } else if (row.kind == RowKind::Pressure) {
            terms = (model->parameters().ideal_gas ? wide_math::abs(sum(at.composition[p])) : 1.0)
                + wide_math::abs(at.volume[p] * at.derivative[p].gradient[1]) + wide_math::abs(at.pressure * at.volume[p] / rt);
            inner = wide_math::abs(at.volume[p]) * gradient_error[p][1];
        } else if (row.kind == RowKind::Moiety) {
            terms = wide_math::abs(balance[species]);
            for (std::size_t a = 0; a < at.shape.size(); ++a)
                for (std::size_t j = 0; j < components; ++j) terms += wide_math::abs(at.amount[a] * B[species][j] * at.composition[a][j]);
            for (const auto& solid : solids) if (solid.column) terms += wide_math::abs(solid.balance[species] * at.point[*solid.column]);
            terms /= feed_total;
        } else if (row.kind == RowKind::Electro || row.kind == RowKind::Normalization) {
            terms = row.kind == RowKind::Normalization ? 1.0 : 0.0;
            for (std::size_t j = 0; j < components; ++j)
                terms += wide_math::abs(at.composition[p][j] * (row.kind == RowKind::Electro ? model->parameters().charge_number[j] : 1));
        }
    }
    return out;
}
template RowUncertaintyValues<double> CompiledSystem::row_uncertainty_at(const PhasePoint&, const std::vector<std::vector<double>>&) const;
template RowUncertaintyValues<WiderScalar> CompiledSystem::row_uncertainty_at(const PhasePointValues<WiderScalar>&, const std::vector<std::vector<WiderScalar>>&) const;

double CompiledSystem::solid_gibbs(const PhasePoint& prepared, std::size_t solid) const {
    const auto& spec = problem.solids[solid];
    double value = solids[solid].log_k + spec.molar_volume_m3_per_mol * (prepared.pressure - *spec.reference->reference_pressure_pa)
        / (GAS_CONSTANT_J_PER_MOL_K * prepared.temperature);
    for (std::size_t i = 0; i < spec.stoichiometry.size(); ++i) value += spec.stoichiometry[i] * prepared.standard.derivative[0][i];
    return value;
}

double CompiledSystem::solid_affinity(const PhasePoint& prepared, std::size_t solid) const {
    double value = -solid_gibbs(prepared, solid);
    for (std::size_t j = 0; j < rank && primal_size + j < prepared.point.size(); ++j)
        value += prepared.point[primal_size + j] * solids[solid].balance[j];
    return value;
}

CompiledSystem classify_solids(CompiledSystem compiled, const std::vector<double>& point) {
    for (auto& solid : compiled.solids)
        solid.present = solid.column && *solid.column < point.size() && point[*solid.column] > 0.0;
    return compiled;
}

std::vector<double> CompiledSystem::classified_residuals(
    const PhasePoint& prepared, const std::vector<double>& values
) const {
    std::vector<double> result = trace_moiety_relative(prepared, values);
    const auto numerical = log_charge_rows.empty() ? result : numerical_rows(prepared, values);
    for (const auto& [row, phase] : log_charge_rows)
        result[row] = std::copysign(std::max(wide_math::abs(result[row]), wide_math::abs(numerical[row])), numerical[row]);
    if (feasibility_cell() || values.size() != rows.size()) return result;
    for (std::size_t row = 0; row < rows.size(); ++row) {
        const auto& item = rows[row];
        if (item.kind == RowKind::Chemical) {
            result[row] *= prepared.amount[item.phase] * prepared.composition[item.phase][item.index] / feed_total;
        } else if (item.kind == RowKind::Pressure) {
            result[row] *= prepared.amount[item.phase] / feed_total;
        }
    }
    return result;
}

std::vector<double> CompiledSystem::trace_moiety_relative(const PhasePoint& at, std::vector<double> values, std::size_t width) const {
    return trace_moiety_relative<double>(at, std::move(values), width);
}
template <typename Base>
std::vector<Base> CompiledSystem::trace_moiety_relative(
    const PhasePointValues<Base>& prepared, std::vector<Base> values, std::size_t width
) const {
    const auto n_feed = feed_at(prepared.point);
    const double feed_total = [&] { if constexpr (std::is_same_v<Base, double>) return feed_total_at(prepared); else return this->feed_total; }();
    const Base trace = policy::EQUILIBRIUM_TRACE_MOIETY_FRACTION.value;
    if (values.size() != rows.size() * width) return values;
    for (std::size_t row = 0; row < rows.size(); ++row) {
        if (rows[row].kind != RowKind::Moiety) continue;
        const auto& b = B[rows[row].index];
        Base inventory = 0.0, terms = 0.0;
        for (std::size_t i = 0; i < b.size(); ++i) inventory += wide_math::abs(b[i]) * n_feed[i];
        if (inventory / feed_total > trace) continue;
        for (std::size_t p = 0; p < prepared.amount.size(); ++p)
            for (std::size_t i = 0; i < b.size(); ++i) terms += wide_math::abs(b[i]) * prepared.amount[p] * prepared.composition[p][i];
        for (const auto& solid : solids)
            if (solid.column) terms += wide_math::abs(solid.balance[rows[row].index] * prepared.point[*solid.column]);
        const Base size = std::max(inventory, terms) / feed_total;
        if (!(size > 0.0) || size > trace) continue;
        for (std::size_t j = 0; j < width; ++j) values[row * width + j] /= size;
    }
    return values;
}

template std::vector<WiderScalar> CompiledSystem::trace_moiety_relative(const PhasePointValues<WiderScalar>&, std::vector<WiderScalar>, std::size_t) const;

// Numerical charge representation only; physical equations and Ipopt derivatives remain unchanged.
std::vector<double> CompiledSystem::numerical_rows(
    const PhasePoint& prepared, std::vector<double> values, std::size_t width
) const {
    values = trace_moiety_relative(prepared, std::move(values), width);
    for (const auto& [row, phase] : log_charge_rows) {
        const auto& z = model->parameters().charge_number;
        const auto& shape = prepared.shape[phase];
        const auto& support = problem.phases[phase].support;
        double logarithm[2];
        for (int sign = 0; sign < 2; ++sign) {
            double largest = -INFINITY, sum = 0.0;
            for (std::size_t j = 0; j < support.size(); ++j)
                if (z[support[j]] && (z[support[j]] > 0) == (sign == 0))
                    largest = std::max(largest, prepared.point[shape.x[j]] + std::log(std::abs(z[support[j]])));
            for (std::size_t j = 0; j < support.size(); ++j)
                if (z[support[j]] && (z[support[j]] > 0) == (sign == 0))
                    sum += std::exp(prepared.point[shape.x[j]] + std::log(std::abs(z[support[j]])) - largest);
            logarithm[sign] = largest + std::log(sum);
        }
        if (width == 1) values[row] = logarithm[0] - logarithm[1];
        else {
            std::fill_n(values.begin() + row * width, width, 0.0);
            for (std::size_t j = 0; j < support.size(); ++j) {
                const int charge = z[support[j]];
                if (charge) values[row * width + shape.x[j]] = std::copysign(std::exp(
                    prepared.point[shape.x[j]] + std::log(std::abs(charge)) - logarithm[charge < 0]), charge);
            }
        }
    }
    return values;
}

double CompiledSystem::classified_residual_norm(const std::vector<double>& point) const {
    const auto prepared = prepare(point);
    return maximum_abs(classified_residuals(prepared, evaluate_at(prepared).values));
}

std::vector<double> CompiledSystem::jacobian(const std::vector<double>& point) const {
    return jacobian_at(prepare(point));
}

std::vector<double> CompiledSystem::jacobian_at(const PhasePoint& at) const { return jacobian_at<double>(at); }

template <typename Base>
std::vector<Base> CompiledSystem::jacobian_at(const PhasePointValues<Base>& prepared) const {
    core::Scope scope("equilibrium_jacobian");
    if constexpr (std::is_same_v<Base, double>) {
    if (problem.loading) {
        const auto n = prepared.point.size();
        std::vector<double> out(rows.size() * n), unit(n);
        for (std::size_t j = 0; j < n; ++j) {
            unit[j] = 1;
            const auto series = row_series(prepared, {prepared.point, unit}, std::vector<double>(control_width()), 1);
            if (series.energy != PropertyActionStatus::Available) throw std::runtime_error("affine row action unavailable");
            for (std::size_t i = 0; i < rows.size(); ++i) out[i * n + j] = series.rows[1][i];
            unit[j] = 0;
        }
        return out;
    }
    if (ideal) return ideal_rows(*this, prepared).Jacobian(prepared.point);
    }
    const auto& point = prepared.point;
    const std::size_t n = point.size();
    const auto& compositions = prepared.composition;
    const auto& amount = prepared.amount;
    const auto& volumes = prepared.volume;
    const auto& packings = prepared.packing;
    const auto& derivatives = prepared.derivative;
    const auto& shape = prepared.shape;
    std::size_t offset = primal_size - free_intensive_count;
    const std::optional<std::size_t> temperature = problem.temperature.free
        ? std::optional<std::size_t>(offset++) : std::nullopt;
    const std::optional<std::size_t> pressure = problem.pressure.free
        ? std::optional<std::size_t>(offset++) : std::nullopt;
    const Base T = prepared.temperature;
    const Base P = prepared.pressure;
    const Base rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto& standard_temperature_derivative = prepared.standard.derivative[1];
    const auto& z = model->parameters().charge_number;
    std::vector<std::optional<std::size_t>> phi(problem.phases.size());
    std::size_t phi_offset = 0;
    for (const auto& row : rows)
        if (row.kind == RowKind::Electro) phi[row.phase] = rank + phi_offset++;
    std::vector<Base> result(rows.size() * n, 0.0);
    auto at = [&](std::size_t row, std::size_t column) -> Base& {
        return result[row * n + column];
    };
    std::vector<std::vector<Base>> columns;  // one first coefficient per coordinate serves every condition row
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const auto& row = rows[row_index];
        if (row.kind == RowKind::Chemical) {
            const auto& derivative = derivatives[row.phase];
            const std::size_t width = derivative.gradient.size();
            const auto& phase = problem.phases[row.phase];
            const std::size_t species = row.index;
            const Base physical_v = derivative.hessian[(2 + species) * width + 1]
                - 1.0 / volumes[row.phase];
            for (std::size_t local = 0; local < phase.support.size(); ++local) {
                const std::size_t component = phase.support[local];
                const Base physical_x = derivative.hessian[
                    (2 + species) * width + (2 + component)]
                    + (component == species ? 1.0 / compositions[row.phase][species] : 0.0);
                const Base dv = packings[row.phase].component[component] / prepared.eta[row.phase];
                at(row_index, shape[row.phase].x[local]) =
                    physical_x * compositions[row.phase][component] + physical_v * dv;
            }
            at(row_index, shape[row.phase].volume) = -volumes[row.phase] * physical_v;
            if (temperature.has_value()) at(row_index, *temperature) =
                derivative.hessian[(2 + species) * width]
                + physical_v * packings[row.phase].temperature_derivative / prepared.eta[row.phase]
                + standard_temperature_derivative[species];
            if (pressure.has_value()) at(row_index, *pressure) =
                prepared.standard.pressure_derivative[species];
            for (std::size_t balance_row = 0; balance_row < rank; ++balance_row)
                at(row_index, primal_size + balance_row) = -B[balance_row][species];
            if (phi[row.phase].has_value()) at(row_index, primal_size + *phi[row.phase]) = -z[species];
        } else if (row.kind == RowKind::Pressure) {
            const auto& derivative = derivatives[row.phase];
            const std::size_t width = derivative.gradient.size();
            const auto& phase = problem.phases[row.phase];
            const Base physical_v = -derivative.gradient[1]
                - volumes[row.phase] * derivative.hessian[width + 1] - P / rt;
            for (std::size_t local = 0; local < phase.support.size(); ++local) {
                const std::size_t component = phase.support[local];
                const Base dv = packings[row.phase].component[component] / prepared.eta[row.phase];
                at(row_index, shape[row.phase].x[local]) =
                    -volumes[row.phase] * derivative.hessian[width + 2 + component]
                        * compositions[row.phase][component] + physical_v * dv
                        + (model->parameters().ideal_gas ? compositions[row.phase][component] : 0.0);
            }
            at(row_index, shape[row.phase].volume) = -volumes[row.phase] * physical_v;
            if (temperature.has_value()) at(row_index, *temperature) = -derivative.hessian[width]
                * volumes[row.phase] + P * volumes[row.phase]
                    / (GAS_CONSTANT_J_PER_MOL_K * T * T)
                + physical_v * packings[row.phase].temperature_derivative / prepared.eta[row.phase];
            if (pressure.has_value()) at(row_index, *pressure) = -volumes[row.phase] / rt;
        } else if (row.kind == RowKind::Moiety) {
            for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
                for (std::size_t local = 0; local < problem.phases[phase].support.size(); ++local) {
                    const std::size_t species = problem.phases[phase].support[local];
                    at(row_index, shape[phase].x[local]) = amount[phase] * B[row.index][species]
                        * compositions[phase][species] / feed_total;
                }
                if (shape[phase].amount.has_value()) {
                    for (const std::size_t species : problem.phases[phase].support)
                        at(row_index, *shape[phase].amount) += B[row.index][species]
                            * compositions[phase][species] / feed_total;
                }
            }
            for (const auto& solid : solids)
                if (solid.column) at(row_index, *solid.column) = solid.balance[row.index] / feed_total;
        } else if (row.kind == RowKind::SolidAffinity) {
            const auto& solid = solids[row.index];
            if (!solid.present) at(row_index, *solid.column) = 1.0;
            for (std::size_t j = 0; solid.present && primal_size + j < n; ++j) at(row_index, primal_size + j) = -solid.balance[j];
        } else if (row.kind == RowKind::Electro) {
            for (std::size_t local = 0; local < problem.phases[row.phase].support.size(); ++local)
                at(row_index, shape[row.phase].x[local]) =
                    z[problem.phases[row.phase].support[local]]
                    * compositions[row.phase][problem.phases[row.phase].support[local]];
        } else if (condition_row(row.kind)) {
            if constexpr (std::is_same_v<Base, double>) { std::vector<std::vector<Base>> state = {point, std::vector<Base>(n, 0.0)};
            const std::vector<Base> control(control_width(), 0.0);
            for (std::size_t c = columns.size(); c < n; ++c) {
                state[1][c] = 1.0;
                columns.push_back(row_series(prepared, state, control, 1).rows[1]);
                state[1][c] = 0.0;
            }
            for (std::size_t c = 0; c < n; ++c) at(row_index, c) = columns[c][row_index]; }
        } else {
            for (std::size_t local = 0; local < problem.phases[row.phase].support.size(); ++local)
                at(row_index, shape[row.phase].x[local]) =
                    compositions[row.phase][problem.phases[row.phase].support[local]];
        }
    }
    scope.probe("rows", static_cast<double>(rows.size()));
    scope.probe("columns", static_cast<double>(n));
    return result;
}

template std::vector<WiderScalar> CompiledSystem::jacobian_at(const PhasePointValues<WiderScalar>&) const;

std::vector<double> CompiledSystem::lever_jacobian(const std::vector<double>& point) const {
    core::Scope scope("equilibrium_lever_jacobian");
    if (!reference_diagnostics.empty()) {
        throw std::invalid_argument(
            "reference_unavailable: source reference T/P/parameter lever derivatives are not admitted");
    }
    auto prepared = prepare(point);
    prepared.point.resize(std::max(prepared.point.size(), unknown_count()), 0.0);
    const std::size_t width = control_width();
    std::vector<double> result(rows.size() * width), control(width);
    std::vector<std::vector<double>> state = {prepared.point, std::vector<double>(prepared.point.size())};
    for (std::size_t column = 0; column < width; ++column) {
        control[column] = 1.0;
        // The callable retains free T/P columns: drive their state coordinate rather than the ignored control.
        if (column == 0 && problem.temperature.free)
            state[1][primal_size - free_intensive_count] = 1.0;
        if (column == 1 && problem.pressure.free)
            state[1][primal_size - 1] = 1.0;
        const auto series = row_series(prepared, state, control, 1);
        for (std::size_t row = 0; row < rows.size(); ++row)
            result[row * width + column] = series.rows[1][row];
        control[column] = 0.0;
        std::fill(state[1].begin(), state[1].end(), 0.0);
    }
    return result;
}

std::vector<double> CompiledSystem::objective_gradient(const std::vector<double>& point) const {
    return objective_gradient_at(prepare(point));
}

std::vector<double> CompiledSystem::objective_gradient_at(const PhasePoint& prepared) const {
    core::Scope scope("equilibrium_objective_gradient");
    const auto& point = prepared.point;
    std::vector<double> result(point.size(), 0.0);
    if (feasibility_cell()) return result;
    const std::size_t components = model->parameters().component_count();
    const double T = prepared.temperature;
    const double P = prepared.pressure;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto& standard_value = prepared.standard.derivative[0];
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        const auto& spec = problem.phases[phase];
        const auto& layout = prepared.shape[phase];
        const auto& x = prepared.composition[phase];
        const double amount = prepared.amount[phase];
        const auto& packing = prepared.packing[phase];
        const double eta = prepared.eta[phase];
        const double volume = prepared.volume[phase];
        const auto& derivative = prepared.derivative[phase];
        double standard = 0.0;
        double ideal = ideal_mechanical_helmholtz(x, volume, 1.0, false);
        for (std::size_t species = 0; species < components; ++species) {
            standard += standard_value[species] * x[species];
        }
        const double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
        const double bracket = derivative.value + ideal + P * volume / rt + standard;
        const double physical_v = derivative.gradient[1] - sum_x / volume + P / rt;
        for (std::size_t local = 0; local < spec.support.size(); ++local) {
            const std::size_t species = spec.support[local];
            const double chemical =
                derivative.gradient[2 + species] + std::log(std::max(
                    x[species] / volume, std::numeric_limits<double>::min()))
                + standard_value[species];
            result[layout.x[local]] = amount * (
                chemical * x[species] + physical_v * packing.component[species] / eta);
        }
        result[layout.volume] = -amount * volume * physical_v;
        result[*layout.amount] = bracket;
    }
    for (std::size_t s = 0; s < solids.size(); ++s)
        if (solids[s].column) result[*solids[s].column] = solid_gibbs(prepared, s);
    return result;
}

std::vector<double> CompiledSystem::hessian(
    const std::vector<double>& point,
    const std::vector<double>& constraint_multipliers
) const {
    return hessian_at(prepare(point), constraint_multipliers);
}

std::vector<double> CompiledSystem::hessian_at(
    const PhasePoint& prepared, const std::vector<double>& constraint_multipliers
) const {
    core::Scope scope("equilibrium_hessian");
    const auto& point = prepared.point;
    const std::size_t n = point.size();
    std::vector<double> result(n * n, 0.0);
    if (feasibility_cell()) return constraint_hessian_at(prepared, constraint_multipliers);
    const double P = prepared.pressure;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * prepared.temperature;
    const auto& standard_value = prepared.standard.derivative[0];
    auto put = [&](std::size_t row, std::size_t column, double value) {
        result[row * n + column] = value;
        result[column * n + row] = value;
    };
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        const auto& spec = problem.phases[phase];
        const auto& layout = prepared.shape[phase];
        const auto& x = prepared.composition[phase];
        const std::size_t volume = layout.volume;
        const double phase_amount_value = prepared.amount[phase];
        const auto& packing = prepared.packing[phase];
        const double eta = prepared.eta[phase];
        const double physical_volume = prepared.volume[phase];
        const auto& derivative = prepared.derivative[phase];
        const std::size_t width = derivative.gradient.size();
        const double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
        // Per-mole chemical potentials and pressure gradient; the extensive ones are N_p times them.
        const double pressure_gradient = derivative.gradient[1] - sum_x / physical_volume + P / rt;
        const double physical_v_gradient = phase_amount_value * pressure_gradient;
        const double physical_vv = phase_amount_value * (
            derivative.hessian[width + 1] + sum_x / (physical_volume * physical_volume));
        std::vector<double> mu(spec.support.size()), h_xv(spec.support.size());
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            const std::size_t si = spec.support[i];
            mu[i] = derivative.gradient[2 + si]
                + std::log(std::max(x[si] / physical_volume, std::numeric_limits<double>::min()))
                + standard_value[si];
            h_xv[i] = phase_amount_value * (derivative.hessian[(2 + si) * width + 1] - 1.0 / physical_volume);
        }
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            const std::size_t si = spec.support[i];
            const double dv_i = packing.component[si] / eta;
            for (std::size_t j = 0; j < spec.support.size(); ++j) {
                const std::size_t sj = spec.support[j];
                const double dv_j = packing.component[sj] / eta;
                // Keep the ideal diagonal analytic and apply residual curvature before trace fractions multiply.
                const double h_xx = phase_amount_value * (x[si]
                    * (x[sj] * derivative.hessian[(2 + si) * width + (2 + sj)]))
                    + (si == sj ? phase_amount_value * x[si] : 0.0);
                put(layout.x[i], layout.x[j],
                    h_xx + x[si] * h_xv[i] * dv_j
                    + dv_i * h_xv[j] * x[sj] + dv_i * physical_vv * dv_j
                    + (si == sj ? phase_amount_value * mu[i] * x[si]
                                      + physical_v_gradient * dv_i : 0.0));
            }
            put(layout.x[i], volume,
                -physical_volume * (x[si] * h_xv[i] + dv_i * physical_vv)
                - physical_v_gradient * dv_i);
        }
        put(volume, volume, physical_volume * physical_volume * physical_vv
            + physical_v_gradient * physical_volume);
        if (layout.amount.has_value()) {
            for (std::size_t i = 0; i < spec.support.size(); ++i)
                put(*layout.amount, layout.x[i], mu[i] * x[spec.support[i]]
                    + pressure_gradient * packing.component[spec.support[i]] / eta);
            put(*layout.amount, volume, -physical_volume * pressure_gradient);
        }
    }
    if (!constraint_multipliers.empty()) {
        const auto weighted = constraint_hessian_at(prepared, constraint_multipliers);
        for (std::size_t index = 0; index < result.size(); ++index) result[index] += weighted[index];
    }
    scope.probe("entries", static_cast<double>(result.size()));
    return result;
}

RowSeries CompiledSystem::row_series(
    const PhasePoint& prepared, const std::vector<std::vector<double>>& state,
    const std::vector<double>& control, std::size_t order
) const {
    core::Scope scope("equilibrium_row_series");
    if (ideal) throw std::invalid_argument("ideal solved-state actions unavailable");
    const std::size_t n = prepared.point.size();
    const std::size_t components = model->parameters().component_count();
    const std::size_t width = control_width(), extended = extended_control_width();
    if (state.size() != order + 1 ||
        (control.size() != width && control.size() != extended - 2 && control.size() != extended) ||
        order > (problem.critical ? 2u : 3u))  // a third critical coefficient needs E6
        throw std::invalid_argument("row series shape");
    for (const auto& item : state)
        if (item.size() != n) throw std::invalid_argument("row series state");
    const bool pressure_jets = problem.pressure.free || control[1] != 0.0;
    const bool temperature_jets = problem.temperature.free || control[0] != 0.0;
    const std::vector<double> active_control(control.begin() + static_cast<std::ptrdiff_t>(2 + components), control.begin() + static_cast<std::ptrdiff_t>(width));
    const bool parameter_jets = std::any_of(active_control.begin(), active_control.end(),
        [](double value) { return value != 0.0; });
    for (std::size_t s = 0; control[0] != 0.0 && s < solids.size(); ++s)
        if (!solids[s].log_k_slope)
            throw std::invalid_argument("reference_unavailable: solid " + problem.solids[s].name + " has no admitted ln K temperature slope");
    if (!reference_diagnostics.empty() &&
        (order > 2 && (parameter_jets || pressure_jets || temperature_jets))) {
        throw std::invalid_argument(
            "reference_unavailable: row series pressure/T/parameter order is unsupported");
    }
    // The prepared term carries the current free-T/P jets. Parameter jets enter the rows only via
    // sum_a ds/dtheta_a (theta_a - theta_a,center): jets along every active parameter serve every direction
    // exactly. If any slope is refused, each direction resolves its own, leaving sibling directions intact.
    const bool reuse = !reference_diagnostics.empty() && control[0] == 0.0 && control[1] == 0.0;
    const auto parameter_cache_ready = [&] {
        return prepared.parameter_standard.has_value()
            && !(order > 1 && pressure_jets && prepared.parameter_standard->pressure_parameter_derivative.empty())
            && !(order > 1 && temperature_jets && prepared.parameter_standard->temperature_parameter_derivative.empty());
    };
    if (reuse && parameter_jets && !parameter_cache_ready()) {
        try {
            prepared.parameter_standard = standard_term(*this, prepared.temperature, prepared.pressure, std::max<std::size_t>(order, 1),
                pressure_jets, temperature_jets, std::vector<double>(active_control.size(), 1.0));
        } catch (const std::invalid_argument&) {}
    }
    const auto& standard = reference_diagnostics.empty()
        ? standard_term(*this, prepared.temperature, prepared.pressure, order, false)
        : !reuse || (parameter_jets && !parameter_cache_ready())
            ? standard_term(*this, prepared.temperature, prepared.pressure, order, pressure_jets, temperature_jets, active_control)
        : parameter_jets ? *prepared.parameter_standard : prepared.standard;
    auto geometry = record_geometry(*this, prepared, standard);
    const std::size_t stride = width + 2 + (has_caloric_condition() ? 3 : 0) + (problem.critical ? 2 * components : 0);
    std::vector<std::vector<double>> output(order + 1);
    for (std::size_t k = 0; k <= order; ++k) {
        auto input = k == 0 ? controls(*this, prepared) : std::vector<double>(n + width + (has_caloric_condition() ? 2 : 0), 0.0);
        std::copy(state[k].begin(), state[k].end(), input.begin());
        if (k == 1) {
            std::copy_n(control.begin(), width, input.begin() + static_cast<std::ptrdiff_t>(n));
            if (has_caloric_condition() && control.size() == extended) {
                input[n + width] = control[control.size() - 2];
                input[n + width + 1] = control.back();
            }
        }
        output[k] = geometry.Forward(k, input);
    }
    RowSeries result;
    result.rows.assign(order + 1, std::vector<double>(rows.size(), 0.0));
    result.temperature.assign(order + 1, 0.0);
    result.pressure.assign(order + 1, 0.0);
    if (problem.loading) for (const auto& coefficient : state) result.loading.push_back(coefficient[loading_column()]);
    result.phase.resize(problem.phases.size());
    const std::size_t row_offset = problem.phases.size() * stride;
    const std::size_t temperature_index = primal_size - free_intensive_count;
    result.solid_amount.assign(solids.size(), std::vector<double>(order + 1));
    result.solid_affinity = result.solid_amount;
    for (std::size_t k = 0; k <= order; ++k) {
        for (std::size_t row = 0; row < rows.size(); ++row) result.rows[k][row] = output[k][row_offset + row];
        for (std::size_t s = 0; s < solids.size(); ++s) {
            result.solid_amount[s][k] = output[k][row_offset + rows.size() + 2 * s];
            result.solid_affinity[s][k] = output[k][row_offset + rows.size() + 2 * s + 1];
        }
        result.temperature[k] = k == 0 ? prepared.temperature
            : problem.temperature.free ? state[k][temperature_index] : k == 1 ? control[0] : 0.0;
        result.pressure[k] = k == 0 ? prepared.pressure
            : problem.pressure.free ? state[k][primal_size - 1] : k == 1 ? control[1] : 0.0;
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
            auto& series = result.phase[phase];
            const auto begin = output[k].begin() + static_cast<std::ptrdiff_t>(phase * stride);
            series.point.emplace_back(begin, begin + static_cast<std::ptrdiff_t>(width));
            series.amount.push_back(begin[static_cast<std::ptrdiff_t>(width)]);
            series.packing.push_back(begin[static_cast<std::ptrdiff_t>(width + 1)]);
        }
    }
    if (order > 0 && control.size() > width && !problem.reactions.empty()) {
        using AD = CppAD::AD<double>;
        std::vector<AD> q = {prepared.temperature, 0.0}, functions;
        CppAD::Independent(q);
        for (const auto &reaction : problem.reactions) {
            const auto r = std::size_t(&reaction - problem.reactions.data());
            const auto start = width + 4 * r;
            const double aa = control[start], bb = control[start + 1], cc = control[start + 2],
                         dd = control[start + 3];
            if (!reaction.correlation && (bb != 0.0 || cc != 0.0 || dd != 0.0))
                throw std::invalid_argument("b/c/d action needs a reaction law");
            const ReactionLogPolynomial law{aa, bb, cc, dd,
                reaction.correlation ? reaction.correlation->reference_temperature : 1.0, 0, 0, {}};
            const AD L = reaction_log_value(law, q[0]);
            const AD LT = reaction_log_value(law, q[0], 1);
            const AD LTT = reaction_log_value(law, q[0], 2);
            functions.push_back(q[1] * L);
            functions.push_back(q[1] * LT);
            functions.push_back(q[1] * LTT);
        }
        CppAD::ADFun<double> source(q, functions);
        for (auto &family : result.reaction_shift)
            family.assign(problem.reactions.size(), std::vector<double>(order + 1, 0.0));
        for (std::size_t k = 0; k <= order; ++k) {
            const auto values =
                source.Forward(k, std::vector<double>{result.temperature[k], k == 1 ? 1.0 : 0.0});
            std::vector<double> targets(problem.reactions.size());
            for (std::size_t r = 0; r < targets.size(); ++r) {
                for (std::size_t j = 0; j < 3; ++j)
                    result.reaction_shift[j][r][k] = values[3 * r + j];
                targets[r] = -values[3 * r];
            }
            const auto species =
                species_standard_slopes(prepared.temperature, components, problem.reactions, targets);
            for (std::size_t row = 0; row < rows.size(); ++row)
                if (rows[row].kind == RowKind::Chemical)
                    result.rows[k][row] += species[rows[row].index];
        }
    }
    // Residual gradient coefficients: the Hessian action of each coefficient plus the model's
    // third and fourth contractions of the lower coefficients (Faa di Bruno through order 3).
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        auto& series = result.phase[phase];
        const auto& derivative = prepared.derivative[phase];
        series.gradient.assign(order + 1, std::vector<double>(width, 0.0));
        series.gradient[0] = derivative.gradient;
        std::vector<std::size_t> needed = {1};
        if (has_caloric_condition()) needed.push_back(0);
        for (const std::size_t species : problem.phases[phase].support) needed.push_back(2 + species);
        const auto contraction = [&](std::vector<std::vector<double>> directions, std::optional<std::size_t> component) {
            if (component) {
                std::vector<double> unit(width, 0.0);
                unit[*component] = 1.0;
                directions.push_back(unit);
            }
            // Central inputs are fixed by the prepared state; preserve raw direction bits and shape.
            std::string key(reinterpret_cast<const char*>(&phase), sizeof(phase));
            for (const auto& direction : directions) {
                const auto size = direction.size();
                key.append(reinterpret_cast<const char*>(&size), sizeof(size));
                key.append(reinterpret_cast<const char*>(direction.data()), size * sizeof(double));
            }
            const auto found = prepared.energy_actions.find(key);
            const bool reused = found != prepared.energy_actions.end();
            const auto& energy = reused ? found->second : prepared.energy_actions.emplace(
                std::move(key), model->residual_energy_action(prepared.temperature,
                    prepared.volume[phase], prepared.composition[phase], directions)).first->second;
            result.maximum_energy_order = std::max(result.maximum_energy_order, directions.size());
            auto& total = result.energy_diagnostics;
            if (!reused) {
                total.forward_sweeps += energy.diagnostics.forward_sweeps;
                total.directional_paths += energy.diagnostics.directional_paths;
                total.site_solves += energy.diagnostics.site_solves;
                total.site_factorizations += energy.diagnostics.site_factorizations;
                total.recordings += energy.diagnostics.recordings;
            }
            total.site_rank = energy.diagnostics.site_rank;
            total.site_pivot_ratio = std::max(total.site_pivot_ratio, energy.diagnostics.site_pivot_ratio);
            total.association_residual = std::max(total.association_residual, energy.diagnostics.association_residual);
            if (energy.status != PropertyActionStatus::Available && result.energy == PropertyActionStatus::Available)
                result.energy = energy.status;
            return energy.action.value_or(std::numeric_limits<double>::quiet_NaN());
        };
        for (std::size_t k = 1; k <= order; ++k) {
            for (std::size_t i = 0; i < width; ++i)
                for (std::size_t j = 0; j < width; ++j)
                    series.gradient[k][i] += derivative.hessian[i * width + j] * series.point[k][j];
            if (k == 2) {
                for (const std::size_t j : needed)
                    series.gradient[2][j] += 0.5 * contraction({series.point[1], series.point[1]}, j);
            } else if (k == 3) {
                for (const std::size_t j : needed) {
                    series.gradient[3][j] += contraction({series.point[1], series.point[2]}, j)
                        + contraction({series.point[1], series.point[1], series.point[1]}, j) / 6.0;
                }
            }
        }
        if (has_caloric_condition()) {
            std::vector<double> energy(order + 1, 0.0);
            energy[0] = derivative.value;
            for (std::size_t k = 1; k <= order; ++k)
                energy[k] = std::inner_product(derivative.gradient.begin(), derivative.gradient.end(), series.point[k].begin(), 0.0);
            for (std::size_t k = 2; k <= order; ++k)
                for (std::size_t i = 0; i < width; ++i)
                    for (std::size_t j = 0; j < width; ++j)
                        energy[k] += derivative.hessian[i * width + j] * series.point[1][i]
                            * series.point[k - 1][j] * (k == 2 ? 0.5 : 1.0);
            if (order == 3)
                energy[3] += contraction({series.point[1], series.point[1], series.point[1]}, std::nullopt) / 6.0;
            for (std::size_t k = 0; k <= order; ++k)
                for (std::size_t j = 0; j <= k; ++j) {
                    const auto offset = phase * stride + width + 2;
                    result.rows[k].back() += output[j][offset] * energy[k - j]
                        + output[j][offset + 1] * series.gradient[k - j][0]
                        + output[j][offset + 2] * series.gradient[k - j][1];
                }
        }
        if (problem.critical) {
            // Residual parts of EqID critical_equations through order 2: d_i [F_nn(p) w]_i (E2-E4) and
            // F_nnn(p)[w,w,w] (E3-E5) along p(t), w(t), d(t); a zero direction contributes exactly zero.
            const std::size_t offset = phase * stride + width + 2, mode = rows.size() - components - 2;
            std::vector<std::vector<double>> d(order + 1, std::vector<double>(components)),
                w(order + 1, std::vector<double>(width, 0.0));
            for (std::size_t k = 0; k <= order; ++k)
                for (std::size_t i = 0; i < components; ++i) {
                    d[k][i] = output[k][offset + 2 * i];
                    w[k][2 + i] = output[k][offset + 2 * i + 1];
                }
            const auto& p = series.point;
            const auto F = [&](const std::vector<std::vector<double>>& v, std::optional<std::size_t> j = std::nullopt) {
                const bool zero = std::any_of(v.begin(), v.end(), [](const std::vector<double>& a) {
                    return std::all_of(a.begin(), a.end(), [](double x) { return x == 0.0; }); });
                return zero ? 0.0 : contraction(v, j);
            };
            for (std::size_t i = 0; i < components; ++i) {
                std::vector<double> g(order + 1, 0.0);
                for (std::size_t k = 0; k <= order; ++k)
                    for (std::size_t j = 0; j < width; ++j) g[k] += derivative.hessian[(2 + i) * width + j] * w[k][j];
                if (order > 0) g[1] += F({p[1], w[0]}, 2 + i);
                if (order > 1) g[2] += F({p[1], w[1]}, 2 + i) + F({p[2], w[0]}, 2 + i) + 0.5 * F({p[1], p[1], w[0]}, 2 + i);
                for (std::size_t k = 0; k <= order; ++k)
                    for (std::size_t j = 0; j <= k; ++j) result.rows[k][mode + i] += d[j][i] * g[k - j];
            }
            std::array<double, 3> c{F({w[0], w[0], w[0]}), 0.0, 0.0};
            if (order > 0) c[1] = 3.0 * F({w[0], w[0], w[1]}) + F({w[0], w[0], w[0], p[1]});
            if (order > 1) c[2] = 3.0 * F({w[0], w[0], w[2]}) + 3.0 * F({w[0], w[1], w[1]}) + 3.0 * F({w[0], w[0], w[1], p[1]})
                + F({w[0], w[0], w[0], p[2]}) + 0.5 * F({w[0], w[0], w[0], p[1], p[1]});
            for (std::size_t k = 0; k <= order; ++k) result.rows[k][mode + components + 1] += c[k];
        }
        for (std::size_t k = 0; k <= order; ++k) {
            for (std::size_t row = 0; row < rows.size(); ++row) {
                if (rows[row].phase != phase) continue;
                if (rows[row].kind == RowKind::Chemical) {
                    result.rows[k][row] += series.gradient[k][2 + rows[row].index];
                } else if (rows[row].kind == RowKind::Pressure) {
                    for (std::size_t j = 0; j <= k; ++j)
                        result.rows[k][row] -= series.point[j][1] * series.gradient[k - j][1];
                }
            }
        }
    }
    scope.probe("order", static_cast<double>(order));
    scope.probe("energy_order", static_cast<double>(result.maximum_energy_order));
    return result;
}

std::vector<double> CompiledSystem::constraint_hessian_at(
    const PhasePoint& prepared, const std::vector<double>& constraint_multipliers
) const {
    if (ideal) {
        if (!constraint_multipliers.empty() && constraint_multipliers.size() != rows.size())
            throw std::invalid_argument("constraint multiplier count");
        return ideal_rows(*this, prepared).Hessian(prepared.point,
            constraint_multipliers.empty() ? std::vector<double>(rows.size(), 0.0) : constraint_multipliers);
    }
    if (feasibility_cell()) return pinned_hessian(*this, prepared, constraint_multipliers);
    const auto& point = prepared.point;
    const std::size_t n = point.size();
    std::vector<double> result(n * n, 0.0);
    if (constraint_multipliers.empty()) return result;
    std::size_t multiplier = 0;
    for (const auto& row : rows) {
        if (stationarity_row(row.kind)) continue;
        if (multiplier >= constraint_multipliers.size()) {
            throw std::invalid_argument("constraint multiplier count");
        }
        const double weight = constraint_multipliers[multiplier++];
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
            if (row.kind != RowKind::Moiety && row.phase != phase) continue;
            const auto& spec = problem.phases[phase];
            const auto& layout = prepared.shape[phase];
            const std::size_t amount = *layout.amount;
            for (std::size_t local = 0; local < spec.support.size(); ++local) {
                const std::size_t species = spec.support[local];
                const std::size_t coordinate = layout.x[local];
                const double x = prepared.composition[phase][species];
                double coefficient = 1.0;
                if (row.kind == RowKind::Electro)
                    coefficient = model->parameters().charge_number[species];
                if (row.kind == RowKind::Moiety)
                    coefficient = B[row.index][species] / feed_total;
                result[coordinate * n + coordinate] += weight * coefficient * x
                    * (row.kind == RowKind::Moiety ? point[amount] : 1.0);
                if (row.kind == RowKind::Moiety) {
                    result[coordinate * n + amount] += weight * coefficient * x;
                    result[amount * n + coordinate] += weight * coefficient * x;
                }
            }
        }
    }
    if (multiplier != constraint_multipliers.size())
        throw std::invalid_argument("constraint multiplier count");
    return result;
}
