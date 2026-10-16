#pragma once
#include "equilibrium/gibbs_minimizer.hpp"
#include "surface/parameter_builder.hpp"
#include <algorithm>
#include <array>
#include <sstream>
#include <type_traits>

namespace epcsaft::equilibrium {
// Explicit caller coordinates; no chemical identity selects a feed map or phase.
struct AslFeedMap {
    std::vector<double> offset, anchor;
    std::vector<std::vector<double>> matrix;
    std::size_t identity_component = 0;
    double minimum_gap = 0.0;
    std::vector<double> amounts(const std::vector<double>& q, bool direction = false) const {
        auto result = direction ? std::vector<double>(offset.size(), 0.0) : offset;
        for (std::size_t i = 0; i < result.size(); ++i)
            for (std::size_t j = 0; j < q.size(); ++j) result[i] += matrix[i][j] * q[j];
        return result;
    }
    void validate(const EquilibriumProblem& problem, const Parameters& parameters) const {
        const auto n = parameters.component_count();
        auto require = [](bool ok) { if (!ok) throw std::invalid_argument("invalid ASL affine feed/phase contract"); };
        require(n > 0 && n <= 64 && !anchor.empty() && anchor.size() <= 64);
        require(offset.size() == n && matrix.size() == n && problem.feed.values.size() == n);
        require(identity_component < n && std::isfinite(minimum_gap) && minimum_gap > 0.0);
        require(!problem.temperature.free && !problem.pressure.free && problem.feed.kind == FeedKind::Amounts);
        require(problem.phases.size() == 2 && problem.topology_declared && !problem.source_reference);
        for (const auto& p : problem.phases)
            require(p.kind == PhaseKind::Liquid && p.is_free() && p.support_declared && !p.support.empty());
        for (double q : anchor) require(std::isfinite(q) && q >= 0.0);
        std::vector<double> charge(anchor.size() + 1, 0.0), norm(charge.size(), 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            require(std::isfinite(offset[i]) && offset[i] >= 0.0 && matrix[i].size() == anchor.size());
            for (std::size_t j = 0; j <= anchor.size(); ++j) {
                const double v = j == anchor.size() ? offset[i] : matrix[i][j];
                require(std::isfinite(v));
                const double z = parameters.charge_number.empty() ? 0 : parameters.charge_number[i];
                charge[j] += z * v; norm[j] += std::abs(z * v);
            }
        }
        // Dot-product roundoff only; this does not relax equilibrium or charge tolerances.
        const double roundoff = std::numeric_limits<double>::epsilon() * (n + anchor.size() + 1);
        for (std::size_t j = 0; j < charge.size(); ++j) require(std::abs(charge[j]) <= roundoff * norm[j]);
        const auto feed = amounts(anchor);
        for (std::size_t i = 0; i < n; ++i)
            require(std::isfinite(feed[i]) && feed[i] >= 0.0 && std::isfinite(problem.feed.values[i])
                    && std::abs(feed[i] - problem.feed.values[i]) <= roundoff * std::max(1.0, feed[i]));
    }
};

class AslRecordWriter {
    std::string text_;
public:
    void put(const std::string& value) {
        if (value.size() > 4096 || value.find('\0') != std::string::npos) throw std::invalid_argument("invalid ASL record atom");
        text_ += std::to_string(value.size()) + ":" + value + ",";
    }
    template<class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0> void put(T value) {
        if constexpr (std::is_floating_point_v<T>) {
            if (!std::isfinite(value)) throw std::invalid_argument("nonfinite ASL record value");
            std::ostringstream out; out << std::hexfloat << value; put(out.str());
        } else put(std::to_string(value));
    }
    template<class T> void vector(const std::vector<T>& values) {
        if (values.size() > 64) throw std::invalid_argument("ASL record vector exceeds 64");
        put(values.size()); for (const auto& v : values) put(v);
    }
    std::string finish() const {
        if (text_.size() > 4 * 1024 * 1024) throw std::invalid_argument("ASL record exceeds 4 MiB");
        return text_;
    }
};

inline void put_document(AslRecordWriter& w, const ParameterDocument& d) {
    w.put("component_ids"); w.vector(d.component_ids);
    w.put("site_component"); w.vector(d.site_component);
    w.put("site_multiplicity"); w.vector(d.site_multiplicity);
    w.put("base_choice"); w.put(d.base_choice);
    w.put("association_choice"); w.put(d.association_choice);
    w.put("electrolyte_choice"); w.put(d.electrolyte_choice);
    w.put("ionic_dispersion_choice"); w.put(d.ionic_dispersion_choice);
    w.put("permittivity_choice"); w.put(d.permittivity_choice);
    w.put("polar_choice"); w.put(d.polar_choice);
    if (d.entries.size() > 4096) throw std::invalid_argument("ASL record exceeds 4096 parameters");
    w.put("parameter_entries"); w.put(d.entries.size());
    for (const auto& e : d.entries) { w.put(e.identity); w.put(e.field); w.vector(e.index); w.put(e.magnitude); w.put(e.unit); }
}

inline std::string asl_record(const ParameterDocument& d, const EquilibriumProblem& p,
                              const ContinuationOptions& options, const AslFeedMap& map, const policy::SolverOptions& solver_options = {}) {
    AslRecordWriter w;
    w.put("epcsaft.asl.equilibrium.v3");
    put_document(w, d);
    w.put("temperature"); w.put(p.temperature.value); w.put("pressure"); w.put(p.pressure.value);
    w.put("feed_kind"); w.put("amounts"); w.put("feed_amounts"); w.vector(p.feed.values);
    w.put("phases"); w.put(p.phases.size());
    for (const auto& phase : p.phases) {
        w.put(phase.name); w.put("liquid"); w.vector(phase.support); w.put(phase.support_declared);
        w.put("free"); w.put(phase.amount.fraction); w.vector(phase.composition_guess);
        w.vector(phase.packing_guess ? std::vector<double>{*phase.packing_guess} : std::vector<double>{});
    }
    w.put("reactions"); w.put(p.reactions.size());
    for (const auto& r : p.reactions) {
        if (r.correlation || r.reference || r.equilibrium_constant.has_value() == r.affinity_j_per_mol.has_value())
            throw std::invalid_argument("ASL requires exactly one scalar reaction target");
        w.put(r.name); w.vector(r.stoichiometry); w.put(r.equilibrium_constant ? "K" : "affinity");
        w.put(r.equilibrium_constant ? *r.equilibrium_constant : *r.affinity_j_per_mol);
    }
    w.put("standard_gibbs_j_per_mol"); w.vector(p.standard_gibbs_j_per_mol);
    w.put("topology_declared"); w.put(p.topology_declared); w.put("solver_options");
    const auto effective = p.solver_options.resolved(solver_options);
    for (std::size_t i = 0; i < effective.values.size(); ++i) w.put(effective.value(i));
    w.put("continuation_options"); w.put(options.initial_step); w.put(options.minimum_step);
    w.put(options.max_steps); w.put(options.max_retries_per_step);
    w.put("feed_offset"); w.vector(map.offset); w.put("feed_matrix"); w.put(map.matrix.size());
    for (const auto& row : map.matrix) w.vector(row);
    w.put("anchor_controls"); w.vector(map.anchor); w.put("phase_identity");
    w.put(map.identity_component); w.put(map.minimum_gap);
    return w.finish();
}

// Thermochemistry interval forms, in IdealCorrelation variant order; each writes its coefficients, then its anchors.
inline const std::array<std::string, 5> ideal_tags = {"polynomial", "shomate", "nasa7", "nasa9", "dippr107"};
inline const std::array<decltype(IdealCorrelation::form), 5> ideal_forms = {IdealPolynomial{}, IdealShomate{}, IdealNASA7{}, IdealNASA9{}, IdealDIPPR107{}};
template <class C> std::vector<double*> ideal_anchors(C& c) {
    if constexpr (std::is_same_v<C, IdealShomate>) return {&c.formation_enthalpy};
    if constexpr (std::is_same_v<C, IdealPolynomial> || std::is_same_v<C, IdealDIPPR107>) return {&c.reference_temperature, &c.reference_enthalpy, &c.reference_entropy};
    return {};
}
// Property callbacks: the parameter document and a thermochemistry record, no state.
inline std::string asl_property_record(const ParameterDocument& d, const ThermochemistryRecord& t, const policy::SolverOptions& solver_options = {}) {
    AslRecordWriter w;
    w.put("epcsaft.asl.property.v2");
    put_document(w, d);
    w.put("solver_options");
    for (std::size_t i = 0; i < solver_options.values.size(); ++i) w.put(solver_options.value(i));
    w.put("standard_pressure"); w.put(t.standard_pressure); w.put("component_intervals"); w.put(t.component_intervals.size());
    for (const auto& intervals : t.component_intervals) {
        w.put(intervals.size());
        for (const auto& item : intervals) {
            w.put(item.lower); w.put(item.upper); w.put(item.lower_inclusive); w.put(item.upper_inclusive);
            w.put(item.correlation.source_pressure); w.put(ideal_tags[item.correlation.form.index()]);
            std::visit([&](auto c) {
                w.vector(std::vector<double>(c.coefficients.begin(), c.coefficients.end()));
                w.put(ideal_anchors(c).size());
                for (const double* value : ideal_anchors(c)) w.put(*value);
            }, item.correlation.form);
        }
    }
    return w.finish();
}
}  // namespace epcsaft::equilibrium
