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
        require(!problem.temperature.free && !problem.pressure.free && problem.feed.kind == FeedKind::Amounts);
        // Two fluid phases with an identity gap, or one liquid beside declared-present (Free) solids with neither.
        const bool pair = problem.phases.size() == 2;
        require(pair ? identity_component < n && std::isfinite(minimum_gap) && minimum_gap > 0.0
                     : problem.phases.size() == 1 && !problem.solids.empty() && !identity_component && minimum_gap == 0.0);
        for (const auto& solid : problem.solids) require(solid.amount.kind == AmountKind::Free);
        require(problem.topology_declared);
        // At least one liquid; a vapor's declared support holds only zero-charge components.
        require(problem.phases[0].kind == PhaseKind::Liquid || (pair && problem.phases[1].kind == PhaseKind::Liquid));
        for (const auto& p : problem.phases) {
            require(p.is_free() && p.support_declared && !p.support.empty());
            for (auto i : p.support)
                require(i < n && (p.kind == PhaseKind::Liquid || (parameters.charge_number.empty() ? 0 : parameters.charge_number[i]) == 0));
        }
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
    template<class T> void optional(const std::optional<T>& v) { vector(v ? std::vector<T>{*v} : std::vector<T>{}); }
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

// The active parameters in model order; their values are the document's.
inline void put_active(AslRecordWriter& w, const std::vector<ActiveParameter>& active) {
    w.put("active_parameters"); w.put(active.size());
    for (const auto& a : active) w.put(static_cast<std::size_t>(a.family)), w.vector(a.indices);
}

// ln K(T) law coefficients and interval in record order; standard_state_id follows them.
template <class Law> auto law_fields(Law& c) {
    return std::array{&c.a, &c.b, &c.c, &c.d, &c.reference_temperature, &c.temperature_min, &c.temperature_max, &c.e, &c.f};
}

// One reaction record: name, nu, its target, then the reference flag, P0 (empty: system pressure) and the
// ReferenceSourceBasis index.
inline void put_reaction(AslRecordWriter& w, const ReactionSpec& r) {
    if (r.equilibrium_constant.has_value() + r.affinity_j_per_mol.has_value() + r.correlation.has_value() != 1)
        throw std::invalid_argument("ASL requires exactly one reaction target: K, affinity or correlation");
    w.put(r.name); w.vector(r.stoichiometry);
    w.put(r.correlation ? "correlation" : r.equilibrium_constant ? "K" : "affinity");
    if (const auto& c = r.correlation) {
        for (const double* v : law_fields(*c)) w.put(*v);
        w.put(c->standard_state_id);
    } else w.put(r.equilibrium_constant ? *r.equilibrium_constant : *r.affinity_j_per_mol);
    w.put(r.reference.has_value());
    if (r.reference) w.optional(r.reference->reference_pressure_pa), w.put(std::size_t(r.reference->source_basis));
}

// Neutral-reference flag; then solvent index (empty: none), solute molality ratios, m0 and SpeciesReference indices.
inline void put_source_reference(AslRecordWriter& w, const std::optional<NeutralReferenceDeclaration>& s) {
    w.put("source_reference"), w.put(s.has_value());
    if (s) {
        w.optional(s->solvent_index), w.vector(s->solute_molality_ratios), w.put(s->standard_molality_mol_per_kg);
        w.put(s->species_reference.size());
        for (auto v : s->species_reference) w.put(std::size_t(v));
    }
}

// Thermochemistry interval forms, in IdealCorrelation variant order; each writes its coefficients, then its anchors.
inline const std::array<std::string, 5> ideal_tags = {"polynomial", "shomate", "nasa7", "nasa9", "dippr107"};
inline const std::array<decltype(IdealCorrelation::form), 5> ideal_forms = {IdealPolynomial{}, IdealShomate{}, IdealNASA7{}, IdealNASA9{}, IdealDIPPR107{}};
template <class C> std::vector<double*> ideal_anchors(C& c) {
    if constexpr (std::is_same_v<C, IdealShomate>) return {&c.formation_enthalpy};
    if constexpr (std::is_same_v<C, IdealPolynomial> || std::is_same_v<C, IdealDIPPR107>) return {&c.reference_temperature, &c.reference_enthalpy, &c.reference_entropy};
    return {};
}
// The presence flag, then the standard pressure and each component's intervals; the property record's block.
inline void put_thermochemistry(AslRecordWriter& w, const std::optional<ThermochemistryRecord>& t) {
    w.put(t.has_value());
    if (!t) return;
    w.put("standard_pressure"); w.put(t->standard_pressure); w.put("component_intervals"); w.put(t->component_intervals.size());
    for (const auto& intervals : t->component_intervals) {
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
}

// The solid block: the count, the neutral reference where `reference` is given (the property record), then each solid's
// ln K(T) law and V_s. Every solid has a T column, so an isothermal K is declared as a constant law on its interval.
inline void put_solids(AslRecordWriter& w, const std::vector<SolidSpec>& solids,
                       const std::optional<NeutralReferenceDeclaration>* reference = nullptr) {
    w.put("solids");
    w.put(solids.size());
    if (reference && !solids.empty() && !*reference) throw std::invalid_argument("ASL solid rows require a neutral reference");
    if (reference && !solids.empty()) put_source_reference(w, *reference);
    for (const auto& solid : solids) {
        if (!solid.correlation) throw std::invalid_argument("ASL solid rows require a ln K(T) law");
        put_reaction(w, solid_law(solid));
        w.put(solid.molar_volume_m3_per_mol);
    }
}

inline std::string asl_record(const ParameterDocument& d, const std::vector<ActiveParameter>& active,
                              const std::optional<ThermochemistryRecord>& t, const EquilibriumProblem& p,
                              const ContinuationOptions& options, const AslFeedMap& map, const policy::SolverOptions& solver_options) {
    AslRecordWriter w;
    w.put("epcsaft.asl.equilibrium.v8");
    put_document(w, d);
    w.put("temperature"); w.put(p.temperature.value); w.put("pressure"); w.put(p.pressure.value);
    w.put("feed_kind"); w.put("amounts"); w.put("feed_amounts"); w.vector(p.feed.values);
    w.put("phases"); w.put(p.phases.size());
    for (const auto& phase : p.phases) {
        w.put(phase.name); w.put(phase.kind == PhaseKind::Liquid ? "liquid" : "vapor"); w.vector(phase.support); w.put(phase.support_declared);
        w.put("free"); w.put(phase.amount.fraction); w.vector(phase.composition_guess);
        w.optional(phase.packing_guess);
    }
    w.put("reactions"); w.put(p.reactions.size());
    for (const auto& r : p.reactions) put_reaction(w, r);
    w.put("standard_gibbs_j_per_mol"); w.vector(p.standard_gibbs_j_per_mol);
    put_source_reference(w, p.source_reference);
    if (!p.solids.empty()) put_solids(w, p.solids);  // a record without solids keeps its layout
    w.put("topology_declared"); w.put(p.topology_declared); w.put("solver_options");
    const auto effective = p.solver_options.resolved(solver_options);
    for (std::size_t i = 0; i < effective.values.size(); ++i) w.put(effective.value(i));
    put_active(w, active);
    put_thermochemistry(w, t);
    w.put("continuation_options"); w.put(options.initial_step); w.put(options.minimum_step);
    w.put(options.max_steps); w.put(options.max_retries_per_step);
    w.put("feed_offset"); w.vector(map.offset); w.put("feed_matrix"); w.put(map.matrix.size());
    for (const auto& row : map.matrix) w.vector(row);
    w.put("anchor_controls"); w.vector(map.anchor); w.put("phase_identity");
    w.put(map.identity_component); w.put(map.minimum_gap);
    return w.finish();
}

// Property callbacks: the parameter document, the active parameters in model order (their values are the document's),
// a thermochemistry record and the declared solids (each one row of l_s), no state.
inline std::string asl_property_record(const ParameterDocument& d, const std::vector<ActiveParameter>& active,
    const std::optional<ThermochemistryRecord>& t, const policy::SolverOptions& solver_options, bool require_phase_identity,
    const std::vector<SolidSpec>& solids = {}, const std::optional<NeutralReferenceDeclaration>& reference = std::nullopt) {
    AslRecordWriter w;
    w.put("epcsaft.asl.property.v5");
    w.put(require_phase_identity);
    put_document(w, d);
    w.put("solver_options");
    for (std::size_t i = 0; i < solver_options.values.size(); ++i) w.put(solver_options.value(i));
    put_active(w, active);
    put_thermochemistry(w, t);
    put_solids(w, solids, &reference);
    return w.finish();
}
}  // namespace epcsaft::equilibrium
