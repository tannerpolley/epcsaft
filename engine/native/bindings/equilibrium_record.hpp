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
        auto require = [](bool ok, const char* condition) {
            if (!ok) throw std::invalid_argument(std::string("invalid ASL affine feed/phase contract: ") + condition);
        };
        require(n > 0 && n <= 64 && !anchor.empty() && anchor.size() <= 64, "1 to 64 components and 1 to 64 anchor controls");
        require(offset.size() == n && matrix.size() == n && problem.feed.values.size() == n, "feed offset, matrix and amounts have one row per component");
        require(!problem.temperature.free && !problem.pressure.free && problem.feed.kind == FeedKind::Amounts, "fixed T and P with an absolute amount feed");
        // Two fluid phases with an identity gap, or one liquid beside declared-present (Free) solids with neither.
        const bool pair = problem.phases.size() == 2;
        require(pair ? identity_component < n && std::isfinite(minimum_gap) && minimum_gap > 0.0
                     : problem.phases.size() == 1 && !problem.solids.empty() && !identity_component && minimum_gap == 0.0,
                "two phases with an identity component and a positive gap, or one phase beside solids with neither");
        for (const auto& solid : problem.solids) require(solid.amount.kind == AmountKind::Free, "every solid amount free");
        require(problem.topology_declared, "declared topology");
        // At least one liquid; a vapor's declared support holds only zero-charge components.
        require(problem.phases[0].kind == PhaseKind::Liquid || (pair && problem.phases[1].kind == PhaseKind::Liquid), "at least one liquid phase");
        for (const auto& p : problem.phases) {
            require(p.is_free() && p.support_declared && !p.support.empty(), "every phase amount free with a declared nonempty support");
            for (auto i : p.support)
                require(i < n && (p.kind == PhaseKind::Liquid || (parameters.charge_number.empty() ? 0 : parameters.charge_number[i]) == 0),
                        "support indices within the components and only zero-charge components in a vapor");
        }
        for (double q : anchor) require(std::isfinite(q) && q >= 0.0, "finite nonnegative anchor controls");
        std::vector<double> charge(anchor.size() + 1, 0.0), norm(charge.size(), 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            require(std::isfinite(offset[i]) && offset[i] >= 0.0 && matrix[i].size() == anchor.size(),
                    "finite nonnegative feed offset and one matrix column per anchor control");
            for (std::size_t j = 0; j <= anchor.size(); ++j) {
                const double v = j == anchor.size() ? offset[i] : matrix[i][j];
                require(std::isfinite(v), "finite feed matrix");
                const double z = parameters.charge_number.empty() ? 0 : parameters.charge_number[i];
                charge[j] += z * v; norm[j] += std::abs(z * v);
            }
        }
        // Dot-product roundoff only; this does not relax equilibrium or charge tolerances.
        const double roundoff = std::numeric_limits<double>::epsilon() * (n + anchor.size() + 1);
        for (std::size_t j = 0; j < charge.size(); ++j)
            require(std::abs(charge[j]) <= roundoff * norm[j], "electroneutral feed offset and matrix columns");
        const auto feed = amounts(anchor);
        for (std::size_t i = 0; i < n; ++i)
            require(std::isfinite(feed[i]) && feed[i] >= 0.0 && std::isfinite(problem.feed.values[i])
                    && std::abs(feed[i] - problem.feed.values[i]) <= roundoff * std::max(1.0, feed[i]),
                    "feed amounts equal to the map at the anchor controls");
    }
};

// The writing archive; the adapter's RecordReader is the reading one. Both offer tag, text, number, whole and count.
class AslRecordWriter {
    std::string text_;
    void put(const std::string& value) {
        if (value.size() > 4096 || value.find('\0') != std::string::npos) throw std::invalid_argument("invalid ASL record atom");
        text_ += std::to_string(value.size()) + ":" + value + ",";
    }
public:
    void tag(const char* name) { put(name); }
    void text(const std::string& value) { put(value); }
    void number(double value) {
        if (!std::isfinite(value)) throw std::invalid_argument("nonfinite ASL record value");
        std::ostringstream out; out << std::hexfloat << value; put(out.str());
    }
    template<class T> void whole(T value, std::size_t = 0) {
        if constexpr (std::is_enum_v<T>) put(std::to_string(static_cast<std::size_t>(value)));
        else put(std::to_string(value));
    }
    void count(std::size_t n, std::size_t maximum) {
        if (n > maximum) throw std::invalid_argument("ASL record vector exceeds " + std::to_string(maximum));
        whole(n);
    }
    bool present(const char* name, bool present) { if (present) put(name); return present; }  // an optional named block
    std::string finish() const {
        if (text_.size() > 4 * 1024 * 1024) throw std::invalid_argument("ASL record exceeds 4 MiB");
        return text_;
    }
};

// One field walk per record block serves both archives, so the writer and the reader cannot disagree on a field. The
// writer walks a copy of its inputs; checks named for one archive test whether A is the writer.
template <class A, class V, class Each> void list(A& a, V& values, Each each, std::size_t maximum = 64) {
    std::size_t n = values.size();
    a.count(n, maximum);
    if constexpr (std::is_same_v<V, std::vector<typename V::value_type>>) values.resize(n);
    else if (n != values.size()) throw std::invalid_argument("invalid ASL interval");  // a fixed coefficient array
    for (auto& value : values) each(value);
}
template <class A, class T, class Each> void optional(A& a, std::optional<T>& value, Each each) {  // a list of 0 or 1
    std::size_t n = value.has_value();
    a.count(n, 1);
    if (n && !value) value.emplace();
    if (value) each(*value);
}
template <class A> auto numbers(A& a) { return [&a](double& value) { a.number(value); }; }
template <class A> auto wholes(A& a, std::size_t maximum = 64) { return [&a, maximum](auto& value) { a.whole(value, maximum); }; }

template <class A> void document(A& a, ParameterDocument& d) {
    a.tag("component_ids"); list(a, d.component_ids, [&a](std::string& id) { a.text(id); });
    a.tag("site_component"); list(a, d.site_component, wholes(a));
    a.tag("site_multiplicity"); list(a, d.site_multiplicity, wholes(a));
    for (auto [name, choice] : {std::pair{"base_choice", &d.base_choice}, {"association_choice", &d.association_choice},
         {"electrolyte_choice", &d.electrolyte_choice}, {"ionic_dispersion_choice", &d.ionic_dispersion_choice},
         {"permittivity_choice", &d.permittivity_choice}, {"polar_choice", &d.polar_choice}})
        a.tag(name), a.text(*choice);
    // Bound parameter-list allocation and parsing work for serialized records.
    if (d.entries.size() > 4096) throw std::invalid_argument("ASL record exceeds 4096 parameters");
    a.tag("parameter_entries");
    list(a, d.entries, [&a](ParameterEntry& e) {
        a.text(e.identity); a.text(e.field); list(a, e.index, wholes(a)); a.number(e.magnitude); a.text(e.unit);
    }, 4096);
}

// The active parameters in model order; their values are the document's.
template <class A> void active_parameters(A& a, std::vector<ActiveParameter>& active) {
    a.tag("active_parameters");
    list(a, active, [&a](ActiveParameter& p) {
        a.whole(p.family, epcsaft::detail::PARAMETER_FAMILY_COUNT - 1); list(a, p.indices, wholes(a, SIZE_MAX));
    });
}

// Each value through SolverOptions::value, so an unset option travels as its default.
template <class A> void solver_values(A& a, policy::SolverOptions& options) {
    a.tag("solver_options");
    for (std::size_t i = 0; i < options.values.size(); ++i) { double v = options.value(i); a.number(v); options.values[i] = v; }
}

// ln K(T) law coefficients and interval in record order; standard_state_id follows them.
template <class Law> auto law_fields(Law& c) {
    return std::array{&c.a, &c.b, &c.c, &c.d, &c.reference_temperature, &c.temperature_min, &c.temperature_max, &c.e, &c.f};
}

// One reaction record: name, nu, its target, then the reference flag, P0 (empty: system pressure) and the
// ReferenceSourceBasis index.
template <class A> void reaction(A& a, ReactionSpec& r) {
    if (std::is_same_v<A, AslRecordWriter> && r.equilibrium_constant.has_value() + r.affinity_j_per_mol.has_value() + r.correlation.has_value() != 1)
        throw std::invalid_argument("ASL requires exactly one reaction target: K, affinity or correlation");
    a.text(r.name); list(a, r.stoichiometry, numbers(a));
    std::string target = r.correlation ? "correlation" : r.equilibrium_constant ? "K" : "affinity";
    a.text(target);
    if (target == "correlation") {
        auto& c = r.correlation ? *r.correlation : r.correlation.emplace();
        for (double* v : law_fields(c)) a.number(*v);
        a.text(c.standard_state_id);
    } else if (target == "K" || target == "affinity") {
        auto& value = target == "K" ? r.equilibrium_constant : r.affinity_j_per_mol;
        a.number(value ? *value : value.emplace());
    } else throw std::invalid_argument("unsupported reaction target");
    optional(a, r.reference, [&a](ReactionReference& reference) {
        optional(a, reference.reference_pressure_pa, numbers(a)); a.whole(reference.source_basis, 1);
    });
}

// Neutral-reference flag; then solvent index (empty: none), solute molality ratios, m0 and SpeciesReference indices.
template <class A> void source_reference(A& a, std::optional<NeutralReferenceDeclaration>& s) {
    a.tag("source_reference");
    optional(a, s, [&a](NeutralReferenceDeclaration& d) {
        optional(a, d.solvent_index, wholes(a)); list(a, d.solute_molality_ratios, numbers(a));
        a.number(d.standard_molality_mol_per_kg); list(a, d.species_reference, wholes(a, 2));
    });
}

// Thermochemistry interval forms, in IdealCorrelation variant order; each holds its coefficients, then its anchors.
inline const std::array<std::string, 5> ideal_tags = {"polynomial", "shomate", "nasa7", "nasa9", "dippr107"};
inline const std::array<decltype(IdealCorrelation::form), 5> ideal_forms = {IdealPolynomial{}, IdealShomate{}, IdealNASA7{}, IdealNASA9{}, IdealDIPPR107{}};
template <class C> std::vector<double*> ideal_anchors(C& c) {
    if constexpr (std::is_same_v<C, IdealShomate>) return {&c.formation_enthalpy};
    if constexpr (std::is_same_v<C, IdealPolynomial> || std::is_same_v<C, IdealDIPPR107>) return {&c.reference_temperature, &c.reference_enthalpy, &c.reference_entropy};
    return {};
}
// The presence flag, then the standard pressure and each component's intervals.
template <class A> void thermochemistry(A& a, std::optional<ThermochemistryRecord>& t) {
    optional(a, t, [&a](ThermochemistryRecord& record) {
        a.tag("standard_pressure"); a.number(record.standard_pressure); a.tag("component_intervals");
        list(a, record.component_intervals, [&a](std::vector<IdealInterval>& intervals) { list(a, intervals, [&a](IdealInterval& item) {
            a.number(item.lower); a.number(item.upper); a.whole(item.lower_inclusive, 1); a.whole(item.upper_inclusive, 1);
            a.number(item.correlation.source_pressure);
            // NOLINTNEXTLINE(performance-unnecessary-copy-initialization): the reader overwrites tag.
            std::string tag = ideal_tags[item.correlation.form.index()];
            a.text(tag);
            const auto form = std::find(ideal_tags.begin(), ideal_tags.end(), tag) - ideal_tags.begin();
            if (form == 5) throw std::invalid_argument("unknown ASL thermochemistry form " + tag);
            if constexpr (!std::is_same_v<A, AslRecordWriter>) item.correlation.form = ideal_forms[form];
            std::visit([&a](auto& c) {
                list(a, c.coefficients, numbers(a));
                const auto anchors = ideal_anchors(c);
                std::size_t n = anchors.size();
                a.whole(n);
                if (n != anchors.size()) throw std::invalid_argument("invalid ASL interval");
                for (double* value : anchors) a.number(*value);
            }, item.correlation.form);
        }); });
    });
}

// The solid block after its "solids" field: the count, the neutral reference where `reference` is given (the property
// record), then each solid's ln K(T) law and V_s, each solid declared present (Free). Every solid has a T column, so an
// isothermal K is declared as a constant law on its interval.
template <class A> void solid_rows(A& a, std::vector<SolidSpec>& solids, std::optional<NeutralReferenceDeclaration>* reference = nullptr) {
    std::size_t n = solids.size();
    a.whole(n);
    if (reference && n > 0) {
        if (std::is_same_v<A, AslRecordWriter> && !*reference) throw std::invalid_argument("ASL solid rows require a neutral reference");
        source_reference(a, *reference);
    }
    solids.resize(n);
    for (auto& solid : solids) {
        auto r = solid_law(solid);
        if (std::is_same_v<A, AslRecordWriter> && !r.correlation) throw std::invalid_argument("ASL solid rows require a ln K(T) law");
        reaction(a, r);
        if (!r.correlation) throw std::invalid_argument("ASL solid rows require a ln K(T) law");
        solid = {r.name, r.stoichiometry, r.equilibrium_constant, r.correlation, r.reference, solid.molar_volume_m3_per_mol, {}};
        a.number(solid.molar_volume_m3_per_mol);
    }
}

template <class A> void equilibrium_record(A& a, ParameterDocument& d, std::vector<ActiveParameter>& active,
    std::optional<ThermochemistryRecord>& t, EquilibriumProblem& p, ContinuationOptions& options, AslFeedMap& map) {
    a.tag("epcsaft.asl.equilibrium.v8");
    document(a, d);
    a.tag("temperature"); a.number(p.temperature.value); a.tag("pressure"); a.number(p.pressure.value);
    std::string feed = "amounts";
    a.tag("feed_kind"); a.text(feed);
    if (feed != "amounts") throw std::invalid_argument("equilibrium callback requires absolute amounts");
    a.tag("feed_amounts"); list(a, p.feed.values, numbers(a));
    a.tag("phases");
    list(a, p.phases, [&a](PhaseSpec& phase) {
        a.text(phase.name);
        std::string kind = phase.kind == PhaseKind::Liquid ? "liquid" : "vapor";
        a.text(kind);
        if (kind != "liquid" && kind != "vapor") throw std::invalid_argument("equilibrium callback requires liquid or vapor phases");
        phase.kind = kind == "liquid" ? PhaseKind::Liquid : PhaseKind::Vapor;
        list(a, phase.support, wholes(a)); a.whole(phase.support_declared, 1);
        std::string amount = "free";
        a.text(amount);
        if (amount != "free") throw std::invalid_argument("equilibrium callback requires free phase amounts");
        a.number(phase.amount.fraction); list(a, phase.composition_guess, numbers(a));
        optional(a, phase.packing_guess, numbers(a));
    });
    a.tag("reactions"); list(a, p.reactions, [&a](ReactionSpec& r) { reaction(a, r); });
    a.tag("standard_gibbs_j_per_mol"); list(a, p.standard_gibbs_j_per_mol, numbers(a));
    source_reference(a, p.source_reference);
    if (a.present("solids", !p.solids.empty())) solid_rows(a, p.solids);  // a record without solids keeps its layout
    a.tag("topology_declared"); a.whole(p.topology_declared, 1);
    solver_values(a, p.solver_options);
    active_parameters(a, active);
    thermochemistry(a, t);
    a.tag("continuation_options"); a.number(options.initial_step); a.number(options.minimum_step);
    a.whole(options.max_steps, SIZE_MAX); a.whole(options.max_retries_per_step, SIZE_MAX);  // continue_equilibrium validates both
    a.tag("feed_offset"); list(a, map.offset, numbers(a));
    a.tag("feed_matrix"); list(a, map.matrix, [&a](std::vector<double>& row) { list(a, row, numbers(a)); });
    a.tag("anchor_controls"); list(a, map.anchor, numbers(a));
    a.tag("phase_identity"); a.whole(map.identity_component); a.number(map.minimum_gap);
}

// Property callbacks: the parameter document, the active parameters in model order (their values are the document's),
// a thermochemistry record and the declared solids (each one row of l_s), no state.
template <class A> void property_record(A& a, bool& require_phase_identity, ParameterDocument& d, std::vector<ActiveParameter>& active,
    std::optional<ThermochemistryRecord>& t, policy::SolverOptions& options, std::vector<SolidSpec>& solids,
    std::optional<NeutralReferenceDeclaration>& reference) {
    a.tag("epcsaft.asl.property.v5");
    a.whole(require_phase_identity, 1);
    document(a, d);
    solver_values(a, options);
    active_parameters(a, active);
    thermochemistry(a, t);
    a.tag("solids"); solid_rows(a, solids, &reference);
}

inline std::string asl_record(ParameterDocument d, std::vector<ActiveParameter> active, std::optional<ThermochemistryRecord> t,
    EquilibriumProblem p, ContinuationOptions options, AslFeedMap map, const policy::SolverOptions& solver_options) {
    p.solver_options = p.solver_options.resolved(solver_options);
    AslRecordWriter w;
    equilibrium_record(w, d, active, t, p, options, map);
    return w.finish();
}

inline std::string asl_property_record(ParameterDocument d, std::vector<ActiveParameter> active, std::optional<ThermochemistryRecord> t,
    policy::SolverOptions solver_options, bool require_phase_identity, std::vector<SolidSpec> solids = {},
    std::optional<NeutralReferenceDeclaration> reference = std::nullopt) {
    AslRecordWriter w;
    property_record(w, require_phase_identity, d, active, t, solver_options, solids, reference);
    return w.finish();
}
}  // namespace epcsaft::equilibrium
