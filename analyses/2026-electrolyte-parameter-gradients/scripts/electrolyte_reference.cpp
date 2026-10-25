// Independent Born-diameter roots for the F3 electrolyte holds: bisection on the Properties
// value owner, reading the declarations scripts/run.py exports. run.py compiles and runs it.
#include "properties/helmholtz_model.hpp"
#include "surface/parameter_builder.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
using epcsaft::ParameterDocument;

ParameterDocument read_parameter_document(std::istream& in) {
    epcsaft::ParameterDocument doc;
    std::size_t n, sites, entries;
    in >> n >> sites >> entries;
    doc.component_ids.resize(n); doc.site_component.resize(sites); doc.site_multiplicity.resize(sites);
    for (auto& id : doc.component_ids) in >> std::quoted(id);
    for (std::size_t i = 0; i < sites; ++i) in >> doc.site_component[i] >> doc.site_multiplicity[i];
    in >> std::quoted(doc.base_choice) >> std::quoted(doc.association_choice) >> std::quoted(doc.electrolyte_choice)
       >> std::quoted(doc.ionic_dispersion_choice) >> std::quoted(doc.permittivity_choice) >> std::quoted(doc.polar_choice);
    for (std::size_t i = 0; i < entries; ++i) {
        epcsaft::ParameterEntry entry;
        std::size_t width;
        in >> std::quoted(entry.identity) >> std::quoted(entry.field) >> width;
        entry.index.resize(width);
        for (auto& index : entry.index) in >> index;
        in >> entry.magnitude >> std::quoted(entry.unit);
        doc.entries.push_back(entry);
    }
    if (!in) throw std::runtime_error("incomplete parameter declarations");
    return doc;
}

// Independent bisection uses the Properties value owner, with no chart or cached Newton root.
int electrolyte_reference_file(const char* filename) {
    std::ifstream in(filename);
    auto parameters = epcsaft::build_parameters(read_parameter_document(in));
    const std::size_t n = parameters.component_count();
    std::size_t count;
    in >> count;
    std::cout << "ion,stored_diameter,independent_diameter,energy_error,bracket_width\n" << std::setprecision(17);
    for (std::size_t h = 0; h < count; ++h) {
        std::size_t ion, counter, solvent;
        double t, p, target;
        in >> ion >> counter >> solvent >> t >> p >> target;
        if (!in) throw std::runtime_error("incomplete reference inputs");
        const double stored = parameters.born_diameter_angstrom.at(ion);
        std::vector<double> x(n, 0.0); x.at(solvent) = 1.0;
        const auto energy = [&](double d) {
            parameters.born_diameter_angstrom[ion] = d;
            epcsaft::HelmholtzModel model(parameters);
            return epcsaft::GAS_CONSTANT_J_PER_MOL_K * t
                * epcsaft::infinite_dilution_reference(model, t, p, x, {ion, counter}, {1, 1}).log_fugacity[0] - target;
        };
        double lower = 1.0, upper = 8.0, low = energy(lower), mid = 0.0, error = 0.0;
        if (std::signbit(low) == std::signbit(energy(upper))) throw std::runtime_error("independent root outside bracket");
        for (int iteration = 0; iteration < 100; ++iteration) {
            mid = 0.5 * (lower + upper); error = energy(mid);
            if (std::abs(error) <= 1e-7 && upper - lower <= 1e-12) break;
            if (std::signbit(error) == std::signbit(low)) { lower = mid; low = error; } else upper = mid;
        }
        std::cout << ion << ',' << stored << ',' << mid << ',' << error << ',' << upper-lower << '\n';
        parameters.born_diameter_angstrom[ion] = stored;
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    return electrolyte_reference_file(argv[1]);
}
