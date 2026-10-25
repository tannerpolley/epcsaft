// Held 2012 Eqs. (11), (22), Table 6: native admitted pair sizes at 298.15 K.
#include "contributions/helmholtz.hpp"
#include "surface/parameter_builder.hpp"
#include <iomanip>
#include <iostream>

int main() {
    std::cout << "series,solvent,ion,native_cross_angstrom,source_cross_angstrom,error_angstrom\n"
              << std::setprecision(17);
    std::string series;
    std::size_t checked=0;
    while (std::cin >> series) {
        epcsaft::ParameterDocument d;
        std::size_t n, count;
        std::cin >> n;
        d.component_ids.resize(n);
        for (auto& id : d.component_ids) std::cin >> id;
        std::cin >> d.base_choice >> d.association_choice >> d.electrolyte_choice
                 >> d.ionic_dispersion_choice >> d.permittivity_choice >> d.polar_choice;
        std::cin >> count;
        d.site_component.resize(count); d.site_multiplicity.resize(count);
        for (std::size_t i=0; i<count; ++i) std::cin >> d.site_component[i] >> d.site_multiplicity[i];
        std::cin >> count;
        d.entries.resize(count);
        for (auto& e : d.entries) {
            std::cin >> std::quoted(e.identity) >> e.field >> count;
            e.index.resize(count);
            for (auto& index : e.index) std::cin >> index;
            std::cin >> e.magnitude >> std::quoted(e.unit);
        }
        const auto p=epcsaft::build_parameters(d);
        const auto sigma=epcsaft::segment_sigma(298.15,p);
        const auto laws=epcsaft::temperature_law_values(298.15,p);
        std::cin >> count;
        for (std::size_t row=0; row<count; ++row) {
            std::size_t s,i; double expected;
            std::cin >> s >> i >> expected;
            if (!std::cin) return 2;
            const double cross=0.5*(sigma[s]+sigma[i])*(1-epcsaft::pair_correction(298.15,p,laws,epcsaft::TemperatureLaw::Lij,s,i));
            std::cout << series << ',' << d.component_ids[s] << ',' << d.component_ids[i] << ',' << cross << ',' << expected << ',' << cross-expected << '\n';
            if (!std::isfinite(cross) || std::abs(cross-expected)>1e-12) return 1;
            ++checked;
        }
    }
    return std::cin.eof() && checked==126 ? 0 : 2;
}
