#pragma once

#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace epcsaft::detail {

enum class ParameterFamily {
    SegmentCount,
    SigmaConstant,
    DispersionEnergyOverK,
    Kij,
    Lij,
    AssociationEnergyOverK,
    AssociationVolume,
    BornDiameter, DebyeHuckelDiameter, PackingDiameter, IonSuppressionCoefficient,
    BornCDielectric, IonicRegionRelativePermittivity,
    KijReciprocalTemperatureSlope,  // b_ij in k_ij(T) = k_ij + b_ij (1/T - 1/T_ref,ij), K
    DielectricVolume, DielectricTemperature,  // Schreckenberg d_V (m3/mol) and d_T (K) of a neutral solvent
    CationPairMixing, CationTripleMixing,
    TemperatureLawCoefficient,  // a_p or s of a declared temperature law, by coefficient slot
    PathwayWeight,  // b_0, b_alpha, b_T of the declared pathway_weight
};
constexpr std::size_t PARAMETER_FAMILY_COUNT = static_cast<std::size_t>(ParameterFamily::PathwayWeight) + 1;

// The one family -> (pointer, size) storage view, for double and AD alike.
template <class ParameterSet>
auto family_values(ParameterSet& parameters, ParameterFamily family, std::size_t size = 0) {
    const auto view = [size](auto& values) {
        if constexpr (!std::is_const_v<ParameterSet>) if (size) values.resize(size);
        return std::pair{values.data(), values.size()};
    };
    switch (family) {
    case ParameterFamily::BornCDielectric: return std::pair{&parameters.born_c_dielectric, std::size_t{1}};
    case ParameterFamily::IonicRegionRelativePermittivity: return std::pair{&parameters.ionic_region_relative_permittivity, std::size_t{1}};
    case ParameterFamily::PackingDiameter: return view(parameters.packing_diameter_angstrom);
    case ParameterFamily::DebyeHuckelDiameter: return view(parameters.debye_huckel_diameter_angstrom);
    case ParameterFamily::BornDiameter: return view(parameters.born_diameter_angstrom);
    case ParameterFamily::IonSuppressionCoefficient: return view(parameters.ion_suppression_coefficient);
    case ParameterFamily::SegmentCount: return view(parameters.segment_count);
    case ParameterFamily::SigmaConstant: return view(parameters.sigma_constant);
    case ParameterFamily::DispersionEnergyOverK: return view(parameters.dispersion_energy_over_k);
    case ParameterFamily::Kij: return view(parameters.k_ij);
    case ParameterFamily::Lij: return view(parameters.l_ij);
    case ParameterFamily::AssociationEnergyOverK: return view(parameters.association_energy_over_k);
    case ParameterFamily::AssociationVolume: return view(parameters.association_volume);
    case ParameterFamily::KijReciprocalTemperatureSlope: return view(parameters.k_ij_reciprocal_temperature_slope);
    case ParameterFamily::DielectricVolume: return view(parameters.dielectric_volume_m3_per_mol);
    case ParameterFamily::DielectricTemperature: return view(parameters.dielectric_temperature_k);
    case ParameterFamily::TemperatureLawCoefficient: return view(parameters.temperature_law_coefficient);
    case ParameterFamily::CationPairMixing: return view(parameters.cation_pair_mixing);
    case ParameterFamily::CationTripleMixing: return view(parameters.cation_triple_mixing);
    case ParameterFamily::PathwayWeight: return view(parameters.pathway_weight);
    }
    throw std::invalid_argument("invalid parameter family");
}

struct ActiveParameter final {
    ParameterFamily family;
    std::vector<std::size_t> indices;
};

using PredicateSignature = std::vector<unsigned char>;

}  // namespace epcsaft::detail
