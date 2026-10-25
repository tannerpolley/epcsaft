#pragma once

// Shared physical constants and unit bridges. Values are SI unless the name
// says otherwise; each is defined here once and reused by consumers.
// Vacuum permittivity retains the CODATA 2018 central value (https://physics.nist.gov/cuu/pdf/wall_2018.pdf), with relative
// standard uncertainty 1.5e-10; updates require replaying affected evidence.

// Measured or defined physical constants, each defined here and nowhere else (the
// architecture checker rejects a second gas constant). Values are the 2019 SI
// definitions (https://www.bipm.org/en/measurement-units/si-defining-constants), exact where the SI fixes them.

namespace epcsaft {

inline constexpr double ELEMENTARY_CHARGE_C = 1.602176634e-19;
inline constexpr double VACUUM_PERMITTIVITY_F_PER_M = 8.8541878128e-12;
inline constexpr double BOLTZMANN_J_PER_K = 1.380649e-23;
inline constexpr double AVOGADRO_PER_MOL = 6.02214076e23;
inline constexpr double GAS_CONSTANT_J_PER_MOL_K =
    AVOGADRO_PER_MOL * BOLTZMANN_J_PER_K;
inline constexpr double SPEED_OF_LIGHT_M_PER_S = 299792458.0;

inline constexpr double PI = 3.141592653589793238462643383279502884;

// Unit bridges. The equations work in angstrom and Debye because the published
// parameter tables do; these are the only places that conversion appears. Debye in C m: https://cccbdb.nist.gov/dipunitsx.asp
inline constexpr double ANGSTROM_M = 1.0e-10;
inline constexpr double DEBYE_C_M = 1.0e-21 / SPEED_OF_LIGHT_M_PER_S;
// mu^2/(4 pi epsilon_0) in J A^3 for mu in Debye; 1e-19 only under the pre-2019 mu_0.
inline constexpr double DEBYE_SQUARED_J_ANGSTROM3 = DEBYE_C_M * DEBYE_C_M
    / (4.0 * PI * VACUUM_PERMITTIVITY_F_PER_M * ANGSTROM_M * ANGSTROM_M * ANGSTROM_M);

// Segment geometry is written in angstrom, so a molar density in mol/m^3 becomes
// a number density per angstrom^3 with this factor.
inline constexpr double NUMBER_DENSITY_PER_ANGSTROM3 =
    AVOGADRO_PER_MOL * ANGSTROM_M * ANGSTROM_M * ANGSTROM_M;

}  // namespace epcsaft
