if(NOT DEFINED EPCSAFT_CPP_DEBUG_EXECUTABLE)
    message(FATAL_ERROR "EPCSAFT_CPP_DEBUG_EXECUTABLE is required")
endif()

execute_process(
    COMMAND
        "${EPCSAFT_CPP_DEBUG_EXECUTABLE}"
        --temperature-k=250
        --molar-density-mol-per-m3=1000
        --mole-fractions=1
        --segment-count=1
        --sigma-constant-angstrom=3.7039
        --sigma-term-component=
        --sigma-term-amplitude=
        --sigma-term-exponent=
        --dispersion-energy-over-k-k=150.03
        --k-ij=0
        --l-ij=0
        --site-component=
        --site-multiplicity=
        --association-active=
        --association-energy-over-k-k=
        --association-volume=
        --molar-mass-kg-per-mol=0.016043
        --charge-number=0
        --relative-permittivity=1
        --dielectric-ion-suppression-coefficient=1
        --born-diameter-angstrom=1
        --solvation-factor=1
        --ionic-region-relative-permittivity=1
        --electrolyte-model=0
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error_output
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "debug harness failed with ${result}: ${error_output}")
endif()

foreach(expected IN ITEMS
    "pressure_pa="
    "residual_helmholtz="
    "residual_chemical_potential_over_rt=["
    "log_fugacity_coefficient=["
    "density_diagnostics=unavailable"
)
    string(FIND "${output}" "${expected}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR "debug harness output lacks '${expected}':\n${output}")
    endif()
endforeach()
