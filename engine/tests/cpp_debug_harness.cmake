if(NOT DEFINED EPCSAFT_CPP_DEBUG_EXECUTABLE)
    message(FATAL_ERROR "EPCSAFT_CPP_DEBUG_EXECUTABLE is required")
endif()
if(NOT DEFINED EPCSAFT_CPP_DEBUG_MODEL)
    message(FATAL_ERROR "EPCSAFT_CPP_DEBUG_MODEL is required")
endif()

execute_process(
    COMMAND
        "${EPCSAFT_CPP_DEBUG_EXECUTABLE}"
        "--model-config=${EPCSAFT_CPP_DEBUG_MODEL}"
        --temperature-k=250
        --molar-density-mol-per-m3=1000
        --mole-fractions=0.4,0.6
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error_output
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "debug harness failed with ${result}: ${error_output}")
endif()

foreach(expected IN ITEMS
    "pressure_pa="
    "helmholtz_over_rt_reference_amount="
    "parameter_fingerprint=sha256:"
)
    string(FIND "${output}" "${expected}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR "debug harness output lacks '${expected}':\n${output}")
    endif()
endforeach()
