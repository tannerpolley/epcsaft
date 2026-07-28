cmake_minimum_required(VERSION 3.20)

project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES CXX)

option(
    EPCSAFT_ENABLE_CPP_DEBUG
    "Build the unsupported direct-C++ debug/custom-state harness"
    OFF
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_FINDPYTHON ON)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Threads REQUIRED)

include(FetchContent)
FetchContent_Declare(
    cppad
    GIT_REPOSITORY https://github.com/coin-or/CppAD.git
    GIT_TAG 20250000.3
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR cmake-do-not-configure
)
FetchContent_MakeAvailable(cppad)
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(nlohmann_json)

set(EPCSAFT_CPPAD_INCLUDE_DIR "${cppad_SOURCE_DIR}/include")
set(EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR "")
if(NOT EXISTS "${EPCSAFT_CPPAD_INCLUDE_DIR}/cppad/configure.hpp")
    set(EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR
        "${CMAKE_CURRENT_BINARY_DIR}/generated_cppad_include"
    )
    file(MAKE_DIRECTORY "${EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR}/cppad")
    set(compiler_has_conversion_warn 0)
    set(cppad_debug_and_release_01 1)
    set(use_cplusplus_2017_ok 1)
    set(cppad_version "20250000.3")
    set(cppad_has_adolc 0)
    set(cppad_has_colpack 0)
    set(cppad_has_eigen 0)
    set(cppad_has_ipopt 0)
    set(cppad_deprecated_01 0)
    set(cppad_boostvector 0)
    set(cppad_cppadvector 0)
    set(cppad_stdvector 1)
    set(cppad_eigenvector 0)
    set(cppad_has_gettimeofday 0)
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles([=[
#include <stdlib.h>
#include <unistd.h>
int main() {
    char pattern[] = "/tmp/fileXXXXXX";
    return mkstemp(pattern) < 0;
}
]=] EPCSAFT_CPPAD_HAS_MKSTEMP)
    if(EPCSAFT_CPPAD_HAS_MKSTEMP)
        set(cppad_has_mkstemp 1)
    else()
        set(cppad_has_mkstemp 0)
    endif()
    set(cppad_has_tmpnam_s 0)
    set(cppad_c_compiler_cmd "${CMAKE_CXX_COMPILER}")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(cppad_c_compiler_gnu_flags 1)
    else()
        set(cppad_c_compiler_gnu_flags 0)
    endif()
    set(cppad_c_compiler_msvc_flags 0)
    set(cppad_is_same_unsigned_int_size_t 0)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(cppad_padding_block_t "")
    else()
        set(cppad_padding_block_t "void* note_used_;")
    endif()
    set(cppad_tape_addr_type "unsigned int")
    set(cppad_is_same_tape_addr_type_size_t 0)
    set(cppad_tape_id_type "unsigned int")
    set(cppad_max_num_threads 48)
    configure_file(
        "${EPCSAFT_CPPAD_INCLUDE_DIR}/cppad/configure.hpp.in"
        "${EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR}/cppad/configure.hpp"
        @ONLY
    )
endif()

add_library(
    epcsaft_native_sdk
    STATIC
    "${cppad_SOURCE_DIR}/cppad_lib/temp_file.cpp"
    src/epcsaft/_native/eos.cpp
    src/epcsaft/_native/native_sdk.cpp
    src/epcsaft/_native/native_model.cpp
)
add_library(epcsaft::native_sdk ALIAS epcsaft_native_sdk)
set_target_properties(
    epcsaft_native_sdk
    PROPERTIES EXPORT_NAME native_sdk POSITION_INDEPENDENT_CODE ON
)
target_include_directories(
    epcsaft_native_sdk
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/include>"
        "$<INSTALL_INTERFACE:epcsaft/include>"
    PRIVATE src/epcsaft/_native
)
target_include_directories(
    epcsaft_native_sdk SYSTEM PRIVATE "${EPCSAFT_CPPAD_INCLUDE_DIR}"
)
target_include_directories(
    epcsaft_native_sdk SYSTEM PRIVATE "${nlohmann_json_SOURCE_DIR}/single_include"
)
if(EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR)
    target_include_directories(
        epcsaft_native_sdk SYSTEM PRIVATE "${EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR}"
    )
endif()
target_compile_features(epcsaft_native_sdk PUBLIC cxx_std_17)
target_link_libraries(epcsaft_native_sdk PUBLIC Threads::Threads)

pybind11_add_module(
    _core
    NO_EXTRAS
    src/epcsaft/_native/module.cpp
    src/epcsaft/_native/native_sdk_pybind.cpp
)
target_include_directories(_core PRIVATE src/epcsaft/_native)
target_compile_features(_core PRIVATE cxx_std_17)
target_link_libraries(_core PRIVATE epcsaft_native_sdk)

install(TARGETS _core LIBRARY DESTINATION epcsaft RUNTIME DESTINATION epcsaft)
install(
    TARGETS epcsaft_native_sdk
    EXPORT epcsaftTargets
    ARCHIVE DESTINATION epcsaft/lib
)
install(DIRECTORY src/epcsaft/include/epcsaft DESTINATION epcsaft/include)
install(
    EXPORT epcsaftTargets
    FILE epcsaftTargets.cmake
    NAMESPACE epcsaft::
    DESTINATION epcsaft/cmake
)
include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/epcsaftConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/epcsaftConfig.cmake"
    INSTALL_DESTINATION epcsaft/cmake
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/epcsaftConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/epcsaftConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/epcsaftConfigVersion.cmake"
    DESTINATION epcsaft/cmake
)

if(EPCSAFT_ENABLE_CPP_DEBUG)
    add_executable(
        epcsaft_cpp_debug
        src/epcsaft/_native/cpp_debug_main.cpp
    )
    target_include_directories(epcsaft_cpp_debug PRIVATE src/epcsaft/_native)
    target_compile_features(epcsaft_cpp_debug PRIVATE cxx_std_17)
    target_link_libraries(epcsaft_cpp_debug PRIVATE epcsaft_native_sdk)

    enable_testing()
    add_test(
        NAME epcsaft_cpp_debug_smoke
        COMMAND
            "${CMAKE_COMMAND}"
            "-DEPCSAFT_CPP_DEBUG_EXECUTABLE=$<TARGET_FILE:epcsaft_cpp_debug>"
            "-DEPCSAFT_CPP_DEBUG_MODEL=${CMAKE_CURRENT_SOURCE_DIR}/tests/data/native_models/gross-methane-ethane.json"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_debug_harness.cmake"
    )
endif()
