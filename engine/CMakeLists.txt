cmake_minimum_required(VERSION 3.24)

project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_SKIP_INSTALL_RPATH TRUE)
set(PYBIND11_FINDPYTHON ON)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM AND NOT CMAKE_CXX_COMPILER_LAUNCHER)
    set(
        CMAKE_CXX_COMPILER_LAUNCHER
        "${CCACHE_PROGRAM}"
        CACHE FILEPATH "Compiler cache used by all Engine native targets"
    )
endif()

option(EPCSAFT_ENABLE_SANITIZERS "Build Engine-owned native code with ASan and UBSan" OFF)

find_package(Python 3.13 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(Threads REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(Ipopt REQUIRED IMPORTED_TARGET ipopt)
find_path(Unwind_INCLUDE_DIR NAMES unwind.h PATH_SUFFIXES libunwind)
find_package(glog CONFIG REQUIRED)
find_package(Ceres 2.2 REQUIRED)

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

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(BUILD_CXX_EXE OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_EXTRAS_LIB OFF CACHE BOOL "" FORCE)
set(ZLIB OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    highs
    URL https://github.com/ERGO-Code/HiGHS/releases/download/v1.15.1/source-archive.tar.gz
    URL_HASH SHA256=30132da2b21899c62f07c4d6566325a3b869a7da4ce61d8b25f833b38e1562b1
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_GetProperties(highs)
if(NOT highs_POPULATED)
    FetchContent_Populate(highs)
    add_subdirectory("${highs_SOURCE_DIR}" "${highs_BINARY_DIR}" EXCLUDE_FROM_ALL)
endif()

FetchContent_Declare(
    boost_headers
    URL https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.bz2
    URL_HASH SHA256=46d9d2c06637b219270877c9e16155cbd015b6dc84349af064c088e9b5b12f7b
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_GetProperties(boost_headers)
if(NOT boost_headers_POPULATED)
    FetchContent_Populate(boost_headers)
endif()

set(NLOPT_CXX ON CACHE BOOL "" FORCE)
set(NLOPT_FORTRAN OFF CACHE BOOL "" FORCE)
set(NLOPT_PYTHON OFF CACHE BOOL "" FORCE)
set(NLOPT_OCTAVE OFF CACHE BOOL "" FORCE)
set(NLOPT_MATLAB OFF CACHE BOOL "" FORCE)
set(NLOPT_GUILE OFF CACHE BOOL "" FORCE)
set(NLOPT_JAVA OFF CACHE BOOL "" FORCE)
set(NLOPT_SWIG OFF CACHE BOOL "" FORCE)
set(NLOPT_LUKSAN OFF CACHE BOOL "" FORCE)
set(NLOPT_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
    nlopt
    URL https://github.com/stevengj/nlopt/archive/refs/tags/v2.11.0.tar.gz
    URL_HASH SHA256=53e552d83e9294d67db37f0f4a23f15933a9ef698485301a18b98b40004cf0de
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_GetProperties(nlopt)
if(NOT nlopt_POPULATED)
    FetchContent_Populate(nlopt)
    add_subdirectory("${nlopt_SOURCE_DIR}" "${nlopt_BINARY_DIR}" EXCLUDE_FROM_ALL)
endif()

set(EPCSAFT_CPPAD_INCLUDE_DIR "${cppad_SOURCE_DIR}/include")
set(EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR "")
if(NOT EXISTS "${EPCSAFT_CPPAD_INCLUDE_DIR}/cppad/configure.hpp")
    set(
        EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR
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
    epcsaft_engine_eos_core
    STATIC
    "${cppad_SOURCE_DIR}/cppad_lib/temp_file.cpp"
    src/epcsaft/_native/eos.cpp
    src/epcsaft/_native/native_transport.cpp
    src/epcsaft/_native/pure_saturation.cpp
)
set_target_properties(epcsaft_engine_eos_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(
    epcsaft_engine_eos_core
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/include"
    PRIVATE src/epcsaft/_native
)
target_include_directories(
    epcsaft_engine_eos_core SYSTEM PRIVATE "${EPCSAFT_CPPAD_INCLUDE_DIR}"
)
target_include_directories(
    epcsaft_engine_eos_core SYSTEM PRIVATE "${nlohmann_json_SOURCE_DIR}/single_include"
)
if(EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR)
    target_include_directories(
        epcsaft_engine_eos_core SYSTEM PRIVATE "${EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR}"
    )
endif()
target_compile_features(epcsaft_engine_eos_core PUBLIC cxx_std_17)
target_link_libraries(
    epcsaft_engine_eos_core
    PUBLIC Threads::Threads
    PRIVATE PkgConfig::Ipopt
)

pybind11_add_module(
    epcsaft_core
    NO_EXTRAS
    src/epcsaft/_native/module.cpp
    src/epcsaft/_native/native_transport_pybind.cpp
    src/epcsaft/_native/saturation_pybind.cpp
    src/epcsaft/regression/native/fit.cpp
)
set_target_properties(epcsaft_core PROPERTIES OUTPUT_NAME _core)
target_include_directories(
    epcsaft_core
    PRIVATE src/epcsaft/_native src/epcsaft/regression/native
)
target_compile_features(epcsaft_core PRIVATE cxx_std_17)
target_link_libraries(epcsaft_core PRIVATE epcsaft_engine_eos_core Ceres::ceres)

add_library(
    epcsaft_engine_equilibrium_core
    STATIC
    cpp/equilibrium/chemical_equilibrium.cpp
    cpp/equilibrium/chemical_equilibrium_compiler.cpp
    cpp/equilibrium/chemical_equilibrium_support.cpp
    cpp/equilibrium/chemical_equilibrium_solver.cpp
    cpp/equilibrium/provider.cpp
    cpp/equilibrium/flash.cpp
    cpp/equilibrium/result_json.cpp
    cpp/equilibrium/held.cpp
    cpp/equilibrium/held_stage_i.cpp
    cpp/equilibrium/held_stage_ii.cpp
    cpp/equilibrium/held_stage_iii.cpp
    cpp/equilibrium/held2.cpp
    cpp/equilibrium/held2_certificates.cpp
    cpp/equilibrium/held2_algorithm.cpp
    cpp/equilibrium/held2_progress.cpp
    cpp/equilibrium/held2_step1.cpp
    cpp/equilibrium/held2_step2.cpp
    cpp/equilibrium/held2_step3.cpp
    cpp/equilibrium/held2_step4.cpp
    cpp/equilibrium/held2_step5.cpp
    cpp/equilibrium/held2_step6.cpp
    cpp/equilibrium/held2_step7.cpp
    cpp/equilibrium/held2_step8.cpp
    cpp/equilibrium/held2_step9.cpp
    cpp/equilibrium/held2_step10.cpp
    cpp/equilibrium/held2_step8_nlp.cpp
)
set_target_properties(epcsaft_engine_equilibrium_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(epcsaft_engine_equilibrium_core PUBLIC cxx_std_17)
target_include_directories(
    epcsaft_engine_equilibrium_core
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/cpp/equilibrium"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/_native"
)
target_include_directories(
    epcsaft_engine_equilibrium_core SYSTEM PRIVATE "${boost_headers_SOURCE_DIR}"
)
target_link_libraries(
    epcsaft_engine_equilibrium_core
    PUBLIC epcsaft_engine_eos_core
    PRIVATE PkgConfig::Ipopt highs nlopt
)

pybind11_add_module(
    epcsaft_engine_equilibrium_native
    MODULE
    cpp/equilibrium/module.cpp
    cpp/equilibrium/chemical_equilibrium_bindings.cpp
)
set_target_properties(epcsaft_engine_equilibrium_native PROPERTIES OUTPUT_NAME _equilibrium)
target_include_directories(
    epcsaft_engine_equilibrium_native
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/_native"
)
target_link_libraries(epcsaft_engine_equilibrium_native PRIVATE epcsaft_engine_equilibrium_core)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(
        EPCSAFT_ENGINE_NATIVE_TARGETS
        epcsaft_engine_eos_core
        epcsaft_core
        epcsaft_engine_equilibrium_core
        epcsaft_engine_equilibrium_native
    )
    foreach(target IN LISTS EPCSAFT_ENGINE_NATIVE_TARGETS)
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endforeach()
    if(EPCSAFT_ENABLE_SANITIZERS)
        foreach(target IN LISTS EPCSAFT_ENGINE_NATIVE_TARGETS)
            target_compile_options(
                ${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer
            )
        endforeach()
        foreach(
            target
            epcsaft_core
            epcsaft_engine_equilibrium_native
        )
            target_link_options(${target} PRIVATE -fsanitize=address,undefined)
        endforeach()
    endif()
endif()

install(
    TARGETS epcsaft_core
    LIBRARY DESTINATION epcsaft
    RUNTIME DESTINATION epcsaft
)
install(
    TARGETS epcsaft_engine_equilibrium_native
    LIBRARY DESTINATION epcsaft/equilibrium
    RUNTIME DESTINATION epcsaft/equilibrium
)
