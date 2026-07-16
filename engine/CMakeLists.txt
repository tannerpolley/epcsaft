cmake_minimum_required(VERSION 3.20)

project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_FINDPYTHON ON)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

include(FetchContent)
FetchContent_Declare(
    cppad
    GIT_REPOSITORY https://github.com/coin-or/CppAD.git
    GIT_TAG 20250000.3
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR cmake-do-not-configure
)
FetchContent_MakeAvailable(cppad)

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

pybind11_add_module(
    _core
    NO_EXTRAS
    src/epcsaft/_native/eos.cpp
    src/epcsaft/_native/module.cpp
)
target_include_directories(_core PRIVATE src/epcsaft/_native)
target_include_directories(_core SYSTEM PRIVATE "${EPCSAFT_CPPAD_INCLUDE_DIR}")
if(EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR)
    target_include_directories(
        _core SYSTEM PRIVATE "${EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR}"
    )
endif()
target_compile_features(_core PRIVATE cxx_std_17)

install(TARGETS _core LIBRARY DESTINATION epcsaft RUNTIME DESTINATION epcsaft)
