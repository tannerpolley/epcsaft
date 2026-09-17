project(
  'epcsaft',
  'c',
  'cpp',
  version: '0.2.0.dev0',
  meson_version: '>=1.8',
  subproject_dir: 'build/subprojects',
  default_options: [
    'b_ndebug=if-release',
    'buildtype=release',
    'cpp_std=c++17',
    'strip=true',
  ],
)
py = import('python').find_installation('python3', pure: false)
cpp = meson.get_compiler('cpp')
native_warnings = cpp.get_supported_arguments('-Wall', '-Wextra', '-Wpedantic')
release_link_args = get_option('buildtype') == 'release' \
  ? cpp.get_supported_link_arguments('-s') \
  : []

threads_dep = dependency('threads')
ipopt_dep = dependency('ipopt', method: 'pkg-config', version: '>=3.11.9')
eigen_dep = dependency('eigen3', method: 'pkg-config', version: '>=3.4')
ceres_dep = dependency(
  'Ceres',
  method: 'cmake',
  modules: ['Ceres::ceres'],
  version: ['>=2.2.0', '<2.3.0'],
)

pybind11_include = run_command(
  py,
  '-c',
  'import pybind11; print(pybind11.get_include())',
  check: true,
).stdout().strip()
pybind11_dep = declare_dependency(
  include_directories: include_directories(pybind11_include, is_system: true),
)

cppad_project = subproject('cppad')
cppad_headers_dep = cppad_project.get_variable('cppad_headers_dep')
cppad_temp_file = cppad_project.get_variable('cppad_temp_file')

boost_project = subproject('boost_headers')
boost_headers_dep = boost_project.get_variable('boost_headers_dep')

highs_project = subproject(
  'highs',
  default_options: [
    'default_library=static',
    'do_install=false',
    'fast_build=true',
    'highsint64=false',
    'use_zlib=disabled',
    'with_pybind11=false',
    'with_tests=false',
  ],
)
highs_dep = highs_project.get_variable('highs_dep')
highs_conf_file = highs_project.get_variable('highs_conf_file')

cmake = import('cmake')
nlopt_options = cmake.subproject_options()
nlopt_options.add_cmake_defines({
  'BUILD_SHARED_LIBS': false,
  'DISABLE_FP_CONTRACT': true,
  'NLOPT_CXX': true,
  'NLOPT_FORTRAN': false,
  'NLOPT_GUILE': false,
  'NLOPT_JAVA': false,
  'NLOPT_LUKSAN': false,
  'NLOPT_MATLAB': false,
  'NLOPT_OCTAVE': false,
  'NLOPT_PYTHON': false,
  'NLOPT_SWIG': false,
  'NLOPT_TESTS': false,
})
nlopt_options.set_install(false)
nlopt_project = cmake.subproject('nlopt', options: nlopt_options)
nlopt_dep = nlopt_project.dependency('nlopt')

# One include root. A layer names its neighbours by path ("properties/density.hpp"),
# which is what makes the include direction checkable.
native_inc = include_directories('native')

# Filled by the build stages; the order is the layer order. Nothing here may
# include pybind11 — that is bindings/, and it is the only place.
# contributions/ is header-only: the terms are templated so one source both
# evaluates a double and records a tape, so there is nothing to compile into
# the archive until the properties layer.
properties_sources = files(
  'native/properties/helmholtz_model.cpp',
  'native/properties/density.cpp',
  'native/properties/infinite_dilution_reference.cpp',
  'native/properties/thermochemistry.cpp',
)
surface_sources = files('native/surface/parameter_builder.cpp')
# Stages 4, 5, 6
equilibrium_sources = files(
  'native/equilibrium/stoichiometry.cpp',
  'native/equilibrium/request_compile.cpp',
  'native/equilibrium/gibbs_minimizer.cpp',
  'native/equilibrium/phase_boundary.cpp',
  'native/equilibrium/implicit_sensitivity.cpp',
)
# Stage 7
regression_sources = files()

epcsaft_native = static_library(
  'epcsaft_native',
  properties_sources,
  surface_sources,
  equilibrium_sources,
  regression_sources,
  cppad_temp_file,
  highs_conf_file,
  cpp_args: native_warnings,
  dependencies: [
    threads_dep,
    ipopt_dep,
    eigen_dep,
    ceres_dep,
    highs_dep,
    nlopt_dep,
    boost_headers_dep,
    cppad_headers_dep,
  ],
  include_directories: native_inc,
  install: false,
  pic: true,
)

asl_library = shared_library(
  'epcsaft_asl',
  'native/bindings/asl_adapter.cpp',
  cpp_args: native_warnings,
  dependencies: [
    threads_dep,
    ipopt_dep,
    eigen_dep,
    ceres_dep,
    highs_dep,
    nlopt_dep,
    boost_headers_dep,
    cppad_headers_dep,
  ],
  include_directories: native_inc,
  link_args: release_link_args + [
    '-Wl,--exclude-libs,ALL',
    '-Wl,--version-script=' + meson.current_source_dir() + '/build/subprojects/packagefiles/asl/asl.exports',
  ],
  link_with: epcsaft_native,
  install: true,
  install_dir: join_paths(py.get_install_dir(), 'epcsaft'),
)

# Stage 3 adds module.cpp and eos_bindings.cpp; stages 4-7 add the rest.
py.extension_module(
  '_core',
  files('native/bindings/module.cpp', 'native/bindings/eos_bindings.cpp',
        'native/bindings/diagnostics_bindings.cpp',
        'native/bindings/equilibrium_bindings.cpp'),
  cpp_args: native_warnings,
  dependencies: [
    pybind11_dep,
    ipopt_dep,
    eigen_dep,
    ceres_dep,
    cppad_headers_dep,
  ],
  include_directories: native_inc,
  install: true,
  link_args: release_link_args,
  link_with: epcsaft_native,
  subdir: 'epcsaft',
)

# Python ships alone: the native tree lives outside src/, so there is nothing
# to exclude but caches.
contributions_checks = executable(
  'contributions_checks',
  'tests/native/contributions_checks.cpp',
  cppad_temp_file,
  cpp_args: native_warnings,
  dependencies: [eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  install: false,
)
test('contributions', contributions_checks)

properties_checks = executable(
  'properties_checks',
  'tests/native/properties_checks.cpp',
  cpp_args: native_warnings,
  dependencies: [eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  install: false,
)
test('properties', properties_checks)

diagnostics_checks = executable(
  'diagnostics_checks',
  'tests/native/diagnostics_checks.cpp',
  cpp_args: native_warnings,
  dependencies: [eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  install: false,
)
test('diagnostics', diagnostics_checks)

properties_floor = executable(
  'properties_floor',
  'tests/native/properties_floor.cpp',
  cpp_args: native_warnings,
  dependencies: [eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  install: false,
)
test(
  'properties_floor', properties_floor,
  suite: 'floor',
  env: {'OMP_NUM_THREADS': '1'},
  timeout: 600,
)

equilibrium_checks = executable(
  'equilibrium_checks',
  'tests/native/equilibrium_checks.cpp',
  cpp_args: native_warnings,
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  install: false,
)
test('equilibrium', equilibrium_checks, timeout: 600)

# A Python-free caller of the native interface; it includes no test support header.
surface_example = executable(
  'surface_example',
  'tests/native/surface_example.cpp',
  cpp_args: native_warnings,
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  install: false,
)
test('surface_example', surface_example)

boundary_checks = executable(
  'boundary_checks',
  'tests/native/boundary_checks.cpp',
  cpp_args: native_warnings,
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  install: false,
)
test('boundary', boundary_checks, timeout: 600)

specification_checks = executable(
  'specification_checks',
  'tests/native/specification_checks.cpp',
  cpp_args: native_warnings,
  dependencies: [eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  install: false,
)
test('specification', specification_checks)

install_subdir(
  'src/epcsaft',
  exclude_directories: [
    '__pycache__',
    'equilibrium/__pycache__',
  ],
  install_dir: py.get_install_dir(),
)
