project(
  'epcsaft',
  'c',
  'cpp',
  version: '1.1.3',
  meson_version: '>=1.8',
  subproject_dir: 'build/subprojects',
  default_options: [
    'b_ndebug=if-release',
    'buildtype=release',
    'cpp_std=c++17',
    'strip=true',
    'b_staticpic=true',
  ],
)
py = import('python').find_installation('python3', pure: false)
cpp = meson.get_compiler('cpp')
native_warnings = cpp.get_supported_arguments('-Wall', '-Wextra', '-Wpedantic')
# Pinned alignment makes key-solve timing less sensitive to code placement on the measured host, at no measurable cost
# and +1.2% .text (#356; build/issue-356/alignment).
add_project_arguments(cpp.get_supported_arguments('-falign-functions=64'), language: 'cpp')
release_link_args = get_option('buildtype') == 'release' \
  ? cpp.get_supported_link_arguments('-s') \
  : []

threads_dep = dependency('threads')
# Exact: gibbs_minimizer.cpp static_asserts Ipopt 3.11.9 internals.
ipopt_dep = dependency('ipopt', method: 'pkg-config', version: '==3.11.9', include_type: 'system')
eigen_dep = dependency('eigen3', method: 'pkg-config', version: ['>=3.4', '<3.5'], include_type: 'system')
ceres_dep = dependency(
  'Ceres',
  method: 'cmake',
  modules: ['Ceres::ceres'],
  version: ['>=2.2.0', '<2.3.0'],
)

nanobind_paths = run_command(
  py, '-c',
  '''
import nanobind
from pathlib import Path
import sys, tomllib
requirements = tomllib.loads(Path(sys.argv[1]).read_text())["project"]["dependencies"]
requirement = next(d for d in requirements if d.startswith("nanobind=="))
assert nanobind.__version__ == requirement.split("==")[1]
for path in (nanobind.include_dir(), nanobind.source_dir(), Path(nanobind.__file__).parent / "ext/robin_map/include"):
    print(path)
''',
  meson.project_source_root() / '../pyproject.toml',
  check: true,
).stdout().strip().split('\n')
nanobind_dep = declare_dependency(
  include_directories: include_directories(nanobind_paths[0], is_system: true),
)
nanobind_support = static_library(
  'nanobind_support', nanobind_paths[1] / 'nb_combined.cpp',
  dependencies: [nanobind_dep, py.dependency()],
  include_directories: include_directories(nanobind_paths[2], is_system: true),
  cpp_args: ['-fno-strict-aliasing', '-ffunction-sections', '-fdata-sections'],
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
  install: false,
)

cppad_project = subproject('cppad')
cppad_headers_dep = cppad_project.get_variable('cppad_headers_dep')
cppad_temp_file = cppad_project.get_variable('cppad_temp_file')

# tools/build-engine.sh supplies a shared HiGHS build with these options (#319).
highs_dep = dependency('epcsaft-highs', method: 'pkg-config', required: false)
highs_conf_file = []
if not highs_dep.found()
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
endif

# One include root. A layer names its neighbours by path ("properties/density.hpp"),
# which is what makes the include direction checkable.
native_inc = include_directories('native')

# Filled by the build stages; the order is the layer order. Nothing here may
# include nanobind — that is bindings/, and it is the only place.
# contributions/ is header-only: the terms are templated so one source both
# evaluates a double and records a tape, so there is nothing to compile into
# the archive until the properties layer.
superancillary_table = custom_target('pcsaft_superancillary',
  input: 'data/superancillary/pcsaft.json',
  output: 'pcsaft_superancillary_table.hpp',
  command: [py, files('tools/pcsaft_superancillary_table.py'), '@INPUT@', '@OUTPUT@'],
)
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
  'native/equilibrium/one_component_two_density.cpp',
  'native/equilibrium/phase_stability.cpp',
  'native/equilibrium/phase_boundary.cpp',
  'native/equilibrium/implicit_sensitivity.cpp',
)
# Stage 7
regression_sources = files(
  'native/regression/coordinates.cpp',
  'native/regression/observations.cpp',
  'native/regression/fit.cpp',
  'native/regression/minimax.cpp',
)

epcsaft_native = shared_library(
  'epcsaft_native',
  'native/core/trace.cpp',
  properties_sources,
  superancillary_table,
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
    cppad_headers_dep,
  ],
  include_directories: native_inc,
  install: true,
  install_dir: join_paths(py.get_install_dir(), 'epcsaft'),
  build_rpath: '$ORIGIN',
  install_rpath: '$ORIGIN',
  gnu_symbol_visibility: 'default',
)

epcsaft_private = static_library(
  'epcsaft_native',
  objects: epcsaft_native.extract_all_objects(recursive: true),
  install: false,
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
    cppad_headers_dep,
  ],
  include_directories: native_inc,
  link_args: release_link_args + [
    '-Wl,--exclude-libs,ALL',
    '-Wl,--version-script=' + meson.current_source_dir() + '/build/subprojects/packagefiles/asl/asl.exports',
  ],
  link_with: epcsaft_private,
  install: true,
  install_dir: join_paths(py.get_install_dir(), 'epcsaft'),
)

# Stage 3 adds module.cpp and eos_bindings.cpp; stages 4-7 add the rest.
py.extension_module(
  '_core',
  files('native/bindings/module.cpp', 'native/bindings/eos_bindings.cpp',
        'native/bindings/diagnostics_bindings.cpp',
        'native/bindings/equilibrium_bindings.cpp',
        'native/bindings/regression_bindings.cpp'),
  cpp_args: native_warnings,
  # Keep binding -O1 for the measured framework comparison; native calculations stay at -O3.
  override_options: ['optimization=1', 'b_lto=false'],
  dependencies: [
    nanobind_dep,
    ipopt_dep,
    eigen_dep,
    ceres_dep,
    cppad_headers_dep,
  ],
  include_directories: native_inc,
  install: true,
  link_args: release_link_args,
  link_with: [epcsaft_native, nanobind_support],
  build_rpath: '$ORIGIN',
  install_rpath: '$ORIGIN',
  subdir: 'epcsaft',
)

# Python ships alone: the native tree lives outside src/, so there is nothing
# to exclude but caches.
contributions_checks = executable(
  'contributions_checks',
  'tests/native/contributions_checks.cpp',
  cppad_temp_file,
  objects: epcsaft_native.extract_objects('native/core/trace.cpp'),
  cpp_args: native_warnings,
  gnu_symbol_visibility: 'inlineshidden',
  dependencies: [eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  install: false,
  build_by_default: false,
)
test('contributions', contributions_checks)

properties_checks = executable(
  'properties_checks',
  files('tests/native/properties_checks.cpp',
        'tests/native/properties_p1_m3_checks.cpp',
        'tests/native/properties_m4_c2_checks.cpp',
        'tests/native/properties_p2_checks.cpp',
        'tests/native/properties_p3_f9_checks.cpp',
        'tests/native/properties_p4_checks.cpp',
        'tests/native/properties_log_fugacity_checks.cpp',
        'tests/native/properties_temperature_law_checks.cpp',
        'tests/native/properties_cation_mixing_checks.cpp'),
  cpp_args: native_warnings,
  gnu_symbol_visibility: 'inlineshidden',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  build_rpath: '$ORIGIN',
  install: false,
  build_by_default: false,
)
test('properties', properties_checks)

diagnostics_checks = executable(
  'diagnostics_checks',
  'tests/native/diagnostics_checks.cpp',
  cpp_args: native_warnings,
  gnu_symbol_visibility: 'inlineshidden',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  build_rpath: '$ORIGIN',
  install: false,
  build_by_default: false,
)
test('diagnostics', diagnostics_checks)

equilibrium_checks = executable(
  'equilibrium_checks',
  files('tests/native/equilibrium_checks.cpp',
        'tests/native/equilibrium_s1_s2_checks.cpp',
        'tests/native/equilibrium_s6_e2_checks.cpp',
        'tests/native/equilibrium_e5_checks.cpp',
        'tests/native/equilibrium_tpd_checks.cpp',
        'tests/native/equilibrium_caloric_specification_checks.cpp',
        'tests/native/equilibrium_critical_envelope_checks.cpp',
        'tests/native/equilibrium_reactive_caloric_checks.cpp'),
  cpp_args: native_warnings,
  gnu_symbol_visibility: 'inlineshidden',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep, cpp.find_library('quadmath')],
  include_directories: native_inc,
  link_with: epcsaft_native,
  build_rpath: '$ORIGIN',
  install: false,
  build_by_default: false,
)
# One test per case, named equilibrium/<execution_policy>/<case> from the case table,
# so the slow group runs in parallel and check-affected.sh can honor the policy.
equilibrium_cases = import('fs').read('tests/native/equilibrium_checks.cpp').split('{{"native.equilibrium.')
foreach index : range(1, equilibrium_cases.length())
  fields = equilibrium_cases[index].split('"')
  assert(fields[4] in ['Automatic_fast', 'Explicit'], 'bad execution_policy for ' + fields[0])
  test('equilibrium/' + fields[4] + '/' + fields[0], equilibrium_checks,
    args: ['native.equilibrium.' + fields[0]], timeout: 600)
endforeach

# A Python-free caller of the native interface; it includes no test support header.
surface_example = executable(
  'surface_example',
  'tests/native/surface_example.cpp',
  cpp_args: native_warnings,
  gnu_symbol_visibility: 'inlineshidden',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  build_rpath: '$ORIGIN',
  install: false,
  build_by_default: false,
)
test('surface_example', surface_example)

boundary_checks = executable(
  'boundary_checks',
  'tests/native/boundary_checks.cpp',
  cpp_args: native_warnings,
  gnu_symbol_visibility: 'inlineshidden',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  build_rpath: '$ORIGIN',
  install: false,
  build_by_default: false,
)
test('boundary', boundary_checks, timeout: 600)

regression_checks = executable(
  'regression_checks',
  files('tests/native/regression_checks.cpp',
        'tests/native/regression_f1_j1_checks.cpp',
        'tests/native/regression_j2_inventory_checks.cpp',
        'tests/native/regression_reaction_checks.cpp',
        'tests/native/regression_measurement_checks.cpp'),
  cpp_args: native_warnings,
  gnu_symbol_visibility: 'inlineshidden',
  dependencies: [ipopt_dep, eigen_dep, ceres_dep, cppad_headers_dep, cpp.find_library('quadmath')],
  include_directories: native_inc,
  link_with: epcsaft_native,
  build_rpath: '$ORIGIN',
  install: false,
  build_by_default: false,
)
# R1 falsifiers and the #86 density actions in the fast gate; the fitter's cases by explicit selector.
test('regression', regression_checks, args: ['R1'])
test('regression_fit', regression_checks, args: ['R2'], timeout: 600)
test('regression/Explicit/reactive-cp', regression_checks, args: ['reactive-cp'], timeout: 600)

specification_checks = executable(
  'specification_checks',
  'tests/native/specification_checks.cpp',
  cpp_args: native_warnings,
  gnu_symbol_visibility: 'inlineshidden',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  build_rpath: '$ORIGIN',
  install: false,
  build_by_default: false,
)
test('specification', specification_checks)

association_comparison_checks = executable(
  'association_comparison_checks',
  'tests/native/association_comparison_checks.cpp',
  cpp_args: native_warnings,
  gnu_symbol_visibility: 'inlineshidden',
  dependencies: [eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  link_with: epcsaft_native,
  build_rpath: '$ORIGIN',
  install: false,
  build_by_default: false,
)
test('association-comparison', association_comparison_checks, timeout: 600)

install_subdir(
  'src/epcsaft',
  exclude_directories: [
    '__pycache__',
    'equilibrium/__pycache__',
  ],
  install_dir: py.get_install_dir(),
)

executable('reaction_reference_e3', '../analyses/2026-reaction-reference-association/scripts/e3.cpp',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep], include_directories: native_inc,
  link_with: epcsaft_native, build_rpath: '$ORIGIN', install: false, build_by_default: false)

executable('wu_chloride_coverage', '../analyses/2026-mixed-chloride-water-activity/scripts/wu_coverage.cpp',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep], include_directories: native_inc,
  link_with: epcsaft_native, build_rpath: '$ORIGIN', install: false, build_by_default: false)
