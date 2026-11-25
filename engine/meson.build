project(
  'epcsaft',
  'c',
  'cpp',
  version: '1.3.3',
  meson_version: '>=1.8',
  subproject_dir: 'build/subprojects',
  default_options: [
    'b_ndebug=if-release',
    'buildtype=release',
    'cpp_std=c++17',
    'strip=true',
    'b_staticpic=true',
    'backend_max_links=1',
  ],
)
py = import('python').find_installation('python3', pure: false)
cpp = meson.get_compiler('cpp')
windows = host_machine.system() == 'windows'
linux = host_machine.system() == 'linux'
native_rpath = {'windows': '', 'darwin': '@loader_path'}.get(host_machine.system(), '$ORIGIN')
# Off Windows, the MSVC linkage model: functions hidden unless EPC_API, types and their typeinfo visible, so an unmarked
# function fails the Linux link as it fails the Windows one (#41).
native_export_args = windows ? ['-DEPCSAFT_NATIVE_BUILD'] : ['-fvisibility-ms-compat']
quadmath_dep = linux ? [cpp.find_library('quadmath')] : []
# One exported symbol hides the statically linked Engine, as the version script does on Linux.
asl_link_args = {'windows': [], 'darwin': ['-Wl,-exported_symbol,_funcadd_ASL']}.get(host_machine.system(), [
  '-Wl,--exclude-libs,ALL',
  '-Wl,--version-script=' + meson.current_source_dir() + '/build/subprojects/packagefiles/asl/asl.exports',
])
native_warnings = windows ? ['/W4', '/EHsc'] : cpp.get_supported_arguments('-Wall', '-Wextra', '-Wpedantic')
# Pinned alignment makes key-solve timing less sensitive to code placement on the measured host, at no measurable cost
# and +1.2% .text (#356; build/issue-356/alignment).
add_project_arguments(windows ? [] : cpp.get_supported_arguments('-falign-functions=64'), language: 'cpp')
release_link_args = linux and get_option('buildtype') == 'release' \
  ? cpp.get_supported_link_arguments('-s') \
  : []

threads_dep = dependency('threads')
# The public amount-step callbacks and equilibrium corpora are qualified on this release.
ipopt_dep = dependency('ipopt', method: 'pkg-config', version: '==3.14.20', include_type: 'system')
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
  cpp_args: windows ? ['/EHsc'] : ['-fno-strict-aliasing', '-ffunction-sections', '-fdata-sections'],
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

# The order is the layer order. Nothing here may include nanobind — that is
# bindings/, and it is the only place. contributions/ is header-only: the terms
# are templated so one source both evaluates a double and records a tape.
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
equilibrium_sources = files(
  'native/equilibrium/stoichiometry.cpp',
  'native/equilibrium/request_compile.cpp',
  'native/equilibrium/row_evaluation.cpp',
  'native/equilibrium/gibbs_minimizer.cpp',
  'native/equilibrium/reaction_activities.cpp',
  'native/equilibrium/one_component_two_density.cpp',
  'native/equilibrium/phase_stability.cpp',
  'native/equilibrium/phase_boundary.cpp',
  'native/equilibrium/implicit_sensitivity.cpp',
  'native/equilibrium/singular_value_decompositions.cpp',
)
regression_sources = files(
  'native/regression/coordinates.cpp',
  'native/regression/observations.cpp',
  'native/regression/fit.cpp',
  'native/regression/minimax.cpp',
)

worker_platform_dep = []
if windows
  worker_platform_dep = [cpp.find_library('advapi32')]
elif linux
  worker_platform_dep = [cpp.find_library('dl')]
endif

native_sources = ['native/core/trace.cpp', properties_sources, superancillary_table, surface_sources,
  equilibrium_sources, regression_sources, cppad_temp_file, highs_conf_file]
native_deps = [threads_dep, ipopt_dep, eigen_dep, ceres_dep, highs_dep, cppad_headers_dep, worker_platform_dep]

epcsaft_native = shared_library(
  'epcsaft_native',
  native_sources,
  cpp_args: native_warnings + native_export_args,
  dependencies: native_deps,
  include_directories: native_inc,
  install: true,
  install_dir: join_paths(py.get_install_dir(), 'epcsaft'),
  build_rpath: native_rpath,
  install_rpath: native_rpath,
)

private_sources = []
if windows
  private_sources = native_sources
  private_options = {
    'cpp_args': native_warnings + ['-DEPCSAFT_NATIVE_STATIC'],
    'dependencies': native_deps,
    'include_directories': native_inc,
  }
else
  private_options = {'objects': epcsaft_native.extract_all_objects(recursive: true)}
endif
epcsaft_private = static_library(
  'epcsaft_native', private_sources,
  kwargs: private_options,
  install: false,
)

asl_export_options = windows ? {'vs_module_defs': files('asl.def')[0]} : {}
asl_library = shared_library(
  'epcsaft_asl',
  'native/bindings/asl_adapter.cpp',
  cpp_args: native_warnings + (windows ? ['-DEPCSAFT_NATIVE_STATIC'] : []),
  dependencies: native_deps,
  include_directories: native_inc,
  link_args: release_link_args + asl_link_args,
  kwargs: asl_export_options,
  link_with: epcsaft_private,
  install: true,
  install_dir: join_paths(py.get_install_dir(), 'epcsaft'),
)

core_extension = py.extension_module(
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
  build_rpath: native_rpath,
  install_rpath: native_rpath,
  subdir: 'epcsaft',
)
if linux
  # An extension module links with unresolved Python symbols, so a binding call to a symbol without EPC_API fails only
  # at import. Relinking its objects with every symbol resolved fails the build the way the Windows link does (#41).
  shared_module('core_export_check', objects: core_extension.extract_all_objects(recursive: false),
    link_with: [epcsaft_native, nanobind_support], dependencies: [py.dependency(embed: true), ipopt_dep, ceres_dep],
    link_args: ['-Wl,-z,defs'], install: false)
endif

contributions_checks = executable(
  'contributions_checks',
  'tests/native/contributions_checks.cpp',
  cppad_temp_file,
  objects: epcsaft_native.extract_objects('native/core/trace.cpp'),
  cpp_args: native_warnings + native_export_args,
  gnu_symbol_visibility: 'inlineshidden',
  dependencies: [eigen_dep, cppad_headers_dep],
  include_directories: native_inc,
  install: false,
  build_by_default: false,
)
test('contributions', contributions_checks)

# Every native test program links the shared library and is built only on request.
check_kwargs = {
  'cpp_args': native_warnings,
  'gnu_symbol_visibility': 'inlineshidden',
  'include_directories': native_inc,
  'link_with': epcsaft_native,
  # CppAD's own temp_file stays hidden in the library; debug builds call it out of line.
  'objects': epcsaft_native.extract_objects(cppad_temp_file),
  'build_rpath': native_rpath,
  'install': false,
  'build_by_default': false,
}

properties_checks = executable(
  'properties_checks',
  files('tests/native/properties_checks.cpp',
        'tests/native/properties_p1_m3_checks.cpp',
        'tests/native/properties_m4_c2_checks.cpp',
        'tests/native/properties_p2_checks.cpp',
        'tests/native/properties_p3_f9_checks.cpp',
        'tests/native/properties_msa_checks.cpp',
        'tests/native/properties_p4_checks.cpp',
        'tests/native/properties_log_fugacity_checks.cpp',
        'tests/native/properties_temperature_law_checks.cpp',
        'tests/native/properties_pathway_weight_checks.cpp',
        'tests/native/properties_cation_mixing_checks.cpp'),
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep] + quadmath_dep,
  kwargs: check_kwargs,
)
# One test per case, named properties/<case> from the case table, as for equilibrium below.
properties_cases = import('fs').read('tests/native/properties_checks.cpp').split('{{"native.properties.')
foreach index : range(1, properties_cases.length())
  name = properties_cases[index].split('"')[0]
  test('properties/' + name, properties_checks, args: ['native.properties.' + name])
endforeach

diagnostics_checks = executable(
  'diagnostics_checks',
  'tests/native/diagnostics_checks.cpp',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  kwargs: check_kwargs,
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
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep] + quadmath_dep,
  kwargs: check_kwargs,
)
# One test per case, named equilibrium/<execution_policy>/<case> from the case table,
# so the slow group runs in parallel and check-affected.sh can honor the policy.
equilibrium_cases = import('fs').read('tests/native/equilibrium_checks.cpp').split('{{"native.equilibrium.')
foreach index : range(1, equilibrium_cases.length())
  fields = equilibrium_cases[index].split('"')
  assert(fields[4] in ['Automatic_fast', 'Explicit'], 'bad execution_policy for ' + fields[0])
  if not linux and (fields[0].startswith('reactive-caloric-E-') or fields[0].startswith('reactive-caloric-I-') or fields[0] == 'reactive-caloric-U' or fields[0].startswith('reactive-caloric-F6-'))
    continue
  endif
  test('equilibrium/' + fields[4] + '/' + fields[0], equilibrium_checks,
    args: ['native.equilibrium.' + fields[0]], timeout: 600)
endforeach

equilibrium_source_iterate_checks = executable(
  'equilibrium_source_iterate_checks',
  'tests/native/equilibrium_source_iterate_checks.cpp',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep] + (linux ? worker_platform_dep : []),
  export_dynamic: linux,
  kwargs: check_kwargs,
)
test('equilibrium/Automatic_fast/mdea-bubble-source-iterate', equilibrium_source_iterate_checks)

# A Python-free caller of the native interface; it includes no test support header.
surface_example = executable(
  'surface_example',
  'tests/native/surface_example.cpp',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  kwargs: check_kwargs,
)
test('surface_example', surface_example)

boundary_checks = executable(
  'boundary_checks',
  'tests/native/boundary_checks.cpp',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  kwargs: check_kwargs,
)
test('boundary', boundary_checks, timeout: 600)

regression_checks = executable(
  'regression_checks',
  files('tests/native/regression_checks.cpp',
        'tests/native/regression_f1_j1_checks.cpp',
        'tests/native/regression_j2_inventory_checks.cpp',
        'tests/native/regression_reaction_checks.cpp',
        'tests/native/regression_measurement_checks.cpp',
        'tests/native/regression_stability_checks.cpp',
        'tests/native/regression_capped_addition_checks.cpp'),
  dependencies: [ipopt_dep, eigen_dep, ceres_dep, cppad_headers_dep] + quadmath_dep,
  kwargs: check_kwargs,
)
# R1 falsifiers and the #86 density actions in the fast gate; the fitter's cases by explicit selector.
test('regression', regression_checks, args: ['R1'])
test('regression/Automatic_fast/stability-requirement', regression_checks, args: ['stability-requirement'])
test('regression/Automatic_fast/nested-support-tangent-plane-action', regression_checks, args: ['nested-support-tangent-plane-action'])
foreach name : ['electrolyte-prediction-origin', 'compiled-result-origin', 'addition-endpoint-admission', 'addition-endpoint-cache', 'supplied-trial-model',
                'equilibrium-refresh-model', 'equilibrium-result-model', 'equilibrium-estimate-admission', 'equilibrium-optimizer-point']
  test('regression/Automatic_fast/' + name, regression_checks, args: [name])
endforeach
test('regression_fit', regression_checks, args: ['R2'], timeout: 600)
test('regression/Explicit/reactive-cp', regression_checks, args: ['reactive-cp'], timeout: 600)

regression_association_site_checks = executable(
  'regression_association_site_checks',
  'tests/native/regression_association_site_checks.cpp',
  dependencies: [ipopt_dep, eigen_dep, ceres_dep, cppad_headers_dep] + quadmath_dep,
  kwargs: check_kwargs,
)
test('regression/Explicit/association-site-fraction', regression_association_site_checks, timeout: 600)

# Tests of the file-local association and stability reuse paths.
foreach native_check : ['association_refinement', 'stability_trial']
  executable_test = executable(
    native_check + '_checks',
    'tests/native/' + native_check + '_checks.cpp',
    cpp_args: native_warnings + (windows ? ['-DEPCSAFT_NATIVE_STATIC'] : []),
    dependencies: native_deps + quadmath_dep,
    include_directories: native_inc,
    link_with: epcsaft_private,
    build_by_default: false,
  )
  if native_check == 'association_refinement'
    test('properties/association-refinement', executable_test)
  else
    foreach name : ['recovery-support', 'relaxed-split-point']
      test('equilibrium/Automatic_fast/' + name, executable_test, args: [name])
    endforeach
  endif
endforeach

asl_pressure_chart_checks = executable(
  'asl_pressure_chart_checks',
  'tests/native/asl_pressure_chart_checks.cpp',
  cpp_args: native_warnings + (windows ? ['-DEPCSAFT_NATIVE_STATIC'] : []),
  dependencies: native_deps + quadmath_dep,
  include_directories: native_inc,
  link_with: epcsaft_private,
  build_by_default: false,
)
test('regression/Explicit/asl-pressure-chart', asl_pressure_chart_checks)

specification_checks = executable(
  'specification_checks',
  'tests/native/specification_checks.cpp',
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep],
  kwargs: check_kwargs,
)
test('specification', specification_checks)

# Python ships alone: the native tree lives outside src/, so there is nothing
# to exclude but caches.
install_subdir(
  'src/epcsaft',
  exclude_directories: [
    '__pycache__',
    'equilibrium/__pycache__',
  ],
  install_dir: py.get_install_dir(),
)


if host_machine.system() == 'windows'
  fit_windows_checks = executable('fit_windows_checks', 'tests/native/fit_windows_checks.cpp',
    cpp_args: native_warnings, include_directories: native_inc, dependencies: cpp.find_library('advapi32'),
    install: false, build_by_default: false)
  test('fit_windows_processes', fit_windows_checks, timeout: 60)
endif
