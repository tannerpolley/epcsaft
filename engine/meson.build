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

eos_inc = include_directories('src/epcsaft/_native')
eos_core = static_library(
  'epcsaft_engine_eos_core',
  cppad_temp_file,
  'src/epcsaft/_native/eos.cpp',
  'src/epcsaft/_native/pure_saturation.cpp',
  cpp_args: native_warnings,
  dependencies: [
    threads_dep,
    ipopt_dep,
    cppad_headers_dep,
  ],
  include_directories: eos_inc,
  install: false,
  pic: true,
)

equilibrium_sources = files(
  'cpp/equilibrium/chemical_equilibrium.cpp',
  'cpp/equilibrium/chemical_equilibrium_compiler.cpp',
  'cpp/equilibrium/phase_direction_compiler.cpp',
  'cpp/equilibrium/chemical_equilibrium_support.cpp',
  'cpp/equilibrium/reactive_finite_phase.cpp',
  'cpp/equilibrium/grepe.cpp',
  'cpp/equilibrium/provider.cpp',
  'cpp/equilibrium/source_reference.cpp',
  'cpp/equilibrium/solved_state_sensitivity.cpp',
  'cpp/equilibrium/phase_equilibrium.cpp',
  'cpp/equilibrium/phase_state.cpp',
  'cpp/equilibrium/coexistence.cpp',
  'cpp/equilibrium/flash.cpp',
  'cpp/equilibrium/phase_stability_state.cpp',
  'cpp/equilibrium/phase_stability_certificates.cpp',
  'cpp/equilibrium/phase_stability.cpp',
  'cpp/equilibrium/phase_stability_stage1.cpp',
  'cpp/equilibrium/phase_stability_stage2.cpp',
  'cpp/equilibrium/phase_stability_stage3.cpp',
  'cpp/equilibrium/phase_stability_stage4.cpp',
  'cpp/equilibrium/phase_stability_stage5.cpp',
  'cpp/equilibrium/phase_stability_stage6.cpp',
  'cpp/equilibrium/phase_stability_stage7.cpp',
)
equilibrium_inc = include_directories('cpp/equilibrium')
equilibrium_core = static_library(
  'epcsaft_engine_equilibrium_core',
  equilibrium_sources,
  highs_conf_file,
  cpp_args: native_warnings,
  dependencies: [
    ipopt_dep,
    highs_dep,
    nlopt_dep,
    eigen_dep,
    boost_headers_dep,
    cppad_headers_dep,
  ],
  include_directories: [equilibrium_inc, eos_inc],
  install: false,
  link_with: eos_core,
  pic: true,
)

interfaces_inc = include_directories('cpp/interfaces')
regression_private_inc = include_directories('cpp/regression')
ceres_kernel_source = files('cpp/regression/ceres_least_squares.cpp')
cpp.compiles(
  ceres_kernel_source[0],
  name: 'owner-neutral Ceres kernel',
  args: native_warnings + [
    '-I' + meson.current_source_dir() / 'cpp/interfaces',
    '-I' + meson.current_source_dir() / 'cpp/regression',
  ],
  dependencies: [ceres_dep, eigen_dep],
  required: true,
)
regression_core = static_library(
  'epcsaft_engine_regression_core',
  'src/epcsaft/regression/native/fit_core.cpp',
  ceres_kernel_source,
  cpp_args: native_warnings,
  dependencies: [
    ceres_dep,
    ipopt_dep,
    eigen_dep,
  ],
  include_directories: [
    interfaces_inc,
    regression_private_inc,
    equilibrium_inc,
    eos_inc,
    include_directories('src/epcsaft/regression/native'),
  ],
  install: false,
  link_with: [eos_core, equilibrium_core],
  pic: true,
)

py.extension_module(
  '_core',
  'src/epcsaft/_native/module.cpp',
  'src/epcsaft/_native/saturation_pybind.cpp',
  'src/epcsaft/regression/native/fit_pybind.cpp',
  'cpp/equilibrium/module.cpp',
  'cpp/equilibrium/chemical_equilibrium_bindings.cpp',
  cpp_args: native_warnings,
  dependencies: [
    pybind11_dep,
    ipopt_dep,
    eigen_dep,
    cppad_headers_dep,
  ],
  include_directories: [
    equilibrium_inc,
    eos_inc,
    include_directories('src/epcsaft/regression/native'),
  ],
  install: true,
  link_args: release_link_args,
  link_with: [eos_core, equilibrium_core, regression_core],
  subdir: 'epcsaft',
)

install_subdir(
  'src/epcsaft',
  exclude_directories: [
    '__pycache__',
    '_native',
    'equilibrium/__pycache__',
    'regression/__pycache__',
    'regression/native',
  ],
  install_dir: py.get_install_dir(),
)
