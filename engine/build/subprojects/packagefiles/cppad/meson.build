project('cppad', 'cpp', version: '20250000.3')

cpp = meson.get_compiler('cpp')
subdir('cppad')
cppad_headers_dep = declare_dependency(
  include_directories: include_directories('.', 'include', is_system: true),
)
cppad_temp_file = files('cppad_lib/temp_file.cpp')
