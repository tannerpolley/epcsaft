cppad_configuration = configuration_data()
cppad_configuration.set('compiler_has_conversion_warn', 0)
cppad_configuration.set('cppad_debug_and_release_01', 1)
cppad_configuration.set('use_cplusplus_2017_ok', 1)
cppad_configuration.set('cppad_version', '20250000.3')
cppad_configuration.set('cppad_has_adolc', 0)
cppad_configuration.set('cppad_has_colpack', 0)
cppad_configuration.set('cppad_has_eigen', 0)
cppad_configuration.set('cppad_has_ipopt', 0)
cppad_configuration.set('cppad_deprecated_01', 0)
cppad_configuration.set('cppad_boostvector', 0)
cppad_configuration.set('cppad_cppadvector', 0)
cppad_configuration.set('cppad_stdvector', 1)
cppad_configuration.set('cppad_eigenvector', 0)
cppad_configuration.set('cppad_has_gettimeofday', 0)
cppad_configuration.set('cppad_tape_addr_type', 'unsigned int')
cppad_configuration.set('cppad_is_same_tape_addr_type_size_t', 0)
cppad_configuration.set('cppad_tape_id_type', 'unsigned int')
cppad_configuration.set('cppad_max_num_threads', 48)
cppad_configuration.set(
  'cppad_has_mkstemp',
  cpp.compiles(
    '''
      #include <stdlib.h>
      #include <unistd.h>
      int main() {
        char pattern[] = "/tmp/fileXXXXXX";
        return mkstemp(pattern) < 0;
      }
    ''',
    name: 'CppAD mkstemp',
  ) ? 1 : 0,
)
cppad_configuration.set('cppad_has_tmpnam_s', 0)
cpp_command = cpp.cmd_array()
cppad_configuration.set('cppad_c_compiler_cmd', cpp_command[cpp_command.length() - 1])
cppad_configuration.set(
  'cppad_c_compiler_gnu_flags',
  cpp.get_id() in ['gcc', 'clang'] ? 1 : 0,
)
cppad_configuration.set('cppad_c_compiler_msvc_flags', 0)
cppad_configuration.set('cppad_is_same_unsigned_int_size_t', 0)
cppad_configuration.set(
  'cppad_padding_block_t',
  cpp.sizeof('void *') == 8 ? '' : 'void* note_used_;',
)
configure_file(
  configuration: cppad_configuration,
  input: '../include/cppad/configure.hpp.in',
  output: 'configure.hpp',
)
