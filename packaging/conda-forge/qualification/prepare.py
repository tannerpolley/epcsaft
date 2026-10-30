from pathlib import Path
import shutil
root=Path.cwd(); engine=root/'engine'; q=root/'packaging/conda-forge/qualification'
shutil.copytree(q,engine/'qualification-939',ignore=shutil.ignore_patterns('candidate.tar.gz'))
p=engine/'meson.build';p.write_text(p.read_text()+'\n'+(q/'windows-qualification.meson').read_text()+'''
executable('qualification-product', 'qualification-939/product.cpp',
  cpp_args: native_warnings + ['-DEPCSAFT_NATIVE_STATIC'],
  include_directories: native_inc, dependencies: [eigen_dep, cppad_headers_dep], install: false, build_by_default: false)
''')

p=engine/'meson.build'
p.write_text(p.read_text()+"""
executable('qualification-tolerance', 'qualification-939/tolerance-caller.cpp',
  cpp_args: native_warnings, include_directories: native_inc,
  dependencies: [ipopt_dep, eigen_dep, cppad_headers_dep], link_with: epcsaft_native, install: false, build_by_default: false)
""")

p=engine/'meson.build'
p.write_text(p.read_text()+"""
executable('qualification-forward-error', 'qualification-939/forward-error.cpp',
 cpp_args: native_warnings + ['-DEPCSAFT_NATIVE_STATIC'],
 include_directories: [native_inc, include_directories('tests/native')],
 dependencies: [threads_dep, ipopt_dep, eigen_dep, ceres_dep, highs_dep, cppad_headers_dep],
 link_with: epcsaft_private, install: false, build_by_default: false)
""")

p=engine/'meson.build'
p.write_text(p.read_text()+"""
executable('qualification-traits', 'qualification-939/traits.cpp', include_directories: native_inc, dependencies: eigen_dep, install: false, build_by_default: false)
""")
p=engine/'meson.build'
p.write_text(p.read_text()+"""
executable('qualification-pure-complete', 'qualification-939/pure-complete.cpp',
 cpp_args: native_warnings + ['-DEPCSAFT_NATIVE_STATIC'],
 include_directories: [native_inc, include_directories('tests/native')],
 dependencies: [threads_dep, ipopt_dep, eigen_dep, ceres_dep, highs_dep, cppad_headers_dep],
 link_with: epcsaft_private, install: false, build_by_default: false)
""")
p=engine/'meson.build'
p.write_text(p.read_text()+"""
executable('qualification-iterate-bounds', 'qualification-939/iterate-bounds.cpp',
 cpp_args: native_warnings + ['-DEPCSAFT_NATIVE_STATIC'],
 include_directories: [native_inc, include_directories('tests/native')],
 dependencies: [threads_dep, ipopt_dep, eigen_dep, ceres_dep, highs_dep, cppad_headers_dep],
 link_with: epcsaft_private, install: false, build_by_default: false)
""")

p=engine/'meson.build'
p.write_text(p.read_text()+"""
executable('qualification-association-reference', 'qualification-939/association-reference.cpp',
 cpp_args: native_warnings + ['-DEPCSAFT_NATIVE_STATIC'],
 include_directories: [native_inc, include_directories('tests/native')],
 dependencies: [threads_dep, ipopt_dep, eigen_dep, ceres_dep, highs_dep, cppad_headers_dep],
 link_with: epcsaft_private, install: false, build_by_default: false)
""")
