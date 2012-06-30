from sysconf.common import *

obj_suffix = '.o'
exe_suffix = ''
lib_prefix = 'lib'
lib_suffix = '.so'


cpp = [ 'g++', '-c', '$flags$', '$path$', '$src$', '-o', '$obj$' ]
cpp_flags = [ '-Wall', '-march=x86-64', '-fPIC', '-Wno-deprecated']
cflags_thread = []

ld = [ 'g++', '$flags$', '$path$', '$obj$', '$mrtrix$', '-o', '$bin$' ]
ld_flags = [ '-lm', '-L/usr/lib/sse2', '-lcblas', '-latlas', '-llapack_atlas', '-lgfortran', '-L/usr/local/lib/' '-lgslcblas', '-lgsl' ]
ld_flags_lib_prefix = '-l'
libs_thread = [ '-lpthread' ]

ld_lib = [ 'g++', '-shared', '$flags$', '$obj$', '-o', '$lib$' ]
ld_lib_flags = []

cpp_flags_debug = cpp_flags + [ '-O0', '-g3', '-D_GLIBCXX_DEBUG' ]
ld_flags_debug = ld_flags + [ '-g3' ]
ld_lib_flags_debug = ld_lib_flags + [ '-g3' ]

cpp_flags += [ '-DNDEBUG', '-DOPTIMISED']

cpp_flags_profile = cpp_flags + [ '-g', '-pg' ]
ld_flags_profile = ld_flags + [ '-g', '-pg' ]
ld_lib_flags_profile = ld_lib_flags + [ '-g', '-pg' ]

cpp_flags += [ '-O3' ]

ld_flags_gl = [] 

pkgconfig = [ 'pkg-config' ]
pkgconfig_env = None

# Qt4 settings:
qt4_path = '/usr/include/qt4'
qt4_include_path = [ qt4_path ] + [ qt4_path + '/Qt' + entry for entry in qt4_modules ]
qt4_lib_flags += [ '-lQt' + entry for entry in qt4_modules ]

