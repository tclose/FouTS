from sysconf.common import *

obj_suffix = '.o'
exe_suffix = ''
lib_prefix = 'lib'
lib_suffix = '.dylib'

cpp = [ 'g++', '-c', '$flags$', '$path$', '$src$', '-o', '$obj$' ]
cpp_flags = [ '-Wall', '-DMACOSX', '-fPIC', '-I/data/home/tclose/Code/Tractography/mrtrix_bundle/include', '-mtune=G5' ]
cflags_thread = []

ld = [ 'g++', '$flags$', '$path$', '$obj$', '$mrtrix$', '-o', '$bin$' ]
ld_flags = [ '-lm', '-L/data/home/tclose/Code/Tractography/mrtrix_bundle/lib', '-latlas_mrtrix', '/usr/local/lib/libgslcblas.a', '/usr/local/lib/libgsl.a' ]
ld_flags_lib_prefix = '-l'
libs_thread = [ '-lpthread' ]

ld_lib = [ 'g++', '-dynamiclib', '$flags$', '$obj$', '-o', '$lib$' ]
ld_lib_flags = ['-L/data/home/tclose/Code/Tractography/mrtrix_bundle/lib', '-latlas_mrtrix' ]

cpp_flags_debug = cpp_flags + [ '-O0', '-g', '-D_GLIBCXX_DEBUG' ]
ld_flags_debug = ld_flags + [ '-g' ]
ld_lib_flags_debug = ld_lib_flags + [ '-g' ]

cpp_flags += [ '-DNDEBUG', '-DOPTIMISED', '-DCALCULATE_GRADIENT' ]

cpp_flags_profile = cpp_flags + [ '-g', '-pg' ]
ld_flags_profile = ld_flags + [ '-g', '-pg' ]
ld_lib_flags_profile = ld_lib_flags + [ '-g', '-pg' ]

cpp_flags += [ '-O3' ]

ld_flags_gl = [] 

pkgconfig = [ 'pkg-config' ]
pkgconfig_env = None

# Qt4 settings:
qt4_path = '/Library/Frameworks'
qt4_include_path = [ qt4_path + '/Qt' + entry + '.framework/Headers' for entry in qt4_modules ]
qt4_lib_flags += [ '-L/Library/Frameworks', '-L/System/Library/Frameworks', '-framework', 'Carbon', '-framework', 'AppKit', '-framework', 'ApplicationServices', '-framework', 'OpenGL', '-lz', '-lGLEW' ]
for entry in qt4_modules:
  qt4_lib_flags += [ '-framework', 'Qt'+entry ]
