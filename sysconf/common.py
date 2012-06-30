import os

bin_dir = 'bin'
cmd_dir = 'cmd'
lib_dir = 'lib'
misc_dir = 'src'
doc_dir = 'doc'
dev_dir = 'dev'
gl = os.path.join ('src', 'opengl', 'gl.h')

cpp_suffix = '.cpp'
h_suffix = '.h'

libname = 'mrtrix'

icon = 'icons/icon.o'
icon_dep = 'icons/icon.rc'

# Qt4 settings:
moc = 'moc'
qt4_modules = [ 'Core', 'Gui', 'OpenGL' ]
qt4_lib_flags = [ '-lGLEW' ]
qt4_cflags = [ '-DQT_SHARED' ] + [ '-DQT_' + entry.upper() + '_LIB' for entry in qt4_modules ]

