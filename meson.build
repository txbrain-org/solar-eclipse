project('subprojects', ['c', 'cpp', 'fortran'], version: '9.0.1')

libtcl_proj = subproject('libtcl', version: '8.4.20')
libtcl_dep = libtcl_proj.get_variable('libtcl')

libricvolume_proj = subproject('libRicVolume')
libricvolume_dep = libricvolume_proj.get_variable('libricvolume_dep')

##libricvolume_dep = dependency('libRicVolume') # this works with using meson.override_dependency

libplinkio_proj = subproject('libplinkio')
libplinkio_dep = libplinkio_proj.get_variable('libplinkio_dep')

libsafe_proj = subproject('libsafe')
libsafe_dep = libsafe_proj.get_variable('libsafe_dep')

subdir('src/solarmain')
