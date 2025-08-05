# Gather C sources
solarmain_cs = []
foreach i : run_command('find', '.', '-name', '*.c', check : true).stdout().split('\n')
  if i != ''
    solarmain_cs += files(i)
  endif
endforeach

# Gather C++ sources (.cc)
solarmain_cpps = []
foreach i : run_command('find', '.', '-name', '*.cc', check : true).stdout().split('\n')
  if i != ''
    solarmain_cpps += files(i)
  endif
endforeach

# Gather Fortran sources
solarmain_fs = []
foreach i : run_command('find', '.', '-name', '*.f', check : true).stdout().split('\n')
  if i != ''
    solarmain_fs += files(i)
  endif
endforeach

# Include directories
solarmain_inc = include_directories('../include')
solarmain_eigen_inc = include_directories('../include/Eigen')

# Compiler flags
cflags = [ '-std=c99', '-m64', '-O2', '-DUSE_SAFELIB', '-D_XOPEN_SOURCE=700', '-fexceptions' ]
fflags = [ '-m64', '-O2', '-fno-second-underscore', '-fexceptions', '-fallow-argument-mismatch' ]
cxxflags = [ '-O2', '-std=c++11', '-fopenmp', '-m64', '-DUSE_SAFELIB', '-fexceptions', '-DTR1', '-DLINUX64' ]

# Dependencies — these should be passed from root meson.build, or obtained here via subproject()
# Example (uncomment if needed):
# libtcl_proj = subproject('libtcl')
# libtcl_dep = libtcl_proj.get_variable('libtcl')

# Replace these placeholders with your actual dependency variables from your root meson.build
#libtcl_dep = dependency('libtcl', required: false)  # fallback if not subproject
#libsafe_dep = dependency('libsafe', required: false)
#libricvolume_dep = dependency('libricvolume', required: false)
#libplinkio_dep = dependency('libplinkio', required: false)

executable('solar',
  solarmain_cs,
  solarmain_cpps,
  solarmain_fs,
  include_directories: [solarmain_inc, solarmain_eigen_inc],
  dependencies: [libtcl_dep, libsafe_dep, libricvolume_dep, libplinkio_dep],
  #dependencies: [libricvolume_dep, libplinkio_dep, libsafe_dep],
  c_args: cflags,
  cpp_args: cxxflags,
  fortran_args: fflags,
)
