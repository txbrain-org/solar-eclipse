solarmain_c = run_command('find', '.', '-name', '*.c', check: true).stdout().split('\n')
solarmain_cs = []
foreach i : solarmain_c # convert to file array
  if i != ''
    solarmain_cs += files(i)
  endif
endforeach

solarmain_cpp = run_command('find', '.', '-name', '*.cc', check: true).stdout().split('\n')
solarmain_cc = run_command('find', '.', '-name', '*.cc', check: true).stdout().split('\n')
solarmain_cpps = []
foreach i : solarmain_cpp
  if i != ''
    solarmain_cpps += files(i)
  endif
endforeach

solarmain_ccs = []
foreach i : solarmain_cc
  if i != ''
    solarmain_ccs += files(i)
  endif
endforeach

solarmain_f = run_command('find', '.', '-name', '*.f', check: true).stdout().split('\n')
solarmain_fs = []
foreach i : solarmain_f
  if i != ''
    solarmain_fs += files(i)
  endif
endforeach

solarmain_h= run_command('find', '.', '-name', '*.h', check: true).stdout().split('\n')
solarmain_inc= include_directories('../include')
solarmain_inc_str = '../include'

# FFLAGS=-m64  -O2 -fno-second-underscore -fexceptions -fallow-argument-mismatch
# CFLAGS=-std=c99 -m64 -O2 -I/home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/include -DUSE_SAFELIB -fexceptions
# CXXFLAGS=-O2  -std=c++11 -fopenmp  -m64 -I../lib -I/home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/include -DUSE_SAFELIB -fexceptions -DTR1 -DLINUX64 -DUSE_SAFELIB  
# LDFLAGS=-std=c++11 -m64 -O2 -L/home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/lib -pthread -fopenmp -fexceptions

cflags = [ '-std=c99', '-m64', '-O2', '-DUSE_SAFELIB', '-D_XOPEN_SOURCE=700', '-fexceptions' ]
fflags = [ '-m64', '-O2', '-fno-second-underscore', '-fexceptions', '-fallow-argument-mismatch' ]
cxxflags = ['-O2', '-std=c++11', '-fopenmp', '-m64', '-DUSE_SAFELIB', '-fexceptions', '-DTR1', '-DLINUX64', '-DUSE_SAFELIB']

#subprojects_dir = meson.current_build_dir() + '/../../../subprojects'
#libtcl_installdir = subprojects_dir + '/libtcl/install/usr/local'
#libtcl_libdir = libtcl_installdir + '/lib'
#libtcl_includedir = libtcl_installdir + '/include'
#
## /home/wb/Work/nGit/projects/solar_build_test/solar/builddir/src/solarmain/src/../../../subprojects/libtcl/install/usr/local/lib/libtcl8.4.so
#message('subprojects_dir=' + subprojects_dir)
#message('libtcl_installdir=' + libtcl_installdir)
#message(libtcl_libdir + '/libtcl8.4.so')

message(meson.current_source_dir() + '/../installdir/usr/local/include')
installdir_inc = include_directories('../installdir/usr/local/include')
installdir_lib = include_directories('../installdir/usr/local/lib')
installdir_so = meson.current_source_dir() + '/../installdir/usr/local/lib/libtcl8.4.so'
installdir_stub = meson.current_source_dir() + '/../installdir/usr/local/lib/libtclstub8.4.a'

installdir_so_file = files(meson.current_source_dir() + '/../installdir/usr/local/lib/libtcl8.4.so')
ls = run_command('ls', meson.current_source_dir() + '/../installdir' ).stdout().split('\n')
message('solarmain_src: ls=', ls)

fs = import('fs')
if fs.exists(installdir_so)
  message('found installdir_so')
endif
if fs.exists(installdir_stub)
  message('found installdir_stub')
endif

#message('creating solarmain...')
executable('solar',
  solarmain_cs,
  solarmain_cpps,
  solarmain_ccs,
  solarmain_fs,
  include_directories: [ solarmain_inc, installdir_inc, installdir_lib ],
  dependencies: libtcl_dep,
  #dependencies: [ libtcl_dep, libtcl_so_dep, libtcl_stub_dep ],
  link_with: [ libsafe_dep, libricvolume_dep, libplinkio_dep ], # NOTE: maybe this should be named name_static
  c_args: cflags,
  cpp_args: cxxflags,
  fortran_args: fflags,
)

