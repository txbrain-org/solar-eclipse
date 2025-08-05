project('libtcl', ['c'], version: '8.4.20')
#subdir('unix')

builddir_root = meson.current_build_dir() + '/../../..'
#message('builddir_root=' + builddir_root)
builddir = meson.current_build_dir()
#message('builddir=' + builddir)
prefix = meson.current_build_dir() + '/installdir/usr/local'
#message('prefix=' + prefix)
srcdir = meson.current_source_dir()
#message('srcdir=' + srcdir)

rc = run_command('bash', '-c',
  'cp -r ' + srcdir + '/*' + ' ' + builddir + '/',
  capture: true, check: true)
message(rc.stdout())

autoconf_tgt = custom_target(
  output: 'configure',
  command: ['bash', '-c', 'cd ' + builddir + '/unix' + ' && autoconf'],
)

configure_opts = [
  '--enable-threads',
  '--prefix=' + prefix,
  '--exec-prefix=' + prefix,
]

configure_tgt = custom_target('Makefile',
  output: 'Makefile',
  command: ['bash', '-c',
    'cd ' + meson.current_build_dir() + '/unix' +
    ' && ' +
    './configure ' + ' '.join(configure_opts)
  ],
  depends: autoconf_tgt,
  build_by_default: true,
)

make_all_tgt= custom_target(
  output: 'make_all', 
  command: ['make', '-C', builddir + '/unix'],
  depends: [configure_tgt],
  build_by_default: true,
)

make_install_output = [
  'tcl.h',
  'tclDecls.h',
  'tclPlatDecls.h',
  'libtcl8.4.so',
  'libtclstub8.4.a'
]

make_install_tgt = custom_target(
  #output: 'installdir',
  output: make_install_output,
  command: ['make', '-C', builddir + '/unix', 'install'],
  depends: make_all_tgt,
  build_by_default: true,
)

cp_installdir_tgt = custom_target(
  output: 'cp_installdir_tgt',
  command: ['cp', '-r', builddir + '/installdir', srcdir + '/../../src/solarmain'],
  depends: make_install_tgt,
)

cp_libtcl_so_tgt = custom_target(
  output: 'cp_libtcl_so_tgt',
  command: ['cp', '-r', builddir + '/installdir/usr/local/lib/libtcl8.4.so', builddir],
  depends: make_install_tgt,
)

cp_libtcl_a_tgt = custom_target(
  output: 'cp_libtcl_a_tgt',
  command: ['cp', '-r', builddir + '/installdir/usr/local/lib/libtclstub8.4.a', builddir],
  depends: make_install_tgt,
)

#mv_installdir_tgt = custom_target(
#  output: 'mv_installdir_tgt',
#  command: [ 'mv', builddir + '/../../src/solarmain/installdir', builddir + '/../../src/solarmain/libtcl' ],
#  depends: cp_installdir_tgt,
#)

# subprojects/libtcl/meson.build:68:0: ERROR: Sandbox violation: Tried to grab file libtcl8.4.so outside current (sub)project.
libtcl = declare_dependency(
  #sources: [ make_install_tgt, cp_installdir_tgt ],
  sources: [ make_install_tgt, cp_libtcl_so_tgt, cp_libtcl_a_tgt ],
)

meson.override_dependency('libtcl', libtcl)

#tgt_make_test = custom_target(
#  output: 'make_test',
#  command: ['make', '-C', builddir + '/unix', 'test'],
#  depends: tgt_make_install,
#  build_by_default: true,
#)
