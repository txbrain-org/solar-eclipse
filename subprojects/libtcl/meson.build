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

cp_tcl_srcs = run_command('bash', '-c',
  'cp -r ' + srcdir + '/*' + ' ' + builddir + '/',
  capture: true, check: true)
message(cp_tcl_srcs.stdout())

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

libtcl = declare_dependency(
  compile_args: ['-I' + meson.current_build_dir() / 'installdir/usr/local/include'],
  link_args: ['-L' + meson.current_build_dir() / 'installdir/usr/local/lib', '-ltcl8.4']
)

meson.override_dependency('libtcl', libtcl)
