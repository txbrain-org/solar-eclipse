srcs = run_command('find', '.', '-name', '*.c', check: true).stdout().split('\n')
srcs = ' '.join(srcs) # convert to string
srcs = srcs.replace('./', '')

hdrs = run_command('find', '.', '-name', '*.h', check: true).stdout().split('\n')
hdrs = ' '.join(hdrs) # convert to string
hdrs = hdrs.replace('./', '')

#message('hdrs: ' + ' ' .join(hdrs))
#message('srcs: ' + ' ' .join(srcs))
message('cwd: ' + meson.current_source_dir())
message('builddir: ' + meson.current_build_dir())

mod = import('unstable-external_project')
p = mod.add_project('configure',
  configure_options: [
    '--enable-threads',
    '--prefix=' + meson.current_build_dir(),
  ],
  verbose: true,
)

tcl_dep = p.dependency('libtcl8.4')
meson.override_dependency('libtcl8.4', tcl_dep)

#libtcl8.4= static_library(
#  'tcl8.4',
#  'safelib.c',
#  include_directories: inc
#)

#$SCRIPT_PATH/unix/configure  --enable-threads  --prefix=$INSTALL_PREFIX
#make -C $SCRIPT_PATH/unix all 
#make -C $SCRIPT_PATH/unix install
