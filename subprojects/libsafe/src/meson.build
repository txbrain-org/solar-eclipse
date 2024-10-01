srcs = run_command('find', '.', '-name', '*.c', check: true).stdout().split('\n')
srcs = ' '.join(srcs) # convert to string
srcs = srcs.replace('./', '')

inc = include_directories('../include')

#hdrs = run_command('find', '.', '-name', '*.h').stdout().split('\n')
#hdrs = ' '.join(hdrs) # convert to string
#hdrs = hdrs.replace('./', '')

#message('hdrs: ' + ' ' .join(hdrs))
message('srcs: ' + ' ' .join(srcs))
message('cwd: ' + meson.current_source_dir())

libsafe= static_library(
  'safe',
  'safelib.c',
  include_directories: inc
)
