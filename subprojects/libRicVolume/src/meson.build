srcs = run_command('find', '.', '-name', '*.c', check: true).stdout().split('\n')
srcs = ' '.join(srcs) # convert to string
srcs = srcs.replace('./', '')

hdrs = run_command('find', '.', '-name', '*.h', check: true).stdout().split('\n')
hdrs = ' '.join(hdrs) # convert to string
hdrs = hdrs.replace('./', '')

message('hdrs: ' + ' ' .join(hdrs))
message('srcs: ' + ' ' .join(srcs))
message('cwd: ' + meson.current_source_dir())

libricvolume= static_library(
  'RicVolume',
  'RicVolume.cpp',
)
