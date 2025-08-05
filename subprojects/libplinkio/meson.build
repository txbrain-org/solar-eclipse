project('libplinkio', ['c'])

c = run_command('find', '.', '-name', '*.c', check: true).stdout().split('\n')

cs = []
foreach i : c # convert to file array
  if i != ''
    cs += files(i)
  endif
endforeach

cpp = run_command('find', '.', '-name', '*.cpp', check: true).stdout().split('\n')
cc = run_command('find', '.', '-name', '*.cc', check: true).stdout().split('\n')

cpps = []
foreach i : cpp
  if i != ''
    cpps += files(i)
  endif
endforeach

ccs = []
foreach i : cc
  if i != ''
    ccs += files(i)
  endif
endforeach

h = run_command('find', '.', '-name', '*.h', check: true).stdout().split('\n')
inc = include_directories('include')
inc_plinkio = include_directories('include/plinkio')

sdict = {
  'srcdir' : meson.current_source_dir(),
  'builddir' : meson.current_build_dir(),
  'c' : c,  # lists of files cannot be passed to summary, so use string
  'h' : h,
  'cpp' : cpp,
  #'include_directory' : inc
}
summary(sdict)

libplinkio = static_library(
  'plinkio',
  cs,
  cpps,
  ccs,
  include_directories: [inc, inc_plinkio],
)

libplinkio_dep = declare_dependency(
  link_with: libplinkio,
  include_directories: [inc, inc_plinkio],
)

libplinkio_dep
