c = run_command('find', '.', '-name', '*.c', check: true).stdout().split('\n')
cs = []
foreach i : c # convert to file array
  if i != ''
    cs += files(i)
  endif
endforeach

cpp = run_command('find', '.', '-name', '*.cpp', check: true).stdout().split('\n')
cc = run_command('find', '.', '-name', '*.cc', check: true).stdout().split('\n')
ccs = []
cpps = []
foreach i : cpp
  if i != ''
    cpps += files(i)
  endif
endforeach
foreach i : cc
  if i != ''
    ccs += files(i)
  endif
endforeach

h = run_command('find', '.', '-name', '*.h', check: true).stdout().split('\n')
hs = []
foreach i : h
  if i != ''
    hs += files(i)
  endif
endforeach

fs = import('fs')
# check if there is an include dir
if fs.is_dir(meson.current_source_dir() + '/../include')
  message('include directory found')
  inc = include_directories('../include')
  inc_str = '../include'
else
  message('include directory not found')
  inc = hs
  inc_str = '../src'
endif

sdict = {
  'srcdir' : meson.current_source_dir(),
  'builddir' : meson.current_build_dir(),
  'includedir' : meson.current_source_dir() + '/../include',
  'c' : c,  # lists of files cannot be passed to summary, so use string
  'h' : h,
  'cpp' : cpp,
  'include_directory' : inc_str
}
summary(sdict)

libsafe = static_library(
  'safe',
  cs,
  cpps,
  ccs,
  include_directories: inc,
)

libsafe_dep = declare_dependency(
  link_with: libsafe,
  include_directories: inc,
)

libsafe_dep
