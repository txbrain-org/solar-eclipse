project('libRicVolume', [ 'c', 'cpp' ])
#subdir('src')

srcdir = meson.current_source_dir()
builddir = meson.current_build_dir()
#message('srcdir=' + srcdir)
#message('builddir=' + builddir)

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

sdict = {
  'srcdir' : meson.current_source_dir(),
  'builddir' : meson.current_build_dir(),
  'c' : c,  # lists of files cannot be passed to summary, so use string
  'h' : h,
  'cpp' : cpp,
}
summary(sdict)

libricvolume = static_library(
  'RicVolume',
  cs,
  cpps,
  ccs,
)

libricvolume_dep = declare_dependency(
  link_with: libricvolume,
)

# meson.override_dependency('libRicVolume', libricvolume_dep) # not needed, see main meson.build for details
libricvolume_dep
