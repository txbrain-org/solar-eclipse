## https://sourceware.org/autobook/autobook/autobook_25.html#Generated-Output-Files
##
## https://www.gnu.org/prep/standards/standards.html#Makefile-Conventions
##   - 7.2.5 Variables for Installation Directories
##   - 7.2.6 Standard Targets for Users
## 
## https://just.systems/man/en/functions.html

set shell := ["bash", "-c"]

# ANSI color codes
green := "\\033[32m"
yellow := "\\033[33m"
blue := "\\033[34m"
red := "\\033[31m"
bold := "\\033[1m"
reset := "\\033[0m"

prefix := '/usr/local'

builddir := 'builddir'
srcdir := 'src'
target := 'solarmain'

# clean [all]
clean param:
  if [[ {{param}} == 'all' ]]; then \
    meson compile -C {{builddir}} --clean; \
    meson setup --wipe {{builddir}}; \
  else \
    meson compile -C {{builddir}} --clean; \
  fi 

# meson compile -C builddir
build:
  meson compile -C {{builddir}}
  @printf "%s\n"

# execute target
run:
  ./{{builddir}}/{{srcdir}}/{{target}}
  @printf "%s\n"

# meson build && meson run
test: (clean "all") build run
