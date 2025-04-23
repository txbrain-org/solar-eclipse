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
cmp_cmds_json := builddir + '/compile_commands.json'


setup:
  @if [[ ! -d {{builddir}} ]]; then \
    meson setup {{builddir}}; \
    printf "%s\n"; \
  fi

# meson compile -C builddir --clean
clean:
  @rm -rf {{builddir}} > /dev/null

# meson compile -C builddir
build: setup
  @meson compile -C {{builddir}} && cp {{cmp_cmds_json}} .
  @printf "%s\n"

# meson compile -v -C builddir
build-verbose: setup
  @meson compile -v -C {{builddir}} && cp {{cmp_cmds_json}} .
  @printf "%s\n"

run: build
  ./{{builddir}}/{{srcdir}}/{{target}}
  @printf "%s\n"

# meson build && meson run
test: clean build run
