include sources.mk
.SUFFIXES: .f .cc .o .h
FFLAGS=-m64  -O2 -fno-second-underscore -fexceptions -fallow-argument-mismatch
CFLAGS=-std=c99 -m64 -O2 -I/home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/include -DUSE_SAFELIB -fexceptions
CXXFLAGS=-O2  -std=c++11 -fopenmp  -m64 -I../lib -I/home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/include -DUSE_SAFELIB -fexceptions -DTR1 -DLINUX64 -DUSE_SAFELIB  
LDFLAGS=-std=c++11 -m64 -O2 -L/home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/lib -pthread -fopenmp -fexceptions
.cc.o:
	/usr/bin/g++ -c $(CXXFLAGS) -o $@ $<
.f.o:
	/usr/bin/gfortran -c $(FFLAGS) -o $@ $<
.c.o:
	/usr/bin/gcc -c $(CFLAGS) -o $@ $<
.PHONY: all
all: solarmain
solarmain : $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS)  $(HEADERS)
	/usr/bin/gfortran -o $@  $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS) $(LDFLAGS) $(LIBRARIES)
install:
	cp solarmain release/bin
clean:
	rm -f $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS)  solarmain
