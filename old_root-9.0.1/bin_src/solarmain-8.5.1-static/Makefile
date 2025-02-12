include sources.mk
.SUFFIXES: .f .cc .o .h
FFLAGS=-m64  -O2 -fno-second-underscore -fexceptions
CFLAGS=-std=c99 -m64 -O2 -I/home/brian/solar_versions/solar-eclipse-9.0.1/bin_src/solarmain-8.5.1-static/include -DUSE_SAFELIB -fexceptions
CXXFLAGS=-O2  -std=c++11 -fopenmp  -m64 -I../lib -I/home/brian/solar_versions/solar-eclipse-9.0.1/bin_src/solarmain-8.5.1-static/include -DUSE_SAFELIB -fexceptions -DTR1 -DLINUX64 -DUSE_SAFELIB  
LDFLAGS=-std=c++11 -m64 -O2 -L/home/brian/solar_versions/solar-eclipse-9.0.1/bin_src/solarmain-8.5.1-static/lib -pthread -fopenmp -fexceptions
.cc.o:
	g++-9 -c $(CXXFLAGS) -o $@ $<
.f.o:
	gfortran-9 -c $(FCFLAGS) -o $@ $<
.c.o:
	gcc-9 -c $(CFLAGS) -o $@ $<
.PHONY: all
all: solarmain
solarmain : $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS)  $(HEADERS)
	gfortran-9 -o $@  $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS) $(LDFLAGS) $(LIBRARIES)
install:
	cp solarmain /home/brian/solar_versions/solar-eclipse-9.0.1/bin_src/solarmain-8.5.1-static
clean:
	rm -f $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS)  solarmain
