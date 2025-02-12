include sources.mk
.SUFFIXES: .f .cc .o .h
FFLAGS=  -fno-second-underscore -fexceptions
CFLAGS=-std=c99   -I/home/brian/solar_versions/solar-eclipse-9.0.1/bin_src/solarmain-9.0.0-test/include -DUSE_SAFELIB -fexceptions
CXXFLAGS=  -std=c++11 -fopenmp  -m64 -I../lib -I/home/brian/solar_versions/solar-eclipse-9.0.1/bin_src/solarmain-9.0.0-test/include -DUSE_SAFELIB -fexceptions -DTR1 -DLINUX64 -DUSE_SAFELIB  
LDFLAGS=-std=c++11   -L/home/brian/solar_versions/solar-eclipse-9.0.1/bin_src/solarmain-9.0.0-test/lib -pthread -fopenmp -fexceptions
.cc.o:
	g++-9 -c $(CXXFLAGS) -o $@ $<
.f.o:
	gfortran-9 -c $(FCFLAGS) -o $@ $<
.c.o:
	gcc-9 -c $(CFLAGS) -o $@ $<
.PHONY: all
all: solarmain
solarmain : $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS)  $(HEADERS)
	gfortran-9 -static -o $@  $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS) $(LDFLAGS) $(LIBRARIES)
install:
	cp solarmain /home/brian/solar_versions/solar-eclipse-9.0.1/bin_src/solarmain-9.0.0-test
clean:
	rm -f $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS)  solarmain
