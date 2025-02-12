include sources.mk
.SUFFIXES: .f .cc .o .h
FCFLAGS=-m64   -fno-second-underscore -fexceptions
CFLAGS=-std=c99  -m64  -I/home/brian/solar_versions/solar-eclipse-9.0.0/bin_src/solarmain-8.5.1/include -DUSE_SAFELIB -fexceptions
CXXFLAGS=   -std=c++11 -fopenmp  -DGPU_TOOLS  -m64 -I../lib -I/home/brian/solar_versions/solar-eclipse-9.0.0/bin_src/solarmain-8.5.1/include -DUSE_SAFELIB -fexceptions -DTR1 -DLINUX64 -DUSE_SAFELIB  
LDFLAGS=-std=c++11 -m64  -static-libgcc -static-libstdc++ -static-libgfortran  -L/home/brian/solar_versions/solar-eclipse-9.0.0/bin_src/solarmain-8.5.1/lib -pthread -fopenmp -fexceptions
.cc.o:
	/usr/bin/g++-7 -c $(CXXFLAGS) -o $@ $<
.f.o:
	/usr/bin/gfortran-7 -c $(FCFLAGS) -o $@ $<
.c.o:
	/usr/bin/gcc-7 -c $(CFLAGS) -o $@ $<
.PHONY: all
all: solarmain
solarmain : $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS)  $(HEADERS)
	/usr/bin/g++ -o $@  $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS) $(LDFLAGS) $(LIBRARIES)
install:
	cp solarmain /home/brian/solar_versions/solar-eclipse-9.0.0/bin_src/solarmain-8.5.1
clean:
	rm -f $(FC_OBJECTS) $(CXX_OBJECTS) $(C_OBJECTS)  solarmain
