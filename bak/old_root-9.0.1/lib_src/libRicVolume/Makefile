include sources.mk
.SUFFIXES: .c .cpp .o .h
.cpp.o:
	/usr/bin/g++ -c -O3 -std=c++11 -I$(SOURCE_PATH) -w -o $@ $<
.c.o:
	/usr/bin/gcc -c -O3 -I$(SOURCE_PATH) -w -o $@ $<
.PHONY: all
all: libRicVolume.a
libRicVolume.a:  $(C_OBJECTS) $(CXX_OBJECTS) $(HEADERS)
	ar rcs $@ $(C_OBJECTS) $(CXX_OBJECTS)
install:
	mkdir -p /home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/include
	mkdir -p /home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/lib
	cp $(HEADERS) /home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/include
	cp libRicVolume.a /home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/lib
clean:
	rm -f $(C_OBJECTS) $(CXX_OBJECTS) libRicVolume.a
