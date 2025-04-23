include sources.mk
.SUFFIXES: .c .o .h
.c.o:
	/usr/bin/gcc -c -O3 -I$(INCLUDE_PATH) -o $@ $<
.PHONY: all
all: libplinkio.a
libplinkio.a :  $(C_OBJECTS) $(HEADERS)
	ar rcs $@ $(C_OBJECTS) 
install:
	mkdir -p /home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/include
	mkdir -p /home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/lib
	cp -R /home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/lib_src/libplinkio/include/* /home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/include
	cp libplinkio.a /home/wburroughs/Work/nGit/projects/solar_build_test/solar/old_root-9.0.1/bin_src/solarmain-9.0.1/lib
clean:
	rm -f $(C_OBJECTS) libplinkio.a
