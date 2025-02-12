include sources.mk
.SUFFIXES: .c .o .h
.c.o:
	gcc-9 -c -O3 -I$(INCLUDE_PATH) -o $@ $<
.PHONY: all
all: libplinkio.a
libplinkio.a :  $(C_OBJECTS) $(HEADERS)
	ar rcs $@ $(C_OBJECTS) 
install:
	mkdir -p /usr/local/include
	mkdir -p /usr/local/lib
	cp -R /home/brian/solar_versions/solar-eclipse-9.0.1/lib_src/libplinkio/include/* /usr/local/include
	cp libplinkio.a /usr/local/lib
clean:
	rm -f $(C_OBJECTS) libplinkio.a
