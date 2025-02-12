include sources.mk
.SUFFIXES: .c .cpp .o .h
.cpp.o:
	g++-9 -c -O3 -std=c++11 -I$(SOURCE_PATH) -w -o $@ $<
.c.o:
	gcc-9 -c -O3 -I$(SOURCE_PATH) -w -o $@ $<
.PHONY: all
all: libRicVolume.a
libRicVolume.a:  $(C_OBJECTS) $(CXX_OBJECTS) $(HEADERS)
	ar rcs $@ $(C_OBJECTS) $(CXX_OBJECTS)
install:
	mkdir -p /usr/local/include
	mkdir -p /usr/local/lib
	cp $(HEADERS) /usr/local/include
	cp libRicVolume.a /usr/local/lib
clean:
	rm -f $(C_OBJECTS) $(CXX_OBJECTS) libRicVolume.a
