include sources.mk
.SUFFIXES: .c .o .h
.c.o:
	gcc-9 -c -O3 -w -I$(INCLUDE_PATH) -o $@ $<
.PHONY: all
all: libsafe.a
libsafe.a :  $(C_OBJECTS) $(HEADERS)
	ar rcs $@ $(C_OBJECTS) 
install:
	mkdir -p /usr/local/include
	mkdir -p /usr/local/lib
	cp $(HEADERS) /usr/local/include
	cp libsafe.a /usr/local/lib
clean:
	rm -f $(C_OBJECTS) libsafe.a
