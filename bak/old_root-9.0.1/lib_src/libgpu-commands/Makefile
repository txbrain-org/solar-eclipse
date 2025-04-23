include sources.mk
.SUFFIXES: .cu .cc .o .h
.cu.o:
	/usr/local/cuda/bin/nvcc -c    -maxrregcount 62  -std=c++11 -g -use_fast_math -Xcompiler -fPIC -Xcompiler -fexceptions -Xcompiler -w -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86 -Xcompiler -I/usr/local/cuda/include -Xcompiler -I/home/brian/solar_versions/solar-eclipse-9.0.1/lib_src/libgpu-commands/include -o $@ $<
.cc.o:
	/usr/bin/g++ -c -g   -I/home/brian/solar_versions/solar-eclipse-9.0.1/lib_src/libgpu-commands/include -pthread -fopenmp -g -std=c++11 -fPIC -w -I/usr/local/cuda/include -o $@ $<
.PHONY: all
all: libgpu-commands.so
libgpu-commands.so :  $(OBJECTS) $(INCLUDES)
	/usr/local/cuda/bin/nvcc  -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_37,code=sm_37 -gencode=arch=compute_50,code=sm_50 -gencode=arch=compute_52,code=sm_52 -gencode=arch=compute_60,code=sm_60 -gencode=arch=compute_61,code=sm_61 -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86 -Xlinker -fexceptions -Xlinker -O3  -shared -o $@ $(OBJECTS)
install:
	mkdir -p /usr/local/lib
	cp libgpu-commands.so /usr/local/lib
clean:
	rm -f $(OBJECTS) libgpu-commands.so
