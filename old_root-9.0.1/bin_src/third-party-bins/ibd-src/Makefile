EXECS = ibdprep ibdmat ibdmc allfreq relate multipnt mrgibd makeped unknown genfreq #dommsibs mkibdf
all: $(EXECS)

# Mac
# Note: -finit-local-zero not available
FC = gfortran-9
CC = gcc-9
SAFELIB_PATH ?= /home/brian/solar_versions/solar-eclipse-9.0.1/bin_src/solarmain-9.0.0-test
CFLAGS =    -I$(SAFELIB_PATH)/include
FFLAGS = 
LDFLAGS =   
ifneq ($(OSTYPE), linux-gnu)
CFLAGS =   -I$(SAFELIB_PATH)/include
LDFLAGS = 
endif
NOVECTOR = 
SUNMATH =
INSTALLER = cp
INSTALLER2 = $(INSTALL_PATH)

ibdprep: ibdprep.o
	$(CC) $(LDFLAGS) -static   -o ibdprep ibdprep.o -L$(SAFELIB_PATH)/lib -lsafe -lm

ibdmat: ibdmat.o
	$(CC) $(LDFLAGS) -o ibdmat ibdmat.o

ibdmc: ibdmc.o
	$(FC) $(LDFLAGS) -o ibdmc ibdmc.o

allfreq: allfreq.o fallfrq.f
	$(FC) $(LDFLAGS)  $(NOVECTOR) -o allfreq allfreq.o fallfrq.f

genfreq: allfreq.o fgenfrq.f
	$(FC) $(LDFLAGS) $(NOVECTOR) -o genfreq allfreq.o fgenfrq.f

relate: relate.o
	$(CC) $(LDFLAGS) -o relate relate.o

multipnt: dgedifa.o multipnt.o pipeback.o safelib.o
	$(CC) $(LDFLAGS) -o multipnt multipnt.o dgedifa.o pipeback.o safelib.o -L$(SAFELIB_PATH)/lib -lsafe $(SUNMATH) -lm

mrgibd: mrgibd.o pipeback.o safelib.o
	$(CC) $(LDFLAGS) -o mrgibd mrgibd.o pipeback.o safelib.o -L$(SAFELIB_PATH)/lib -lsafe

mkibdf: mkibdf.o pipeback.o safelib.o
	$(CC) $(LDFLAGS) -o mkibdf mkibdf.o pipeback.o safelib.o -L$(SAFELIB_PATH)/lib -lsafe

makeped: makeped.c
	$(CC) -static  -o makeped makeped.c

unknown: unknown.o
	$(CC) $(LDFLAGS) -o unknown unknown.o

install:
	for i in $(EXECS) ; do ( $(INSTALLER) $$i $(INSTALLER2) ) ; done

test:
	for i in $(EXECS) ; do ( /usr/sbin/install -f ../newbin -o $$i ) ; done

clean:
	rm -f $(EXECS) *.o
