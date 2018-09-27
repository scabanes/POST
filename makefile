netcdfpath=/planeto/emlmd/netcdf64-4.0.1_ifort
spherepackpath=./spherepack3.2
FC=ifort


FFLAGS=-I${netcdfpath}/include -I${spherepackpath}/lib
LDFLAGS=-L${netcdfpath}/lib -lnetcdf -L${spherepackpath}/lib -lspherepack

SRCS= $(wildcard *.f90)
OBJS=$(SRCS:.f90=.o)
EXEC=$(SRCS:.f90=)
TMP=$(SRCS:=~)

all: spectra_analysis test_harmonic

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%: %.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXEC) $(TMP)
