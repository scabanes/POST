netcdfpath=/planeto/emlmd/netcdf64-4.0.1_ifort
FC=ifort


FFLAGS=-I${netcdfpath}/include
LDFLAGS=-L${netcdfpath}/lib -lnetcdf

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
