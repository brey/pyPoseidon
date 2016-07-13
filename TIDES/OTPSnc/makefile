ARCH = $(shell uname -s)
ifeq ($(ARCH),SunOS)
 FC = f90
 NCLIB = /home/server/local/lib.solaris2x
 NCINCLUDE = /usr/local/include
endif
ifeq ($(ARCH),Linux)
 FC = gfortran #pgf90
 NCLIB = /usr/local/netcdf4/lib
 NCINCLUDE = /usr/local/include
 NCLIBS= -lnetcdf
endif

predict_tide: predict_tide.f subs.f constit.h
	$(FC) -o predict_tide predict_tide.f subs.f -L$(NCLIB) $(NCLIBS) -I$(NCINCLUDE) 
	#rm *.o
extract_HC:  extract_HC.f subs.f
	$(FC) -o  extract_HC extract_HC.f subs.f -L$(NCLIB) $(NCLIBS) -I$(NCINCLUDE)
	#rm *.o
