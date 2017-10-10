EXEC1 = oeorbit
##############################################################################
# 
#
#
#MACHINE = LINUX
FC = ifort

FFLAGS =  -g -fPIC -D$(MACHINE) -DLANGUAGE_FORTRAN -O2 -DD #-check
#FFLAGS =  -DLANGUAGE_FORTRAN -O2 -DD

LIBES =  -ljpeg -lz
MATHLIB  = -lm
CRTM_INC = /usr/local/CRTM/crtm_v2.2.1/include
CRTM_LIB = /usr/local/CRTM/crtm_v2.2.1/lib
#CRTM_INC = /usr/local/CRTM/include
#CRTM_LIB = /usr/local/CRTM/lib
INCFLAGS = -I/usr/local/include \
	    -I$(CRTM_INC)
LIBFLAGS = -L$(CRTM_LIB) -lcrtm #-lCRTM
FLAGS = $(INCFLAGS) $(LIBFLAGS) #$(LIBES) $(MATHLIB)

LIBDIR2 = -L/usr/lib64
##############################################################################
SOURCECODE= \
	$(EXEC1).F

OBJECTSC=\
	RSS_RTM_module.o \
	define_oeorbit.o \
	subr_oeorbit.o \
	crtm.o \
	nr_subrs.o \
	oeorbit.o 

############################################################################## 

$(EXEC1):       $(OBJECTSC)
	        $(FC) $(FFLAGS) $(OBJECTSC) $(FLAGS) -o $(EXEC1) #-check all
		rm -f *.o *.mod

tags:   $(SOURCECODE)
	ctags $(SOURCECODE)

clean:
	rm -f *.o *.mod
##############################################################################
#
.suffixes:  .f90 .f .o
.f90.o:
	$(FC) -c  $(FFLAGS) $(INCFLAGS) $<
.f.o:
	$(FC) -c  $(FFLAGS) $(FLAGS) $<
crtm.o: crtm.f90
	$(FC) -c  $(FFLAGS) $(FLAGS) $< 
RSS_RTM_module.o: RSS_RTM_module.f90
	$(FC) -c  $(LIBES) $(LIBDIR2) $<
##############################################################################


