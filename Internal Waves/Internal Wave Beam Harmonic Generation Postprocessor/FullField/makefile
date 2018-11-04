SHELL=/bin/sh

F77_SRC = fftmod.f

F90_SRC =  FieldFilter.f90  Interpolation2D.f90 Pread.f90 filter.f90
	

F77_OBJS = $(patsubst %.f, %.o, $(F77_SRC))
F90_OBJS = $(patsubst %.f90, %.o, $(F90_SRC))
OBJS = $(F90_OBJS) $(F77_OBJS) 

FC90 = f90
FC77 = f90
#LIBS = -lblas -llapack

LIBS =  -L/opt/absoft11.1/lib64 -lreflapack -lrefblas -lU77 -lV77 -lgomp
CMD = Amin_Oglu_Esteban

#
# ABSOFT FORTRAN COMPILER FLAGSET:
FFLAGS77 = -c  -YEXT_NAMES=LCS -YEXT_SFX=_ -N113 -O3 -s
FFLAGS90 = -YEXT_NAMES=LCS -YEXT_SFX=_ -c -N113 -m64 -s -g -f free -O3 

#
# LINKER FLAGS.
LDFLAGS += 


all: $(CMD)

# Link the executable. 
$(CMD):	$(OBJS) 
	$(FC) $(LDFLAGS) -o $(@) $(PROF) $(OBJS) $(LIBS) 

%.o: %.f90 
	$(FC90) -o $@ $(FFLAGS90) $< 

%.o: %.f 
	$(FC77) -o $@ $(FFLAGS77) $< 

clean clena: 
	-rm -f $(OBJS) $(CMD)


