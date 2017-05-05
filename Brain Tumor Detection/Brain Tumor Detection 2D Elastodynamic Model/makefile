FCOMP = gfortran
FCFLAGS = -g -O2
#LDFLAGS = -L/usr/local/lib-L/usr/local/lib
#LDFLAGS = -L/usr/local/liblas/build -lblas
LDFLAGS = -lblas -llapack 
PROGRAM =  BrainTumor
SRCS = BrainTumor.f90 GLLPoints.f90 NumericalIntegral.f90 Interpolation.f90 MatrixOperations.f90 LinearSolvers.f90 \
   IntegralTransforms.f90 fftmod.f MeshTransform.f90 NumericalDerivative.f90 Elastodynamic.f90 BoundaryConditions.f90 \
   RHS.f90 TimeIntegration.f90


OBJECTS = $(SRCS:.f90=.o)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FCOMP) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM)


