############################################################

# Makefile for vemPTC
# Author : Sam Yang
# Last updated : 2016/10
# GNU general public license
# Required libraries : refprop, BLAS, LAPACK, fortranposix
#			gnuplotfortran, lnopt, lc++, lCoolProp

############################################################

EXE       	= vemPTC
FC           	= gfortran 
#FCFLAG      	= -O3
FCFLAG		= -O0 -fcheck=all -fbacktrace
LIBS		= -ldl -lc++ -lCoolProp
INCLD		= -I/usr/local/include
MOD		= -I/usr/local/lib

############################################################

OBJECTS =              		\
	etc/cpinterface.o		\
	modules/mod_var.o		\
	modules/mod_funcs.o		\
	mesh/MeshGen.o			\
	mesh/NeighSearch.o		\
	mesh/EleProp.o			\
	mesh/ViewFactor.o		\
	mesh/VisitVTK.o			\
	solver/fcn.o			\
	solver/minpack.o		\
	solver/opkdmain.o		\
	solver/opkda1.o			\
	solver/opkda2.o			\
	etc/printM.o			\
	etc/vtkOut.o			\
	mainvalidation.o		\
	
############################################################

$(EXE) : $(OBJECTS)
	$(FC) $(FCFLAG) -o $@ $^ $(LIBS) $(INCLD) $(MOD)

etc/cpinterface.o : etc/cpinterface.f90
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o etc/cpinterface.f90

modules/mod_var.o : modules/mod_var.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o modules/mod_var.f95

modules/mod_funcs.o : modules/mod_funcs.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o modules/mod_funcs.f95

mesh/MeshGen.o : mesh/MeshGen.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o mesh/MeshGen.f95

mesh/NeighSearch.o : mesh/NeighSearch.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o mesh/NeighSearch.f95

mesh/EleProp.o : mesh/EleProp.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o mesh/EleProp.f95

mesh/ViewFactor.o : mesh/ViewFactor.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o mesh/ViewFactor.f95

mesh/VisitVTK.o : mesh/VisitVTK.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o mesh/VisitVTK.f95

solver/fcn.o : solver/fcn.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o solver/fcn.f95

solver/minpack.o : solver/minpack.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o solver/minpack.f95

solver/opkdmain.o : solver/opkdmain.f
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o solver/opkdmain.f

solver/opkda1.o : solver/opkda1.f
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o solver/opkda1.f

solver/opkda2.o : solver/opkda2.f
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o solver/opkda2.f

etc/printM.o : etc/printM.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o etc/printM.f95

etc/vtkOut.o : etc/vtkOut.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o etc/vtkOut.f95

mainvalidation.o : mainvalidation.f95
	$(FC) $(FCFLAG) $(INCLD) $(MOD) -c -o $*.o mainvalidation.f95


clean :
	rm *.o *.mod modules/*.o solver/*.o solver/newton/*.o \
	 	mesh/*.o etc/*.o ${EXE}

cleanall :
	rm	output/*.vtk

run :
	./vemESRDC