#
## Makefile for MDION
#
## Byoungseon Jeon
# Dept. Applied Science, UC.Davis
# # Dec. 15, 2005
# #
# # suffix definition for compiling
# .SUFFIXES: .o .f90
# #
# #FLAGS = ${MPI_COMPILE_FLAGS} -g   -Mbounds -Mchkstk -Mchkptr -Mchkfpstk -Minform=inform -C  #
#   -g  -C
# option for intel fortran compiler
#${MPI_COMPILE_FLAGS} -O3 -parallel -openmp -xW -ipo -no-prec-div -V  
#-warn all,nodec, interfaces -gen_interfaces -traceback -fpe0 -fps
#LIB =  $(MPI_LD_FLAGS) -lmpi # -lmpichfarg  # -lmpi
#LIB = -I${MPI_ROOT}/include -L${MPI_ROOT}/lib -lmpi
.SUFFIXES: .o .f90

#F90 = ~/sw-local/g95-install/bin/powerpc-apple-darwin6.8-g95
F90 = ~/sw-local/ompi_ppc/bin/mpif90
FLAGS = -g
#LIB = -L~/hw/VASP/G5/mpich2-1.0.3_64bit/lib -lmpich
#INCLUDE = -I~/hw/VASP/G5/mpich2-1.0.3_64bit/include
TARGET = octree_run
#
## Object files
OBJTS = data.o tree.o main.o mpienv.o print.o vverlet.o sample.o
#
## generation of executable
${TARGET}:${OBJTS}
	${F90} ${FLAGS} -o ${TARGET} ${OBJTS} ${LIB} ${INCLUDE}
#
## generation of object files
.f90.o:
	${F90} ${FLAGS} -c $< ${LIB} ${INCLUDE}
#
## clean object files and executable
clean:
	rm -rf *.f90~ *.o *.mod core ${TARGET}


