#
# Makefile
# ========
#
# This Makefile contains the dependencies within the code,
# and the instructions for building the data modules, the main
# model and the various related utility/testing programs.
#
# Machine-specific details are now included via "make.macro", which
# should be linked to a file appropriate to the system being used,
# including details of the compiler and its flags, and the locations
# of the required software libraries: netCDF

include make.macro

# Order of dependencies: 
#    source file (+ any other *.f/*.F files or code segments incorporated via an INCLUDE statement)
#    then this Makefile and the system specific make.macro, parameters.mod
#    then any other modules referenced, 
#    then the netcdf.inc file (if relevant)
#    then source code of any routines called by the (sub)programs in the source file.
# It might be desirable to (eventually) automate this using the "makedepend" utility.

TARGET = CoarseGraining.out

MAKEFILE = Makefile make.macro

MODS =  class_precision.mod nmlsRD_subs.mod ncdfRD_subs.mod ncdfWR_subs.mod NEMO_nmconv.mod task_data.mod task_operations.mod

SRCS =	class_precision.F90 nmlsRD_subs.F90 ncdfRD_subs.F90 ncdfWR_subs.F90 NEMO_nmconv.F90 main.F90 task_data.F90 task_operations.F90

OBJS =  class_precision.o nmlsRD_subs.o ncdfRD_subs.o ncdfWR_subs.o NEMO_nmconv.o main.o task_data.o task_operations.o

main: ${MODS} ${OBJS}
	${FC} ${FFLAGS} -o ../main.out ${OBJS} ${NCLINK} -lm

# Modules:
# -------

# Precision attributes
class_precision.mod: class_precision.F90 ${MAKEFILE}
	${FC} ${FFLAGS} -c class_precision.F90

# Namelist operations
nmlsRD_subs.mod: nmlsRD_subs.F90 ${MAKEFILE} class_precision.mod
	${FC} ${FFLAGS} -c nmlsRD_subs.F90

# netCDF read subroutines
ncdfRD_subs.mod: ncdfRD_subs.F90 ${MAKEFILE} class_precision.mod
	${FC} ${FFLAGS} ${NCINCL} -c ncdfRD_subs.F90

# netCDF write subroutines (not general)
ncdfWR_subs.mod: ncdfWR_subs.F90 ${MAKEFILE} ncdfRD_subs.mod NEMO_nmconv.mod \
        task_data.mod ${NCDINC} 
	${FC} ${FFLAGS} ${NCINCL} -c ncdfWR_subs.F90

# Nemo naming conventions (for easy access to variables in netCDF files)
NEMO_nmconv.mod: NEMO_nmconv.F90 ${MAKEFILE} class_precision.mod
	${FC} ${FFLAGS} -c NEMO_nmconv.F90

# Data structures and initialization for the Task implemented
task_data.mod: task_data.F90 ${MAKEFILE} class_precision.mod
	${FC} ${FFLAGS} -c task_data.F90

# Operations to perform for the task 
task_operations.mod: task_operations.F90 ${MAKEFILE}
	${FC} ${FFLAGS} -c task_operations.F90


# Object files:
# -------------

main.o: main.F90 \
               ${MAKEFILE} ${MODS} ${NCDINC} 
	       ${FC} ${FFLAGS} ${NCINCL} -c main.F90

fullprog.F: ${SRCS}
	cat ${SRCS} > fullprog.F90

# Can also use ${FFLAGS}, which are usually chosen for optimisation
# May sometimes want the ${DEBUGS} options, of course
# Note that when building straight from source files in this way,
# the files need to be ordered consistently with any dependencies

EXECS = main
INTERMEDIATES = *.mod *.o *.l *.lst *__genmod.* *.list *.chek

clean:
	rm -f ${INTERMEDIATES} ${EXECS}

tidy:
	rm -f ${INTERMEDIATES}
