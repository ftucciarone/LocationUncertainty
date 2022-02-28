#
# Makefile of coarse-grain
# ========================
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

MAKEFILE = Makefile make.macro

MODS =  class_precision.mod ionml_subs.mod param.mod state.mod ionc_subs.mod coarse_subs.mod

SRCS =	class_precision.F90 ionml_subs.F90 param_data.F90 main.F90 state_data.F90 \
	ionc_subs.F90 coarse_subs.F90

OBJS =  class_precision.o ionml_subs.o param_data.o main.o state_data.o \
	ionc_subs.o coarse_subs.o 

main: ${MODS} ${OBJS}
	${FC} ${FFLAGS} -o main.out ${OBJS} ${NCLINK} -lm

# Modules:
# ========
#


class_precision.mod: class_precision.F90 ${MAKEFILE}
	${FC} ${FFLAGS} -c class_precision.F90

ionml_subs.mod: ionml_subs.F90 ${MAKEFILE}
	${FC} ${FFLAGS} -c ionml_subs.F90


param.mod: param_data.F90 ${MAKEFILE}
	${FC} ${FFLAGS} -c param_data.F90

# Need to permit 80 character lines to allow for explanatory inline comments
# Add this to the FFLAGS for each compiler, as the
# option syntax for this varies between compilers

# Data (storage) modules:
# -----------------------
# Make these first as they have fewest dependencies,
# but many of the program modules depend on them

state.mod: state_data.F90 ${MAKEFILE} param.mod
	${FC} ${FFLAGS} -c state_data.F90

# (Mainly) program modules:
# -------------------------

coarse_subs.mod: coarse_subs.F90 ${MAKEFILE} param.mod
	${FC} ${FFLAGS} -c coarse_subs.F90

ionc_subs.mod: ionc_subs.F90 ${MAKEFILE} param.mod \
        state.mod ${NCDINC} 
	${FC} ${FFLAGS} ${NCINCL} -c ionc_subs.F90

# Object files:
# -------------
# Keep the library routines as objects rather thatn modules for
# ease of swapping between the library and source 

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
