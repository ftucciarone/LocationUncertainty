#!/bin/bash



MONTH=March23
CONFIGURATION=Complete
TARGET_IGRIDA=/srv/tempdd/ftucciar/$MONTH/release-4.0/cfgs/$CONFIGURATION/MY_SRC


echo
echo "            Updating Location Uncertainty routines on configuration: "
echo $TARGET_IGRIDA
echo



# mainfiles
scp *.F90 ftucciar@igrida-frontend:$TARGET_IGRIDA

# Diagnostics
scp DGN/tludgns.F90 ftucciar@igrida-frontend:$TARGET_IGRIDA

# Dynamics
scp DYN/*.F90 ftucciar@igrida-frontend:$TARGET_IGRIDA

# Stochastic Pressure Module
scp STP/*.F90 ftucciar@igrida-frontend:$TARGET_IGRIDA

# Tracer
scp TRA/*.F90 ftucciar@igrida-frontend:$TARGET_IGRIDA

# Noise Generation
scp ../srcNOI/NOI/tlu_*.F90 ftucciar@igrida-frontend:$TARGET_IGRIDA













