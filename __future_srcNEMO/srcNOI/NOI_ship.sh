#!/bin/bash



MONTH=March23
CONFIGURATION=svd_noise
TARGET_IGRIDA=/srv/tempdd/ftucciar/$MONTH/release-4.0/cfgs/$CONFIGURATION/MY_SRC


echo
echo "            Updating Location Uncertainty routines on configuration: "
echo $TARGET_IGRIDA
echo



# mainfiles
scp *.F90 ftucciar@igrida-frontend:$TARGET_IGRIDA

# Noise Generation
scp NOI/*.F90 ftucciar@igrida-frontend:$TARGET_IGRIDA
