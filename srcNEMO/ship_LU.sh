#!/bin/bash

# directory to be shipped
TARGET_DIR=LU_routines
TARGET_PATH=/Users/ftucciar/

# name of the package created
TARGET_PKG=LU_routines

echo
echo "            Updating Location Uncertainty routines on IGRIDA "
echo

if [ -d "$TARGET_PATH/$TARGET_DIR" ]
then

    echo "                                                 Compressing "
    cd $TARGET_PATH
    zip -r  /Users/ftucciar/$TARGET_PKG.zip  $TARGET_DIR -x '*.DS_Store' '*.sh.swp*' 
    echo "                                                    Shipping "
    scp /Users/ftucciar/$TARGET_PKG.zip ftucciar@igrida-frontend:./
    echo "                                                        Done "

else

    echo "Error: Directory $TARGET_PATH$TARGET_DIR does not exist."

fi



