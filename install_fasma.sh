#!/bin/bash

tar -zxvf FASMA/models/apogee_kurucz.tar.gz -C FASMA/models/
tar -zxvf FASMA/models/marcs.tar.gz -C FASMA/models/
echo "Atmosphere models installed in dir: models"
echo "Installing dependencies..."
pip install --user .
echo "Dependencies installed"
echo ""

echo -n "Press the type of system you have: 'rh64', 'rh', 'maclap', 'macdesk' "
read answer
if [ "$answer" == "rh64" ] ;then
    make -C FASMA/MOOG/ -f Makefile.rh64silent clean ; make -C FASMA/MOOG/ -f Makefile.rh64silent
elif [ "$answer" == "rh" ] ;then
    make -C FASMA/MOOG/ -f Makefile.rhsilent clean ; make -C FASMA/MOOG/ -f Makefile.rhsilent
elif [ "$answer" == "maclap" ] ;then
    make -C FASMA/MOOG/ -f Makefile.maclapsilent clean ; make -C FASMA/MOOG/ -f Makefile.maclapsilent
elif [ "$answer" == "macdesk" ] ;then
    make -C FASMA/MOOG/ -f Makefile.macdesksilent clean ; make -C FASMA/MOOG/ -f Makefile.macdesksilent
else
    echo "The input is not valid. Do it again."
fi

if ! command -v FASMA/MOOG/MOOGSILENT &> /dev/null ; then
    echo "MOOGSILENT is not installed properly!"
    exit 1
fi
    echo "MOOGSILENT is installed."
    echo "FASMA is successfully installed!"
