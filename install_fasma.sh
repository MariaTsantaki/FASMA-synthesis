#!/bin/bash

tar -zxvf FASMA/models/apogee_kurucz.tar.gz -C FASMA/models/
tar -zxvf FASMA/models/marcs.tar.gz -C FASMA/models/
echo "Atmosphere models installed in dir: models"
echo -n "Press the type of system you have: 'rh64', 'rh', 'maclap', 'macdesk' "
read answer

MOOGDATAPATH=$PWD/FASMA/MOOG/data/
echo $MOOGDATAPATH

# backup original file
cp -f FASMA/MOOG/Moogsilent.f FASMA/MOOG/_Moogsilent.f

head -n 21 FASMA/MOOG/_Moogsilent.f > FASMA/MOOG/Moogsilent.f

printf $MOOGDATAPATH | awk 'BEGIN{RS="/";c=""} {if($RN>1) {printf ("     . %s\"/%s\"\n",c,$1);c="//"}}' >> FASMA/MOOG/Moogsilent.f

tail -n +23 FASMA/MOOG/_Moogsilent.f>> FASMA/MOOG/Moogsilent.f

#for debiugging
cp -f FASMA/MOOG/Moogsilent.f FASMA/MOOG/__Moogsilent.f

if [ "$answer" == "rh64" ] ;then
    make -C FASMA/MOOG/ -f Makefile.rh64silent clean ; make -C FASMA/MOOG/ -f Makefile.rh64silent
elif [ "$answer" == "rh" ] ;then
    make -C FASMA/MOOG/ -f Makefile.rhsilent clean ; make -C FASMA/MOOG/ -f Makefile.rhsilent
elif [ "$answer" == "maclap" ] ;then
    sed -i '29s/pcl/mac/' FASMA/MOOG/Moogsilent.f
    make -C FASMA/MOOG/ -f Makefile.maclapsilent clean ; make -C FASMA/MOOG/ -f Makefile.maclapsilent
    sed -i '29s/mac/pcl/' FASMA/MOOG/Moogsilent.f
elif [ "$answer" == "macdesk" ] ;then
    sed -i '29s/pcl/mac/' FASMA/MOOG/Moogsilent.f
    make -C FASMA/MOOG/ -f Makefile.macdesksilent clean ; make -C FASMA/MOOG/ -f Makefile.macdesksilent
    sed -i '29s/mac/pcl/' FASMA/MOOG/Moogsilent.f
else
    echo "The input is not valid. Do it again."
fi

#restore original file
cp -f FASMA/MOOG/_Moogsilent.f FASMA/MOOG/Moogsilent.f

if ! command -v FASMA/MOOG/MOOGSILENT &> /dev/null ; then
    echo "MOOGSILENT is not installed properly!"
    exit 1
fi
    echo "MOOGSILENT is installed."

echo "Installing dependencies..."
pip install .
echo "Dependencies installed."
echo "FASMA is successfully installed!"
