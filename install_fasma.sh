#!/bin/bash

tar -zxvf FASMA/models/apogee_kurucz.tar.gz -C FASMA/models/
tar -zxvf FASMA/models/marcs.tar.gz -C FASMA/models/
echo "Atmosphere models installed in dir: models"
echo -n "Press the type of system you have: 'rh64', 'rh', 'maclap', 'macdesk' "
read answer

# Insert users pwd. Note: Path can not be too long!
sed -i "22s/data/$(pwd | sed 's/\//#/g')\/FASMA\/MOOG\/data/" FASMA/MOOG/Moogsilent.f
sed -i '22s/#/\//g' FASMA/MOOG/Moogsilent.f

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

sed -i "22s/.*/     .  'data\/'/" FASMA/MOOG/Moogsilent.f


if ! command -v FASMA/MOOG/MOOGSILENT &> /dev/null ; then
    echo "MOOGSILENT is not installed properly!"
    exit 1
fi
    echo "MOOGSILENT is installed."

echo "Installing dependencies..."
pip install .
echo "Dependencies installed."
echo "FASMA is successfully installed!"
