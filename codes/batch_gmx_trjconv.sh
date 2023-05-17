#!/bin/bash

basedir=$(pwd)

#cd ../3_with_unit_2/bilayer
#gmx trjconv -f md_2.xtc -o md_2_skip10.xtc -skip 10
#cd $basedir

#cd ../4_with_unit_5/bilayer
#gmx trjconv -f md_2.xtc -o md_2_skip10.xtc -skip 10
#cd $basedir

#cd cd 80/4/bilayer/
#gmx trjconv -f md_2.xtc -o md_2_skip10.xtc -skip 10
#cd $basedir

#cd ../5/bilayer
#gmx trjconv -f md_2.xtc -o md_2_skip10.xtc -skip 10
#cd $basedir

for state in 70 80 90 100
	     
do

    cd $state

    cd 5/bilayer/4x_bilayer/
    gmx trjconv -f md_2.xtc -o md_2_skip10.xtc -skip 10
    
    cd $basedir
    
done

