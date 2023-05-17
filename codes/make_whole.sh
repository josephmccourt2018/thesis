#!/bin/bash

basedir=$(pwd)
for state in 10 20 30 40 60 70 80 90

do
    cd $state
    gmx trjconv -f combined.xtc -s npt.tpr -o combined_whole.xtc -pbc whole<<EOF
0
EOF
    
    cd $basedir
done
