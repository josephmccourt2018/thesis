#!/bin/bash

basedir=$(pwd)
for state in 0 10 20 30 40 50 60 70 80 90 100

do
    cd $state
    gmx trjcat -f md_semi.xtc md_semi_cont.xtc md_semi_cont_2.xtc md_semi_cont_3.xtc md_semi_cont_4.xtc -settime -o combined.xtc<<EOF
0
50000
100000
150000
200000
EOF
    
    cd $basedir
done
