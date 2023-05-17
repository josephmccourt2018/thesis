#!/bin/bash

basedir=$(pwd)


for tails in 3 4 5
	     
do
    cd 10
    cd $tails

    cd bilayer/4x_bilayer/
    gmx make_ndx -f md_2.gro -o dens_groups.ndx<<EOF
a CA | a A1 | a A2
name 6 HEADS
a C1 | a C2 | a C3 | a C4 | a C5
name 7 TAILS
q
EOF

    gmx density -f md_2_skip10.xtc -s md_2.tpr -n dens_groups.ndx -sl 100 -o pw.xvg -symm yes -center yes<<EOF
0
4
EOF
    
    gmx density -f md_2_skip10.xtc -s md_2.tpr -n dens_groups.ndx -sl 100 -o ions.xvg -symm yes -center yes<<EOF
0
5
EOF

    gmx density -f md_2_skip10.xtc -s md_2.tpr -n dens_groups.ndx -sl 100 -o heads.xvg -symm yes -center yes<<EOF
0
6
EOF

    gmx density -f md_2_skip10.xtc -s md_2.tpr -n dens_groups.ndx -sl 100 -o tails.xvg -symm yes -center yes<<EOF
0
7
EOF
    
    cd $basedir
    
done
