#!/bin/bash

basedir=$(pwd)
for state in 10 20 50 #70 80
	     
do

    cd $state

    cp ../prod_semi_trr.mdp .
    gmx grompp -f prod_semi_trr.mdp -c md_semi.gro -t md_semi.cpt -n index.ndx -p topol.top -o md_semi_cont.tpr
    gmx mdrun -deffnm md_semi_cont -ntmpi 4 -ntomp 6
    
    cd $basedir
    
done

	  
