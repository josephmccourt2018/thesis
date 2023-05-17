#!/bin/bash

#basedir=$(pwd)

#SETTING UP BILAYER SYSTEM FOR THE z-NORMAL MEMBRANE

basedir=$(pwd)

molecules=20 #number of molecules in 1 bilayer

#states of charge percent
for state in 10 20 30 40 50 60 70 80 90

do
    mkdir $state
    cd $state
    ProteinL=$(( state*molecules/100 ))
    ProteinLN=$(( molecules - $ProteinL ))
    #copy all necessary files for the proceeding commands
    #itp
    cp ../../$state/*.itp .
    #mdp
    cp ../../$state/*.mdp .
    cp ../../prod_semi_trr_2.mdp .
    #pdb
    cp ../../$state/*.pdb .
    #charmm36 ff
    cp -r ../../$state/*.ff .
    #top
    cp ../../$state/*.top .
    #pl
    cp ../../$state/*.pl .
    #bilayer stucture files
    cp ../../$state/bilayer.gro .
  

    #make a new box that for the system that will be solvated
    echo 0 0 | gmx trjconv -f bilayer.gro -s bilayer.gro -o bilayer_newbox.gro -box -1 -1 10
    
    #solvate, delete waters in tail region, fix topology minimize, generate index file
    gmx editconf -f bilayer_newbox.gro -o bilayer_centered.gro -c

    sed -i '$d' topol.top
    sed -i '$d' topol.top
    
    gmx solvate -cp bilayer_centered.gro -p topol.top -o solvated.gro
    perl water_deletor.pl -in solvated.gro -out solvated_fixed.gro -ref N -middle CA -nwater 3
    
    #######################################################################
    #need to enter in deleted waters manually (fix topology file)
    #######################################################################
    echo Did you fix topology file? Enter y when ready
    read ready
    echo $ready

    
    gmx grompp -f ions.mdp -c solvated_fixed.gro -p topol.top -o ions.tpr
    gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
    gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
    gmx mdrun -deffnm em
 
    gmx make_ndx -f em.gro -o index.ndx
    
    cd $basedir
done
