#!/bin/bash

basedir=$(pwd)
for state in 0

do
    cd 10nm_box
    cd $state
    rm *\#
    gmx make_ndx -f md_semi_cont.gro -o lys_groups_cont.ndx<<EOF
r LSN & a CA
r LSN & a NZ
q
EOF
    
    gmx rdf -f md_semi_cont.xtc -n lys_groups_cont.ndx -o rdf_lsn_Nz_cont.xvg -ref 16 -sel 16
    
    cd $basedir
done
