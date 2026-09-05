#!/bin/bash 

grom=/gromacs/gromacs-2021.2/bin/gmx

index=../../01_min/index.ndx
tpr=../../03_prod_r1/001.tpr
traj=../../04_analysis/concatenated_r1.xtc


file=SELECTED_distances_rep1.dat

## for d=0.0 another frame was used, from pulling (file SELECTED_distances_pull.dat)

i=0

for frame in $(awk '{print $1}' ${file})

do

i=$((i+1))
x=$(printf "%03d" $i)

$grom trjconv -f ${traj} \
      -s ${tpr} \
      -n ${index} \
      -o conf_${x}.gro \
      -dump ${frame} <<EOF
0      
EOF
 
done

