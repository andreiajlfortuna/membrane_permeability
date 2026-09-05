#!/bin/bash

grom=/gromacs/gromacs-2021.2/bin/gmx

ndx=../../01_min/index.ndx
tpr=../../03_prod_r5/001.tpr


cat conf_*.gro > aux.gro

#################### center and pbc ####################

###### Onetail - last carbon atom of the tail but ONLY one

$grom trjconv -f aux.gro -s $tpr -o aux_2.xtc -n $ndx -center -pbc mol <<EOF
Onetail
System
EOF

################### center and pbc tails ####################

$grom trjconv -f aux_2.xtc -s $tpr -o aux_3.xtc -n $ndx -center -pbc mol <<EOF
C118
System
EOF

###### remove the xy translation of the small molecule ######

$grom trjconv -f aux_3.xtc -s $tpr -o aux_4.xtc -n $ndx -fit transxy -pbc mol <<EOF
MOL
System
EOF


###### Put everything inside the box ######


## create trajectory

$grom trjconv -f aux_4.xtc -s $tpr -o traj.pdb -n $ndx -pbc mol <<EOF
System
EOF

## write gro files

$grom trjconv -f aux_4.xtc -s $tpr -o traj_.gro -sep -n $ndx -pbc mol <<EOF
System
EOF

## rename files from sequencial to distance value for them to be read by slurm_MD.dat file


for i in $(seq 0 1 18)

do 

file=traj_${i}.gro

if [ ${i} -eq 0 ]

then
a=$(echo "0.0")

else

result=$(echo "${a}+0.2" | bc)
a=$(awk -v num=$result 'BEGIN { printf "%.1f", num }')
 

fi

mv traj_${i}.gro traj_${a}.gro

done 
