#!/bin/bash 
# Variable to invoke the program GROMACS
grom=/gromacs/gromacs-2021.2/bin/gmx

cp ../01_min/index.ndx .
### create index with P1 and P2

gro=../02_init_r1/init3.gro

np=`grep " P31" ${gro} | wc -l | awk '{print $1/2}'`

echo "[ P1 ]" >> index.ndx
grep " P31" ${gro} | tail -n ${np} | awk '{print substr($0,16,5)}' >> index.ndx

echo "[ P2 ]" >> index.ndx
grep " P31" ${gro} | head -n ${np} | awk '{print substr($0,16,5)}' >> index.ndx


##For each replicate

for rep in $(seq 1 1 5)

do

#gzip -d ../03_prod_r${rep}/001.tpr.gz 



## concatenate trajectory

        $grom trjcat -f ../03_prod_r${rep}/*.xtc \
                    -o concatenated_r${rep}.xtc \
                    -n index.ndx <<EOF
System              
EOF

## get distances of MOL

mkdir rep_${rep}

cd rep_${rep}/

tpr=../../03_prod_r${rep}/001.tpr
index=../index.ndx

${grom} traj -f ../concatenated_r${rep}.xtc \
	     -s ${tpr} \
	     -n ${index} \
	     -ox MOL_${rep}.xvg \
	     -com \
	     -pbc \
	     -nox \
	     -noy \
	     -dt 10 <<EOF
MOL
EOF

sed -i '/@/d' MOL_${rep}.xvg 
sed -i '/#/d' MOL_${rep}.xvg 

## get distances of P1 and P2

${grom} traj -f ../concatenated_r${rep}.xtc \
	     -s ${tpr} \
	     -n ${index} \
	     -ox P1_${rep}.xvg \
	     -com \
	     -pbc \
	     -nox \
	     -noy \
	     -dt 10 <<EOF
P1
EOF

sed -i '/@/d' P1_${rep}.xvg 
sed -i '/#/d' P1_${rep}.xvg 
awk -i inplace '{print $2}' P1_${rep}.xvg 

${grom} traj -f ../concatenated_r${rep}.xtc \
	     -s ${tpr} \
	     -n ${index} \
	     -ox P2_${rep}.xvg \
	     -com \
	     -pbc \
	     -nox \
	     -noy \
	     -dt 10 <<EOF
P2
EOF

sed -i '/@/d' P2_${rep}.xvg 
sed -i '/#/d' P2_${rep}.xvg 
awk -i inplace '{print $2}' P2_${rep}.xvg 

echo "rep ${rep} is done"

dir=/home/afortuna/permeability/01_sims/ligands/literature_compounds/clonidine/EP1/04_analysis

$dir/distance_2_membrane MOL_${rep}.xvg P1_${rep}.xvg P2_${rep}.xvg dist-rep_${rep}



cd ../

done
