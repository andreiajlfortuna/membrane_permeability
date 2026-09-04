
grom=/gromacs/gromacs-2021.2/bin/gmx


for rep in $(seq 1 1 5)

do

cd rep_${rep}/

tpr=../../03_prod_r${rep}/001.tpr

index=../index.ndx

$grom distance -f ../concatenated_r${rep}.xtc \
               -s ${tpr} \
               -n ${index} \
               -oxyz distances_C118_MOL_rep_${rep}.xvg \
               -select 'com of group "C118" plus cog of group "MOL"'


sed -i '/^[#@]/d' distances_C118_MOL_rep_${rep}.xvg  
awk '{print $1, $4}' distances_C118_MOL_rep_${rep}.xvg >> dist_C118_MOL_clean_${rep}.xvg 

#for rep in 1 2 3 4 5; do ./sel_confs.py -f rep_${rep}/dist_C118_MOL_clean_${rep}.xvg -win 0.2 -low 0 -lar 3.5 > SELECTED_distances_rep${rep}.dat
#rep=4
#$grom trjconv -f concatenated_r4.xtc -s ${tpr} -n ${index} -o conf_55550.gro -dump 55550

cd ../
done