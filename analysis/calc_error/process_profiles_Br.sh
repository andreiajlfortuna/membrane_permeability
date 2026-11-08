#!/bin/bash -e

dir=/home/andreg/Documents/Cobimetinib/MD/POPC_CpHMD/us_CpHMD/EP

rm -f inser_aux* temp* aux_* new_aux* XBs_* inser_*

for accept in O12-O22 O33-O34
do
    for rep in 1 2 3
    do
        for umb in 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4
        do

            ####   Br1 ($3<3.50 && $4>140) // H15 ($3<2.50 && $4>150)
            awk 'NR>1 {if ($3<3.37 && $4>140) {print $3, $4, "1"} else {print $3, $4, "0"}}' ../XBs/MD_angle_dist_BR1_${accept}_rep_${rep}_1_dists_cap_${umb}.csv > aux_${accept}_r${rep}_${umb}

            #add the XBs=0 for umbrella 3.6 and 3.8 cuz there are no NCI's w/ the phospates in bulk water
            awk '{$3="0"}1' aux_${accept}_r${rep}_0.0 > aux_${accept}_r${rep}_3.6 #there is no problem for the dist ($1) and angle ($2) to be a copy, only the zero ($3) matters, which means no XB
           
           
            done
        done
    done

for accept in O33-O34 O12-O22
do
    for rep in 1 2 3
    do
        for umb in 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4
        do

            #match XBs and inser entries
            awk 'NR>1 {print $1, $2}' ../insertion_files/inser-Br1-${umb}-r${rep}_insertion.xvg | head -n -1 > temp1_${accept}_r${rep}_${umb}
            #awk '{if (NR % 10 == 1) {print $3}}' aux_${accept}_r${rep}_${umb} | head -n -1 > temp2_${accept}_r${rep}_${umb}
            awk '{print $3}' aux_${accept}_r${rep}_${umb} | head -n -1 > temp2_${accept}_r${rep}_${umb}
            paste temp1_${accept}_r${rep}_${umb} temp2_${accept}_r${rep}_${umb} > temp3_${accept}_r${rep}_${umb}

            awk '{print $1, $3}' temp3_${accept}_r${rep}_${umb} > XBs_${accept}_r${rep}_${umb}
            awk '{print $1, $2}' temp3_${accept}_r${rep}_${umb} > inser_aux_${accept}_r${rep}_${umb}
            
            awk '!seen[$1]++' inser_aux_${accept}_r${rep}_${umb} > inser_${accept}_r${rep}_${umb}

            #cat all umbrellas per rep
            cat XBs_${accept}_r${rep}_?.? >> XBs_${accept}_r${rep}.dat
            cat inser_${accept}_r${rep}_?.? >> inser_${accept}_r${rep}.dat

            #cat reps
            cat inser_${accept}_r${rep}.dat >> new_aux1_${accept}
            cat XBs_${accept}_r${rep}.dat >> new_aux2_${accept}
            paste new_aux1_${accept} new_aux2_${accept} >> aux_${accept}

            #remove MembIT empty data
            egrep -v "NaN" aux_${accept} > clean_BR_${accept}.dat


        done

###############################################################################################################

min=-20
max=20


#insertion min; bin increment; max
    for sli in $(seq $min 0.1 $max )
    do
        slicei=`echo $sli | sed "s/,/./g"`
        sliceo=`echo $slicei | awk '{print $1+1}'` #bin size
        slicem=`echo $slicei $sliceo | awk '{print ($1+$2+0.01)/2}'`

        awk -v slicei=$slicei -v sliceo=$sliceo -v m=$slicem '$2>=slicei && $2<=sliceo{s+=$4;n++}END{print m, s, n}' clean_BR_${accept}.dat >> profile_EP_${accept}_r${rep}_raw.xvg
        awk 'NF!=1{print $1, $2/$3, $3}' profile_EP_${accept}_r${rep}_raw.xvg > profile_EP_${accept}_r${rep}_proc.xvg

        #number of points min: cutoff to plot profile
        awk '$3>150' profile_EP_${accept}_r${rep}_proc.xvg > profile_EP_BR_${accept}_r${rep}_final.xvg
        done
    done
done

rm -f aux* temp* inser* XBs* profile_EP_O1A-O2A_r?_raw.xvg profile_EP_OB-OG_r?_raw.xvg profile_EP_O1A-O2A_r?_proc.xvg profile_EP_OB-OG_r?_proc.xvg clean_EP_*

echo "done"

exit 0

