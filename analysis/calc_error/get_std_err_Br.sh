#!/bin/bash

for accept in O12-O22 O33-O34
do

cp ../profile_EP_BR_${accept}_final.xvg .

n=`cat profile_EP_BR_${accept}_final.xvg | wc -l`
rm -f aux*

    for i in $(seq 1 $n)
    do
        for rep in {1..3}
        do
            slab=`awk -v i=$i 'NR == i { print $1}' profile_EP_BR_${accept}_final.xvg` 

            value_rep=`awk -v slab="$slab" '{if ($1 == slab) print $2}' profile_EP_BR_${accept}_r${rep}_final.xvg`
            if [ -z "$value_rep" ]
            then 
            echo "0" >> aux1
            else 
            echo "$value_rep" >> aux1
            fi
        done

    ./standard-error aux1 >> aux2.dat
    rm aux1*

    done

awk '{print $4}' aux2.dat > std_error_${accept}.dat

paste profile_EP_BR_${accept}_final.xvg std_error_${accept}.dat > profile_w_error_BR_${accept}.xvg

rm aux2.dat

done

exit 0

