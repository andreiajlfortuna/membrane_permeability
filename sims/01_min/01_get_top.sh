#!/bin/sh


#
#grom=/gromacs/gromacs-5.1.2/bin/gmx
grom=/gromacs/gromacs-2020.6-GPU/bin/gmx
export AMBERHOME=/home/pjcosta/bin/amber/amber16

######################################
# create top and gro for ligand
######################################

mol2=$(ls *.mol2)
name=$(basename ${mol2} -opt-RESP_EP.mol2 )
lig=$(echo "${name}_EP")


$AMBERHOME/bin/parmchk -i ${mol2} -f mol2 -o ${lig}.frcmod

frcmod=$(ls *.frcmod)

cat << EOF > leap.in
source leaprc.gaff #Source leaprc file for gaff
mods = loadamberparams ${frcmod} #Source the missing parameters

mol = loadmol2 ${mol2}
saveAmberParm mol ${lig}.top ${lig}.crd
quit
EOF

$AMBERHOME/bin/tleap -f leap.in

cp /home/pjcosta/bin/acpype.py . 

python3 ./acpype.py -p ${lig}.top -x ${lig}.crd
mv MOL_GMX.gro ${lig}_GMX.gro
mv MOL_GMX.top ${lig}_GMX.top

cp ${lig}_GMX.top ${lig}_GMX.top.bak

#######################################
## find halogens
#######################################

mol2=`ls *.mol2`

ncl=`cat ${mol2} | awk '$6 == "cl"' | wc -l`
nbr=`cat ${mol2} | awk '$6 == "br"' | wc -l`
ni=` cat ${mol2} | awk '$6 == "i"'  | wc -l`
totX=`echo ${ncl}+${nbr}+${ni} | bc -l`

if [ ${ncl} -gt 0 -a ${nbr} -eq 0 -a ${ni} -eq 0 ]
then
echo "Molecule ${bname} has ${ncl} Chlorine atoms." >&2
Xtype1=CL 
dEP1=1.948
nX1=${ncl}
nX2=0
dEP2=0
elif [ ${ncl} -eq 0 -a ${nbr} -gt 0 -a ${ni} -eq 0 ]
then
echo "Molecule ${bname} has ${nbr} Bromine atoms." >&2
Xtype1=BR
dEP1=2.02
nX1=${nbr}
nX2=0
dEP2=0
elif [ ${ncl} -eq 0 -a ${nbr} -eq 0 -a ${ni} -gt 0 ]
then
echo "Molecule ${bname} has ${ni} Iodine atoms." >&2
Xtype1=I
nX1=${ni}
dEP1=2.15
nX2=0
dEP2=0
elif [ ${ncl} -gt 0 -a ${nbr} -gt 0 -a ${ni} -eq 0 ]
then
echo "Molecule ${bname} has ${ncl} Chlorine atoms and ${nbr} Bromine atoms." >&2
Xtype1=CL
nX1=${ncl}
dEP1=1.948
Xtype2=BR
nX2=${nbr}
dEP2=2.02
elif [ ${ncl} -eq 0 -a ${nbr} -gt 0 -a ${ni} -gt 0 ]
then
echo "Molecule ${bname} has ${nbr} bromine atoms and ${ni} Iodine atoms." >&2
Xtype1=BR
nX1=${nbr}
dEP1=2.02
Xtype2=I
nX2=${ni}
dEP2=2.15
elif [ ${ncl} -gt 0 -a ${nbr} -eq 0 -a ${ni} -gt 0 ]
then
Xtype1=CL
nX1=${ncl}
dEP1=1.948
Xtype2=I
nX2=${ni}
dEP2=2.15
echo "Molecule ${bname} has ${ncl} Chlorine atoms and ${ni} Iodine atoms." >&2
fi
echo "totX is ${totX}"


#######################################################
# Remove bonds, angles, dihedrals for EPs
#######################################################
sed -i '/ - EP/d;/EP-/d;/-    EP/d;/EP -/d' ${lig}_GMX.top

########################################################
# In GMX, EPs are represented as virtual sites
# Add virtual site for each halogen
#########################################################

# Identifies the halogens


for nX in $(sed -n '/bonds/,/pairs/p' ${lig}_GMX.top.bak | awk -v Xtype1=${Xtype1} '{if ($7~Xtype1) print $1}')
do

nC=`sed -n '/bonds/,/pairs/p' ${lig}_GMX.top | awk -v nX=$nX '{if ($1 == nX) print $2; else if ($2 == nX) print $1}'`
nEP=`sed -n '/bonds/,/pairs/p' ${lig}_GMX.top.bak | grep EP | awk -v nX=$nX '{if ($1 == nX) print $2; else if ($2 == nX) print $1}'`
Cdist=`sed -n '/bonds/,/pairs/p' ${lig}_GMX.top | awk -v nX=$nX '{if ($1 == nX || $2 == nX) print $4}'`
epdist=`awk -v xep=${dEP1} 'BEGIN {print -xep/10}'`

echo ${nEP} ${nX} ${nC} 2 ${epdist} >>  VSITE.${lig}
done 

if [ ${nX2} -gt 0 ]
then

for nX in $(sed -n '/bonds/,/pairs/p' ${lig}_GMX.top.bak | awk -v Xtype2=${Xtype2} '{if ($7~Xtype2) print $1}')
do

nC=`sed -n '/bonds/,/pairs/p' ${lig}_GMX.top | awk -v nX=$nX '{if ($1 == nX) print $2; else if ($2 == nX) print $1}'`
nEP=`sed -n '/bonds/,/pairs/p' ${lig}_GMX.top.bak | grep EP | awk -v nX=$nX '{if ($1 == nX) print $2; else if ($2 == nX) print $1}'`
Cdist=`sed -n '/bonds/,/pairs/p' ${lig}_GMX.top | awk -v nX=$nX '{if ($1 == nX || $2 == nX) print $4}'`
epdist=`awk -v xep=${dEP2} 'BEGIN {print -xep/10}'`
 
echo ${nEP} ${nX} ${nC} 2 ${epdist} >>  VSITE.${lig}
done 
fi
# write virtual site block
# we will futuraly use the new 2fd type

sed -i '1 i\[ virtual_sites2 ]' VSITE.${lig}
sed -i '1 i\ ' VSITE.${lig}

sed -i "/qtot 0.000/r VSITE.${lig}" ${lig}_GMX.top
sed -i "/qtot -0.000/r VSITE.${lig}" ${lig}_GMX.top
rm -rf VSITE.${lig}
rm -rf em.mdp md.mdp

#create ligand.itp file and gaff_atom_type.itp file
cp ${lig}_GMX.top ${lig}.itp

awk '/atomtypes/,/moleculetype/' ${lig}_GMX.top | head -n -1 > ${lig}_gaff_atom_types.itp

n=$(awk '/defaults/,/moleculetype/' ${lig}_GMX.top  | wc -l)
nline=$((n+1))
sed -i "3,${nline}d" ${lig}.itp 
NLf=`wc -l ${lig}.itp | awk '{print $1}'`
NLi=$((NLf-6))
sed -i "${NLi},${NLf}d" ${lig}.itp

echo "" >> ${lig}.itp 
echo "; Include Position restraint file" >> ${lig}.itp
echo "#ifdef POSRES_MOL" >> ${lig}.itp
echo "#include \"posres_${lig}.itp\"" >> ${lig}.itp
echo "#endif" >> ${lig}.itp
echo "" >> ${lig}.itp

######################################
# create index
######################################

/bin/rm -f index.ndx

gro=${lig}_POPC.gro

${grom} make_ndx -f ${gro} <<EOF
"PA" | "PC" | "OL"
name 9 POPC
"POPC" | "MOL"
q
EOF

#
echo "[ P ]" >> index.ndx
grep ' P' ${gro} | awk '{print substr($0,16,5)}' >> index.ndx
 
echo "[ C118 ]" >> index.ndx
grep ' C118' ${gro} | awk '{print substr($0,16,5)}' >> index.ndx

echo "[ One-tail ]" >> index.ndx
echo "131" >> index.ndx

#######################################
# change system top
#######################################

cp POPC_template.top ${lig}_POPC.top

sed -i "s/YYYY/${lig}/g" ${lig}_POPC.top
 

### correct the posre file with the atom numbers if needed
