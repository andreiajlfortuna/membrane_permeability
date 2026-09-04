#!/bin/sh

grom=/gromacs/gromacs-2020.6/bin/gmx
#grom=/gromacs/gromacs-2020.6-GPU/bin/gmx
export AMBERHOME=/home/pjcosta/bin/amber/amber16
#
CPUs=4
#
ndx=../01_min/index.ndx
top=$(ls ../01_min/*_EP_POPC.top)
prev=../01_min/min2.gro
curr=init1
#

######################################
# Generate tpr file for minimization #
######################################

${grom} grompp -f ${curr}.mdp -po ${curr}_out.mdp -c ${prev} -r ${prev} -p ${top} -pp TMP_processed.top -o ${curr}.tpr -n ${ndx} -maxwarn 1000

##################
# Run initiation #
##################

${grom} mdrun -nt ${CPUs} -pin on -s ${curr}.tpr -x ${curr}.xtc -c ${curr}.gro -e ${curr}.edr -g ${curr}.log -v -nice 19 -rcon 0

####################
# Extract pdb file #
####################

${grom} trjconv -f ${curr}.gro -s ${curr}.tpr -n ${ndx} -o ${curr}.pdb -pbc mol<<EOF
System
EOF

#
prev=./init1.gro
curr=init2
#

######################################
# Generate tpr file for minimization #
######################################

${grom} grompp -f ${curr}.mdp -po ${curr}_out.mdp -c ${prev} -r ${prev} -p ${top} -pp TMP_processed.top -o ${curr}.tpr -n ${ndx} -maxwarn 1000

##################
# Run initiation #
##################

${grom} mdrun -nt ${CPUs} -pin on -s ${curr}.tpr -x ${curr}.xtc -c ${curr}.gro -e ${curr}.edr -g ${curr}.log -v -nice 19 -rcon 0

####################
# Extract pdb file #
####################

${grom} trjconv -f ${curr}.gro -s ${curr}.tpr -n ${ndx} -o ${curr}.pdb -pbc mol<<EOF
System
EOF

#
prev=./init2.gro
curr=init3
#

######################################
# Generate tpr file for minimization #
######################################

${grom} grompp -f ${curr}.mdp -po ${curr}_out.mdp -c ${prev} -r ${prev} -p ${top} -pp TMP_processed.top -o ${curr}.tpr -n ${ndx} -maxwarn 1000

##################
# Run initiation #
##################

${grom} mdrun -nt ${CPUs} -pin on -s ${curr}.tpr -x ${curr}.xtc -c ${curr}.gro -e ${curr}.edr -g ${curr}.log -v -nice 19 -rcon 0

####################
# Extract pdb file #
####################

${grom} trjconv -f ${curr}.gro -s ${curr}.tpr -n ${ndx} -o ${curr}.pdb -pbc mol<<EOF
System
EOF

#
exit 0

