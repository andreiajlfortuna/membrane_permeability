#! /bin/bash -e
#
## Working Dir:
Dir=`pwd`
#
# Source the .dat file that has all the settings
#runname=`basename $0`; source ${runname%.*}.dat
#
# Override Partition chosen in case it's a GPU job
#
## Beginning of QEXE file
 
 
source run_init.dat

cat <<EOF >$Name.slurm
#! /bin/bash -e
#
#. /etc/profile
#. ~/.bashrc
#. ~/.bash_profile
export EMAIL=$Email
InitDate=\`date\`
#
LocalDir=\`pwd\`
#


 
CPUs=8
#
ndx=../01_min/index.ndx
top=../01_min/*_POPC_NOCHARGES.top
 
 
 
curr=init3
#

######################################
# Generate tpr file for minimization #
######################################

${grom} grompp -f \${curr}.mdp -po \${curr}_out.mdp -c conf_104040.gro -r conf_104040.gro -p \${top} -pp TMP_processed.top -o \${curr}.tpr -n \${ndx} -maxwarn 1000

##################
# Run initiation #
##################
# Using CPU
#${grom} mdrun -nt \${CPUs} -pin on -s \${curr}.tpr -x \${curr}.xtc -c \${curr}.gro -e \${curr}.edr -g \${curr}.log -v -nice 19 -rcon 0

# Using GPU

OMP_NUM_THREADS=\${CPUs}

if [[ \`nvidia-smi -L | wc -l\` > 1 ]]
    then
       if [[ \`ps -aux | awk '{for(i=1;i<=NF;i++){ if(\$i=="-gpu_id"){print 1-\$(i+1)}}}' | wc -l\` > 1 ]]
       then
              echo "Error: no GPU is available for GPU job" >> ${Dir}/\${j}.info
              exit 1
       elif [[ \`ps -aux | awk '{for(i=1;i<=NF;i++){ if(\$i=="-gpu_id"){print 1-\$(i+1)}}}' | wc -l\` < 1 ]]
       then
	       GPUid=0
	else
              GPUid=\`ps -aux | awk '{for(i=1;i<=NF;i++){ if(\$i=="-gpu_id"){print 1-\$(i+1)}}}'\`
       fi

    else
        GPUid=0
    fi

${grom} mdrun -ntomp \${CPUs} -ntmpi 1 -gpu_id \$GPUid -pin on -s \${curr}.tpr -x \${curr}.xtc -c \${curr}.gro -e \${curr}.edr -g \${curr}.log -v -nice 19 -rcon 0


EOF

chmod +x $Name.slurm

requeue=1

if [ $requeue == 1 ]
   then 
       sbatch --requeue -p $Partition -N 1 -n $ncpus -o $Name.sout -e $Name.serr $Name.slurm
else
       sbatch -p $Partition -N 1 -n $ncpus -o $Name.sout -e $Name.serr $Name.slurm
fi
echo ""
echo "Job submitted to Partition(s): $Partition with $ncpus Processors"
#
## End of Script
#
