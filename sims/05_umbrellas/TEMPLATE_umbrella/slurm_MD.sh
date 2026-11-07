#! /bin/bash -e
#
## Working Dir:
Dir=`pwd`
#
# Source the .dat file that has all the settings
runname=`basename $0`; source ${runname%.*}.dat
#
# Override Partition chosen in case it's a GPU job
if [ `echo $grom|grep "GPU"` ]; then
    if [ `echo $Partition|grep -v "GPU"` ]; then
        echo "GMX will run with GPU, but no GPU partition was selected"
        exit 1
    fi
fi
#
## Beginning of executable file
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
ncpus=$ncpus
#
# Finds current cycle
i=1
j=001
while [ -f \${j}.gro ] ; do
        i=\$((i+1))
        j=\`printf "%03d\n" \${i}\`
done
#
# Get the first .gro
if [ \${i} -eq 1 ]; then rm -f 000.gro; ln -s $gro 000.gro; fi
#
k=\$((i-1))
l=\`printf "%03d\n" \${k}\` 
#
# Info about the Job:
echo "Job executed in Machine \$HOSTNAME" > \${j}.info
echo "Job executed with \$ncpus processors" >> \${j}.info
echo "Job executed in DIR: /tmp/\${USER}_MD\$$ " >> \${j}.info
echo "" >> \${j}.info
echo -e "Job started on: " >> \${j}.info
date >> \${j}.info
#
# Run job locally:
# Copy important files for local directory
mkdir -p /tmp/\${USER}_MD\$$
cp -f \${l}.gro /tmp/\${USER}_MD\$$/\${l}.gro
tinit=\`awk -v K=\$k '\$1=="nsteps"{a=\$3};\$1=="dt"{b=\$3}END{print a*b*K}' $mdp\`
sed "s/tinit\s*=\s*0.0/tinit               =  \${tinit}.0/g" $mdp > /tmp/\${USER}_MD\$$/\${j}.mdp
#
# Enter local DIR
cd /tmp/\${USER}_MD\$$/
#
# Run MD segment:
#
## Make the tpr file
$grom grompp -f \${j}.mdp \
    -po \${j}_out.mdp \
    -c \${l} \
    -n ${Dir}/${index} \
    -p ${Dir}/${top} \
    -pp TMP_processed.top \
    -o \${j}.tpr \
    -maxwarn 1000 > ${Dir}/${Name}.out 2>${Dir}/${Name}.err
#
# Check if its CPU or GPU
if [ ! \`echo $grom|grep "GPU"\` ]; then
   #
   # Run MD
   $grom mdrun -nt \$ncpus -pin auto\
        -s \${j}.tpr \
        -x \${j}.xtc \
        -c \${j}.gro \
        -e \${j}.edr \
        -g \${j}.log \
        -rcon $Rcon \
        -nice 19 >> ${Dir}/${Name}.out 2>>${Dir}/${Name}.err
    
else
    #
    OMP_NUM_THREADS=$ncpus
    #
    # Check which GPU should be used
    #
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
    #
    # Run MD
    $grom mdrun -ntomp $ncpus -ntmpi 1 -pin auto -gpu_id \$GPUid \
        -s \${j}.tpr \
        -x \${j}.xtc \
        -c \${j}.gro \
        -e \${j}.edr \
        -g \${j}.log \
        -rcon $Rcon \
        -nice 19 >> ${Dir}/${Name}.out 2>>${Dir}/${Name}.err
    #
fi
#
# Copy the Segment to remote DIR and check if everything is OK.
cp -df \${j}*.{tpr,gro,xtc,edr,log} ${Dir}
# Copy pull/US related files
if [ -f pullx.xvg ]; then for pull in x f; do cp -f pull\${pull}.xvg $Dir/\${j}\${pull}.xvg ; done; fi 
if (for f in \${j}*.{tpr,gro,xtc,edr,log}; do diff \$f ${Dir}/\$f; done); then
cd ${Dir}
gzip -9 \${j}*.{tpr,edr,log}
rm -rf /tmp/\${USER}_MD\$$
else
echo "Error in file copy... please check local files" >> ${Dir}/\${j}.info
exit 1
fi
#
# Info on the end of segment:
echo "" >> \${j}.info
echo -e "Job finished on: " >> \${j}.info
date >> \${j}.info
#
#
# Send email 
if [[ \$EMAIL != "" ]]; then
echo -e "Job $Name just finished on \$HOSTNAME with $ncpus CPU cores\n--Running since: \$InitDate \n--Finished on:   \`date\`" | mutt -s "BioISI Cluster Report: $Name" \$EMAIL
fi
#
if [ \${i} -lt ${Segments} ] # of cycles 
then
    ./`basename $0`
fi
#
exit 0
    
EOF
#
chmod +x $Name.slurm
#
if [ $requeue == 1 ]
   then 
       sbatch --requeue -p $Partition -N 1 -n $ncpus --mem=1 -o $Name.sout -e $Name.serr $Name.slurm
else
       sbatch -p $Partition -N 1 -n $ncpus --mem=1 -o $Name.sout -e $Name.serr $Name.slurm
fi
echo ""
echo "Job submitted to Partition(s): $Partition with $ncpus Processors"
#
## End of Script
#
