#!/bin/bash 

#PBS -N MPI
#PBS -o stencil_output
#PBS -q teaching
#PBS -l nodes=1:ppn=16,walltime=00:10:00

#! Mail to user if job aborts
#PBS -m a

#! application name
application="stencil"

#! Run options for the application
options="1024 1024 100"

###############################################################
### You should not have to change anything below this line ####
###############################################################

#! change the working directory (default is home directory)

cd $PBS_O_WORKDIR

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This jobs runs on the following machines:
echo `cat $PBS_NODEFZZILE | uniq`
 
#! Create a machine file for MPI
cat $PBS_NODEFILE > machine.file.$PBS_JOBID

numnodes=`wc $PBS_NODEFILE | awk '{ print $1 }'`

#! Run the parallel MPI executable (nodes*ppn)
mpirun -np $numnodes -machinefile machine.file.$PBS_JOBID $application $options