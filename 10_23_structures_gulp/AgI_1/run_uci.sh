#!/bin/bash
#$ -N forstgulp
#$ -q pub8i
#$ -m beas

# Notes:
# -----------------
# For SINGLE node MPI jobs, use the parallel envrironment "one-node-mpi"
# and not the "mpi".
#
# Use "-R y" to request job reservation otherwise single 1-core jobs
# may prevent this multicore MPI job from running.   This is called
# job starvation.


# Module load OpenMPI.
module load openmpi-1.8.3/gcc-4.8.3



echo "Our MPI job has the following node/cores allocated:"  
cat $PE_HOSTFILE

# Run the program with the ouput going to file out
# No need to specify the number of $CORES to run with, but it is here
# for documentaion.
 


/data/apps/user_contributed_software/szare/gulp/4.2/gulp-4.2.0/Src/gulp <mgo.gin> output.txt
