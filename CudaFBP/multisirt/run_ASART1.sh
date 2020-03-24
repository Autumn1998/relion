#!/bin/bash

#$ -V
#$ -cwd

### Combine output and error messages
#$ -j y
#$ -o testing.outerr.$JOB_ID
##$ -e testing.outerr.$JOB_ID
#$ -q all.q


# Requests  number of cpus and environment
#$ -pe orte 24
#$ -S /bin/bash


 
export LD_LIBRARY_PATH=/opt/openmpi/gnu/mx/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/cuda/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/x5wan/software/OpenCV-2.2.0/lib:$MPI/lib:$LD_LIBRARY_PATH


#
# Execute the job
#
MPIRUN=/opt/openmpi/gnu/mx/bin/mpirun
time $MPIRUN --mca btl_tcp_if_include eth0 -np $NSLOTS ./main Ani12Gr3Syn1.st Ani12Gr3Syn1_SIRT1-0.2.mrc Ani12Gr3Syn1_fbp Ani12Gr3Syn1.txbr SIRT -n 1 -t 0.2
#time $MPI/bin/mpirun -np $NSLOTS ./main Ani12Gr3Syn1.st Ani12Gr3Syn1_ASART3-0.2.mrc Ani12Gr3Syn1_fbp Ani12Gr3Syn1.txbr ASART -n 3 -t 0.2
exit

