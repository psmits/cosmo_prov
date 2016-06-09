#!/bin/sh

### Set the name of the job, where jobname is a unique name for your job
#PBS -N coping

### Select the shell you would like the script to execute within
#PBS -S /bin/sh

### Inform the scheduler of the expected runtime, where walltime=HH:MM:SS
#PBS -l walltime=500:00:00

### Inform the scheduler of the number of CPU cores for your job.
### This example will allocate four cores on a single node.
#PBS -l nodes=1:ppn=4

### Inform the scheduler of the amount of memory you expect to use.
### Use units of 'b', 'kb', 'mb', or 'gb'
#PBS -l mem=32gb

### Set the destination for your program's output.
#PBS -o $HOME/myjob.out
#PBS -e $HOME/myjob.err


# Execute the program
bash /psmits/cosmo_prov/src/extinction_model_wrap.sh
