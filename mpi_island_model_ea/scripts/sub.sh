#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --exclude node052

module load openmpi

WORK="/home/morgaia/mpi_island_model_ea/mpi_island_model_ea"
cd $WORK

MODE="evo"
if [ $# -gt 0 ]; then
   echo "changing ea mode to $1"
   sed -i "s/^ea_mode:[0|1]$/ea_mode:${1}/" config.txt
   MODE="bench"
fi

if [[ "$MODE" == "bench" ]]; then
COSCRIPT="${WORK}/sub_${MODE}_${SLURM_NTASKS}.sh"
echo "#!/bin/bash" > $COSCRIPT
echo "#SBATCH -n $SLURM_NTASKS" >> $COSCRIPT
echo "#SBATCH --exclude node052" >> $COSCRIPT
echo "#SBATCH --chdir=$WORK" >> $COSCRIPT
echo "#SBATCH --nodelist $SLURM_NODELIST" >> $COSCRIPT
echo "cd $WORK" >> $COSCRIPT
echo "sed -i 's/^ea_mode:1/ea_mode:0/' ${WORK}/config.txt" >> $COSCRIPT
echo "module load openmpi" >> $COSCRIPT
echo "mpirun -n $SLURM_NTASKS -mca btl_tcp_if_include ib0 ${WORK}/island" >> $COSCRIPT
cd $WORK
chmod 750 $COSCRIPT
ssh morgaia@login sbatch -d afterok:$SLURM_JOBID $COSCRIPT
fi

mpirun -n $SLURM_NTASKS -mca btl_tcp_if_include ib0 ./island

exit 0
