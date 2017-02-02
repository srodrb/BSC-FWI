#!/usr/bin
module load extrae/3.4.1

rm fwi.log
./generateSchedule.bin fwi_params.txt.singlenode fwi_frequencies.txt
./generateModel.bin fwi_schedule.txt

export OMP_NUM_THREADS=8
source ../Scripts/extrae_configs/trace.omp.sh
./fwi.bin fwi_schedule.txt
