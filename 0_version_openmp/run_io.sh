#!/usr/bin
rm fwi.log
../utils/generateSchedule.bin fwi_params.txt.mpi_local fwi_frequencies.txt
../utils/generateModel.bin fwi_schedule.txt
# gdb --args ./fwi.bin fwi_schedule.txt
./fwi.bin fwi_schedule.txt
