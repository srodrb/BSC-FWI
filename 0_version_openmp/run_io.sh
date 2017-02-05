#!/usr/bin
rm fwi.log
./generateSchedule.bin fwi_params.txt.mpi_local fwi_frequencies.txt
./generateModel.bin fwi_schedule.txt
# gdb --args ./fwi.bin fwi_schedule.txt
# ./fwi.bin fwi_schedule.txt
