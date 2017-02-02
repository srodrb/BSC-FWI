#!/usr/bin
rm fwi.log
./generateSchedule.bin fwi_params.txt.singlenode fwi_frequencies.txt
./generateModel.bin fwi_schedule.txt
./fwi.bin fwi_schedule.txt
