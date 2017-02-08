#!/bin/bash
export FWIDIR=$PWD/..
./generateSchedule.bin fwi_params.txt fwi_frequencies.txt
./generateModel.bin fwi_schedule.txt
