cmake_minimum_required(VERSION 2.8.8)

project(FWI C)

#--------------------------------------
#   Common compilation flags, if any.
#--------------------------------------
set (CMAKE_C_FLAGS "-g -std=c99 -O2")
add_definitions ("-DSHARED_MEMORY_IMPLEMENTATION")

#--------------------------------------
#   Configure release
#--------------------------------------
set (debug "Debug mode enabled" OFF)

if (debug)
	message(STATUS "Debug mode enabled")
	add_definitions ("-DDEBUG")
endif()

#--------------------------------------
#       Configure Input/Output
#--------------------------------------
option (perform-io "Perform local I/O" ON)

if (perform-io)
	message(STATUS "Input/Output is enabled.")
else()
	message(STATUS "Input/Output is disabled.")
	add_definitions("-DDO_NOT_PERFORM_IO")
endif()

#--------------------------------------
#        Use NMVe local storage
#--------------------------------------
option (use-nmve "Use NMVe devices for local I/O" OFF)

if (use-nmve)
	add_definitions("-DUSE_NMVE")
	message(STATUS "NMVe devices will be used for local I/O")
else()
	message(STATUS "GPFS will be used for local I/O")
endif()

#--------------------------------------
#    Create model generator utility
#--------------------------------------
include_directories(
	"../common/"
	"../0_version_openmp/"
	"./")
add_executable(generateModel.bin
	"../common/fwi_constants.c"
	"../common/fwi_common.c"
	#"../common/fwi_comms.c"
	"../0_version_openmp/fwi_propagator.c"
	"../common/fwi_kernel.c"
	"../common/fwi_sched.c"
	"../common/fwi_generatemodel.c")

target_link_libraries(generateModel.bin "-lm")
install (TARGETS generateModel.bin DESTINATION "${CMAKE_SOURCE_DIR}/")

#--------------------------------------
#   Create model scheduler utility
#---------------------------------------
include_directories(
	"../common/"
	"../0_version_openmp/"
	"./")
add_executable(generateSchedule.bin
	"../common/fwi_constants.c"
	"../common/fwi_common.c"
	#"../common/fwi_comms.c"
	"../0_version_openmp/fwi_propagator.c"
	"../common/fwi_kernel.c"
	"../common/fwi_sched.c"
	"../common/fwi_schedule.c")

target_link_libraries(generateSchedule.bin "-lm")
install (TARGETS generateSchedule.bin DESTINATION "${CMAKE_SOURCE_DIR}/")