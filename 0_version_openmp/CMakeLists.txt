cmake_minimum_required(VERSION 2.8.8)

project(FWI C)

# set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${})

#--------------------------------------
# Configure high-level compiler flags
#--------------------------------------
set (architecture "generic" CACHE STRING "Architecture descriptor")

if ("${architecture}" STREQUAL "generic")
	add_definitions ("-DGENERIC")
	set (CMAKE_C_FLAGS "-O0")
elseif ("${architecture}" STREQUAL "KNL")
	add_definitions ("-DKNL")
	set (CMAKE_C_FLAGS "-xMIC-AVX512")
else()
	message (FATAL_ERROR "Architecture descriptor is needed!")
endif()

#--------------------------------------
#   Common compilation flags, if any.
#--------------------------------------
# set (CMAKE_C_FLAGS "-std=c99 -qopenmp -fma -align -finline-functions")
set (CMAKE_C_FLAGS "-g -std=c99 -fopenmp")
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
#       Configure statistics gathering
#--------------------------------------
option (collect-stats "Collect performance statistics" OFF)

if (collect-stats AND perform-io)
	message(STATUS "Kernel will collect performance statistics")
	add_definitions("-DCOLLECT_STATS")
elseif (collect-stats AND NOT perform-io)
	message(WARNING "Kernel cant collect performance statistics, perform-io flag is disabled")
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
#          Create fwi binary
#--------------------------------------
include_directories(
	"../common/"
	"./")

add_executable(fwi.bin
 	"./fwi_main.c"	
	"../common/fwi_constants.c"
	"../common/fwi_common.c"
	"../common/fwi_comms.c"
	"../common/fwi_kernel.c"
	"../common/fwi_sched.c"
	"./fwi_propagator.c")

target_link_libraries(fwi.bin "-lm")

install (TARGETS fwi.bin DESTINATION "${CMAKE_SOURCE_DIR}/")
#--------------------------------------
#    Create model generator utility
#--------------------------------------

include_directories(
	"../common/"
	"./")
add_executable(generateModel.bin
	"../common/fwi_constants.c"
	"../common/fwi_common.c"
	#"../common/fwi_comms.c"
	"./fwi_propagator.c"
	"../common/fwi_kernel.c"
	"../common/fwi_sched.c"
	"../common/fwi_generatemodel.c")

target_link_libraries(generateModel.bin "-lm")
set_target_properties(generateModel.bin PROPERTIES CMAKE_C_FLAGS "-Wall -O2")
install (TARGETS generateModel.bin DESTINATION "${CMAKE_SOURCE_DIR}/")

#--------------------------------------
#   Create model scheduler utility
#---------------------------------------
include_directories(
	"../common/"
	"./")
add_executable(generateSchedule.bin
	"../common/fwi_constants.c"
	"../common/fwi_common.c"
	#"../common/fwi_comms.c"
	"./fwi_propagator.c"
	"../common/fwi_kernel.c"
	"../common/fwi_sched.c"
	"../common/fwi_schedule.c")

target_link_libraries(generateSchedule.bin "-lm")
set_target_properties(generateSchedule.bin PROPERTIES CMAKE_C_FLAGS "-Wall -O2")
install (TARGETS generateSchedule.bin DESTINATION "${CMAKE_SOURCE_DIR}/")
