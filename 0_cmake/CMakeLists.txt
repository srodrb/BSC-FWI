cmake_minimum_required(VERSION 2.8.8)

project(FWI C)

# set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${})

#--------------------------------------
# Configure high-level compiler flags
#--------------------------------------
set (architecture "generic" CACHE STRING "Architecture descriptor")

if ("${architecture}" STREQUAL "generic")
	add_definitions ("-DGENERIC")
	set (CMAKE_C_FLAGS "-O2")
elseif ("${architecture}" STREQUAL "KNL")
	add_definitions ("-DKNL")
	set (CMAKE_C_FLAGS "-xMIC-AVX512")
else()
	message (FATAL_ERROR "Architecture descriptor is needed!")
endif()

#--------------------------------------
#   Common compilation flags, if any.
#--------------------------------------
set (CMAKE_C_FLAGS "-std=c99 -qopenmp -fma -align -finline-functions")


#--------------------------------------
#       Configure Input/Output
#--------------------------------------
option (perform-io "Perform local I/O" ON)

if (perform-io)
	add_definitions("-DLOG_IO_STATS")
	message(STATUS "Real Input/Output is enabled. Statistics will be collected")
else()
	add_definitions("-DDO_NOT_PERFORM_IO")
	message(STATUS "Input/Output is disabled.")
endif()

#--------------------------------------
#        Use NMVe local storage
#--------------------------------------
option (use-nmve "Use NMVe devices for local I/O" OFF)

if (use-nmve)
	add_definitions("-DUSE_NMVE")
	message(STATUS "NMVe devices will be used for local I/O")
else()
	message(STATUS "Global file system will be used for local I/O")
endif()


#--------------------------------------
#          Create fwi binary
#--------------------------------------
include_directories(
	"../common/"
	"./")

add_executable(fwi.bin
 	"./fwi_main.c"	
	"../common/fwi_common.c"
	"../common/fwi_kernel.c"
	"./fwi_propagator.c")

#--------------------------------------
#    Create model generator utility
#--------------------------------------
include_directories(
	"../common/"
	"./")
add_executable(ModelGenerator
	"../common/fwi_common.c"
	"../common/fwi_kernel.c"
	"./fwi_propagator.c"
	"../common/fwi_generatemodel.c")

#--------------------------------------
#   Create model scheduler utility
#---------------------------------------
# include_directories(
# 	"../common/"
# 	"./")
# add_executable(Scheduler
# 	"../common/fwi_common.c"
# 	"../common/fwi_kernel.c"
# 	"./fwi_propagator.c"
# 	"../common/fwi_schedule.c")






