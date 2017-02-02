cmake_minimum_required(VERSION 2.8.8)

project(FWI C)

#--------------------------------------
# Find and configure MPI compiler
#--------------------------------------
find_package(MPI REQUIRED QUIET)

message(STATUS "--------------------------------------")
message(STATUS "----[MPI ENABLED] setting:            ")
message(STATUS "---- MPI_C_COMPILER:      ${MPI_C_COMPILER}")
message(STATUS "---- MPI_C_COMPILE_FLAGS: ${MPI_C_COMPILE_FLAGS}")
message(STATUS "---- MPI_C_INCLUDE_PATH:  ${MPI_C_INCLUDE_PATH}")
message(STATUS "---- MPI_C_LINK_FLAGS:    ${MPI_C_LINK_FLAGS}")
message(STATUS "---- MPI_C_LIBRARIES:     ${MPI_C_LIBRARIES}")
message(STATUS "--------------------------------------")

include_directories(${MPI_C_INCLUDE_PATH})
add_definitions("-DDISTRIBUTED_MEMORY_IMPLEMENTATION")


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
	message (FATAL_ERROR "Invalid architecture descriptor!")
endif()

#--------------------------------------
#   Common compilation flags, if any.
#--------------------------------------
# set (CMAKE_C_FLAGS "-std=c99 -qopenmp -fma -align -finline-functions")
set (CMAKE_C_FLAGS "-g -std=c99 -qopenmp")


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
	message(STATUS "GPFS system will be used for local I/O")
endif()


#--------------------------------------
#          Create fwi binary
#--------------------------------------
include_directories(
	"../common/"
	"./")

add_executable(fwi.intel64
	"../common/fwi_constants.c"
	"../common/fwi_common.c"
	"../common/fwi_comms.c"
	"./fwi_propagator.c"
	"../common/fwi_kernel.c"
	"../common/fwi_sched.c"
 	"./fwi_main.c")

target_link_libraries(fwi.intel64 "-lm")
install (TARGETS fwi.intel64 DESTINATION "${CMAKE_SOURCE_DIR}/")
#--------------------------------------
#    Create model generator utility
#--------------------------------------
# include_directories(
# 	"../common/"
# 	"./")
# add_executable(generateModel.bin
# 	"../common/fwi_constants.c"
# 	"../common/fwi_common.c"
# 	"../common/fwi_kernel.c"
# 	"../common/fwi_sched.c"
# 	"./fwi_propagator.c"
# 	"../common/fwi_generatemodel.c")
# 
# set_target_properties(generateModel.bin PROPERTIES CMAKE_C_FLAGS "-Wall -O2")
# install (TARGETS generateModel.bin DESTINATION "${CMAKE_SOURCE_DIR}/")
# 
# #--------------------------------------
# #   Create model scheduler utility
# #---------------------------------------
# include_directories(
# 	"../common/"
# 	"./")
# add_executable(generateSchedule.bin
# 	"../common/fwi_constants.c"
# 	"../common/fwi_common.c"
# 	"../common/fwi_kernel.c"
# 	"../common/fwi_sched.c"
# 	"./fwi_propagator.c"
# 	"../common/fwi_schedule.c")
# 
# set_target_properties(generateSchedule.bin PROPERTIES CMAKE_C_FLAGS "-Wall -O2")
# install (TARGETS generateSchedule.bin DESTINATION "${CMAKE_SOURCE_DIR}/")



