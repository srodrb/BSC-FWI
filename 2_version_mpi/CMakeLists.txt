cmake_minimum_required(VERSION 2.8.8)

project(FWI C)

#--------------------------------------
# Find and configure MPI compiler
#--------------------------------------
find_package(MPI REQUIRED QUIET)

message(STATUS "-------------------------------------------------")
message(STATUS "---- MPI setting:                                ")
message(STATUS "---- MPI_C_COMPILER:      ${MPI_C_COMPILER}      ")
message(STATUS "---- MPI_C_COMPILE_FLAGS: ${MPI_C_COMPILE_FLAGS} ")
message(STATUS "---- MPI_C_INCLUDE_PATH:  ${MPI_C_INCLUDE_PATH}  ")
message(STATUS "---- MPI_C_LINK_FLAGS:    ${MPI_C_LINK_FLAGS}    ")
message(STATUS "---- MPI_C_LIBRARIES:     ${MPI_C_LIBRARIES}     ")
message(STATUS "-------------------------------------------------")

include_directories(${MPI_C_INCLUDE_PATH})
add_definitions("-DDISTRIBUTED_MEMORY_IMPLEMENTATION")

#--------------------------------------
# Configure high-level compiler flags
#--------------------------------------
set (architecture "generic" CACHE STRING "Architecture descriptor")

if ("${architecture}" STREQUAL "generic")
	add_definitions ("-DGENERIC")
	set (CMAKE_C_FLAGS "-xHost")
elseif ("${architecture}" STREQUAL "KNL")
	add_definitions ("-DKNL")
	set (CMAKE_C_FLAGS "-xMIC-AVX512")
else()
	message (FATAL_ERROR "Invalid architecture descriptor!")
endif()

message(STATUS "Compiler Id ${CMAKE_C_COMPILER_ID}")

#---------------------------------------------------
# Configure compiler flags for Intel/GCC compilers
#--------------------------------------------------
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -qopenmp -fma -align -finline-functions")
else()
	set (CMAKE_C_FLAGS "-O3 -std=gnu99 -fopenmp")
endif()

#--------------------------------------
#   Configure release
#--------------------------------------
option (debug "Debug mode enabled" OFF)

if (debug)
	message(STATUS "Debug mode enabled")
	set (CMAKE_C_FLAGS "-O0 -Wunused-result")
	add_definitions ("-DDEBUG")
else()
	message(STATUS "Debug mode disabled")
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
	message(WARNING "Kernel cant collect all performance statistics, perform-io flag is disabled")
	add_definitions("-DCOLLECT_STATS")
endif()


#--------------------------------------
#        Use NMVe local storage
#--------------------------------------
option (use-nvme "Use NVMe devices for local I/O" OFF)

if (use-nvme)
	add_definitions("-DUSE_NVME")
	message(STATUS "NMVe devices will be used for local I/O")
else()
	message(STATUS "Default global file system will be used for local I/O")
endif()

#--------------------------------------
#        Use NMVe local storage
#--------------------------------------
message(STATUS "----------------------------------------")
message(STATUS "CMAKE_C_COMPILER:               ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_C_FLAGS:                  ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_INCLUDE_PATH:           ${CMAKE_C_INCLUDE_PATH}")
message(STATUS "CMAKE_C_LIBRARIES:              ${CMAKE_C_LIBRARIES}")
message(STATUS "----------------------------------------")


#--------------------------------------
#          Create fwi binary
#--------------------------------------
include_directories(
	"../common/"
	"./")

add_executable(fwi.intel64
	"../common/fwi_constants.c"
	"../common/fwi_common.c"
	"./fwi_propagator.c"
	"../common/fwi_kernel.c"
	"../common/fwi_sched.c"
 	"./fwi_main.c")

target_link_libraries(fwi.intel64 "-lm")
install (TARGETS fwi.intel64 DESTINATION "${CMAKE_SOURCE_DIR}/")
