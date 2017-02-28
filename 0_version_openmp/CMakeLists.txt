cmake_minimum_required(VERSION 2.8.8)
project(FWI C)

#--------------------------------------
# Configure high-level compiler flags
#--------------------------------------
set (architecture "generic" CACHE STRING "Architecture descriptor")

if ("${architecture}" STREQUAL "generic")
	add_definitions ("-DGENERIC")
	set (CMAKE_C_FLAGS "-xHost -std=c99 -qopenmp -restrict -align -finline-functions -fma -prec-div -prec-sqrt -fp-model=precise")
elseif ("${architecture}" STREQUAL "KNL")
	message(STATUS "Configuring for KNL")
	add_definitions ("-DKNL")
	set (CMAKE_C_FLAGS "-O3 -std=c99 -qopenmp -xMIC-AVX512 -restrict -align -fma -finline-functions")
else()
	message (FATAL_ERROR "Architecture descriptor is needed!")
endif()

message(STATUS "CMAKE_C_FLAGS are ${CMAKE_C_FLAGS}")

#--------------------------------------
#   Common compilation flags, if any.
#--------------------------------------
# if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
# 	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -qopenmp -fma -align -finline-functions")
# else()
# 	set (CMAKE_C_FLAGS "-O3 -g -std=gnu99 -fopenmp")
# endif()

add_definitions ("-DSHARED_MEMORY_IMPLEMENTATION")

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
#--------------------------------------option (use-nmve "Use NMVe devices for local I/O" OFF)
option (use-nvme "Use NVMe devices for local I/O" OFF)

if (use-nvme)
	add_definitions("-DUSE_NVME")
	message(STATUS "NVMe devices will be used for local I/O")
else()
	message(STATUS "Default global file system will be used for local I/O")
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
