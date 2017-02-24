cmake_minimum_required(VERSION 2.8.8)

project(FWI C)


#--------------------------------------
# Configure high-level compiler flags
#--------------------------------------
set (architecture "generic" CACHE STRING "Architecture descriptor")
message(STATUS "Architecture is set to ${architecture}")

if ("${architecture}" STREQUAL "generic")
	add_definitions ("-DGENERIC")
	set (BACKEND_COMPILER_FLAGS "--Wn,-xHost")
elseif ("${architecture}" STREQUAL "KNL")
	add_definitions ("-DKNL")
	set (BACKEND_COMPILER_FLAGS "--Wn,-xMIC-AVX512")
else()
	message (FATAL_ERROR "Invalid architecture descriptor!")
endif()

add_definitions ("-DDISTRIBUTED_MEMORY_IMPLEMENTATION")

#--------------------------------------
#   Common compilation flags, if any.
#--------------------------------------
set (CMAKE_C_FLAGS "--ompss --no-copy-deps -g")

#---------------------------------------------------
# Configure compiler flags for Intel/GCC compilers
#--------------------------------------------------
message(STATUS "Compiler Id ${CMAKE_C_COMPILER_ID}")

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BACKEND_COMPILER_FLAGS} --Wn,-fma --Wn,-align --Wn,-finline-functions")
else()
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BACKEND_COMPILER_FLAGS} -std=gnu99")
endif()

message(STATUS "Compiler CFLAGS are: ${CMAKE_C_FLAGS}")


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
#   Configure release
#--------------------------------------
option (debug "Debug mode enabled" OFF)

if (debug)
	message(STATUS "Debug mode enabled")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wunused-result")
	add_definitions ("-DDEBUG")
else()
	message(STATUS "Debug mode disabled")
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
	"../common/fwi_kernel.c"
	"../common/fwi_sched.c"
	"./fwi_offload.c"
	"./fwi_propagator.c"
 	"./fwi_main.c")

target_link_libraries(fwi.intel64 "-lm")
install (TARGETS fwi.intel64 DESTINATION "${CMAKE_SOURCE_DIR}/")
