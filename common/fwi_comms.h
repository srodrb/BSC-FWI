#pragma once

#if defined(DISTRIBUTED_MEMORY_IMPLEMENTATION)
	#include "fwi_propagator.h"
	#include "mpi.h"
	
	/*
	NAME:exchange_boundaries
	PURPOSE: data exchanges between the boundary layers of the analyzed volume
	
	v                   (in) struct containing velocity arrays (4 points / cell x 3 components / point = 12 arrays)
	numElement          (in) Number of elements to exchange
	rank                (in) MPI process identifier
	numTasks            (in) MPI number of tasks
	idxt                (in) identifier related to the folder
	nyf                 (in) final plane to be exchanged
	ny0                 (in) intial plane to be exchanged
	
	RETURN none
	*/

	#define EXCHANGE(sendbuf, recvbuf, dst, src, count) {  exchange_buffer((sendbuf),(recvbuf),(dst),(src),(count), __FILE__, __LINE__); }
	
	integer exchange_buffer (const real*   sendbuf, 
	                               real*   recvbuf, 
	                         const integer dst, 
	                         const integer src, 
	                         const integer message_size,
	                         const char*   file,
	                         const integer line);

	/* --------------- BOUNDARY EXCHANGES ---------------------------------------- */

	/*
	NAME:exchange_boundaries
	PURPOSE: data exchanges between the boundary layers of the analyzed volume
	
	v                   (in) struct containing velocity arrays (4 points / cell x 3 components / point = 12 arrays)
	plane_size          (in) Number of elements per plane to exchange
	rank                (in) rank id (CPU id)
	nranks              (in) number of CPUs
	nyf                 (in) final plane to be exchanged
	ny0                 (in) intial plane to be exchanged
	
	RETURN none
	*/
	void exchange_velocity_boundaries ( v_t v,
                                    const int plane_size,
                                    int rank,
                                    int nranks,
                                    int nyf,
                                    int ny0);
	/*
	NAME:exchange_stress_boundaries
	PURPOSE: data exchanges between the boundary layers of the analyzed volume
	
	s                   (in) struct containing stress arrays (4 points / cell x 6 components / point = 24 arrays)
	plane_size          (in) Number of elements per plane to exchange
	rank                (in) rank id (CPU id)
	nranks              (in) number of CPUs
	nyf                 (in) final plane to be exchanged
	ny0                 (in) intial plane to be exchanged
	
	RETURN none
	*/
	void exchange_stress_boundaries   ( s_t s,
                                    const int plane_size,
                                    const int rank,
                                    const int nranks,
                                    const int nyf,
                                    const int ny0);

#endif /* end of DISTRIBUTED MEMORY definition */

