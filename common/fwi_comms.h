#pragma once

#include "fwi_propagator.h"

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





#if defined(DISTRIBUTED_MEMORY_IMPLEMENTATION)
	/*
	 * This region implements MPI functions
	 */

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

	#define EXCHANGE(sendbuf, recvbuf, dst, src, count) {                           \
		exchange_buffer((sendbuf),(recvbuf),(dst),(src),(count), __FILE__, __LINE__); \
	}
	
	inline void exchange_buffer (const real*   sendbuf, 
	                          real*   recvbuf, 
	                    const integer dst, 
	                    const integer src, 
	                    const integer message_size,
	                    const char*   file,
	                    const integer line)
	{
			int err;
	    int tag = 100;
	    
	    print_debug( "         [BEFORE]MPI sendrecv [count:%d][dst:%d][src:%d] %s : %d", 
					message_size,  dst, src, file, line);
	
	    MPI_Status  statuses[2];
	    MPI_Request requests[2];
	    
	    MPI_Irecv( recvbuf, message_size, MPI_FLOAT, dst, tag, MPI_COMM_WORLD, &requests[0] );
	    MPI_Isend( sendbuf, message_size, MPI_FLOAT, dst, tag, MPI_COMM_WORLD, &requests[1] );
	    err = MPI_Waitall(2, requests, statuses);
	
	    print_debug( "         [AFTER ]MPI sendrecv                          %s : %d", 
					file, line);    
	
	    if ( err != MPI_SUCCESS )
			{
				print_error("MPI error %d!", err);
				abort();
			}
	};

#else
	#define EXCHANGE(sendbuf, recvbuf, dst, src, count)
#endif /* end of DISTRIBUTED MEMORY definition */
		

