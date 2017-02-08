/*
 * =====================================================================================
 *
 *       Filename:  fwi_main.c
 *
 *    Description:  Main file of the FWI mockup
 *
 *        Version:  1.0
 *        Created:  10/12/15 10:33:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */

#include "fwi_kernel.h"
#include "fwi_offload.h"
#include "fwi_sched.h"

/*
 * In order to generate a source for injection,
 * /system/support/bscgeo/src/wavelet.c
 * functions can be used.
 */




/*
 * propagator: kernel mode
 * waveletFreq: wavelet frequency in Hz.
 * shotid: shot identificator (integer)
 * outputfolder:
 * nworkers: number of workers required to compute the shot.
 * ppw: number of y-planes per worker.
 */
void kernel ( propagator_t propagator, 
							real waveletFreq, 
							int shot, 
							char* outputfolder, 
							int nworkers )
{
	/* allocate slave nodes */
		booster_alloc_t workers = allocate_workers( nworkers, 1, shot );

 		/* load simulation parameters */		
    real dt,dz,dx,dy;
    integer dimmz, dimmx, dimmy;
    int stacki, forw_steps, back_steps, MaxYPlanesPerWorker;

		const int FIRSTWORKER = 0;
		const int LASTWORKER  = nworkers -1;

    char shotfolder[200];
    sprintf(shotfolder, "%s/shot.%2.2fHz.%03d", outputfolder, waveletFreq, shot);
    load_shot_parameters( shot, &stacki, &dt, 
													&forw_steps, &back_steps, 
													&dz, &dx, &dy, 
													&dimmz, &dimmx, &dimmy, 
													&MaxYPlanesPerWorker,
													outputfolder,
													waveletFreq);

		for(int worker = 0; worker < nworkers; worker++)
		{
			/* Compute the integration limits in order to load the correct slice from the input
			 * velocity model. These are not the limits for the wave propagator! (they are local,
			 * i.e. starts at zero!) */
			const integer y0 = (worker == FIRSTWORKER) ? 0     : (MaxYPlanesPerWorker * worker) - HALO;
			const integer yf = (worker == LASTWORKER) ? dimmy : y0 + MaxYPlanesPerWorker;

			/*
			 * Compute integration limits for the wave propagator. It assumes that the volume
			 * is local, so the indices start at zero
			 */
	    const integer edimmy = (yf - y0);
			const integer nz0 = 0;
			const integer nx0 = 0;
			const integer ny0 = 0;
			const integer nzf = dimmz;
			const integer nxf = dimmx;
			const integer nyf = edimmy;
	    const integer numberOfCells = dimmz * dimmx * edimmy;
		
			print_debug("number of cells in kernel() %d\n", numberOfCells);
	    print_debug("The length of local arrays is " I " cells", numberOfCells);
	
			#pragma omp task onto(workers.intercomm, worker) in(propagator, shot, [200]shotfolder) copy_deps
			{
				real    *rho;
    		v_t     v;
    		s_t     s;
    		coeff_t coeffs;

				/* allocate shot memory */
    		alloc_memory_shot  ( dimmz, dimmx, (nyf - ny0), &coeffs, &s, &v, &rho);

				/* load initial model from a binary file */
				load_local_velocity_model ( waveletFreq, dimmz, dimmx, y0, yf, &coeffs, &s, &v, rho);

				/* Allocate memory for IO buffer */
				// real* io_buffer = (real*) __malloc( ALIGN_REAL, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

    		/* inspects every array positions for leaks. Enabled when DEBUG flag is defined */
    		check_memory_shot  ( dimmz, dimmx, (nyf - ny0), &coeffs, &s, &v, rho);

				/* some variables for timming */
				double start_t, end_t;

   		 	switch( propagator )
   		 	{
      		case( RTM_KERNEL ):
      		{
						start_t = dtime(); 
						propagate_shot ( FORWARD,
                        v, s, coeffs, rho,
                        forw_steps, back_steps -1,
                        dt,dz,dx,dy,
                        nz0, nzf, nx0, nxf, ny0, nyf,
                        stacki,
                        shotfolder,
												NULL,
                        dimmz, dimmx, (nyf - ny0));

						end_t = dtime();
        
						print_stats("Forward propagation finished in %lf seconds", end_t - start_t );

						start_t = dtime();
        		propagate_shot ( BACKWARD,
                        v, s, coeffs, rho,
                        forw_steps, back_steps -1,
                        dt,dz,dx,dy,
                        nz0, nzf, nx0, nxf, ny0, nyf,
                        stacki,
                        shotfolder,
                        NULL,
                        dimmz, dimmx, (nyf - ny0));

						end_t = dtime();

        		print_stats("Backward propagation finished in %lf seconds", end_t - start_t );
						
						/* store gradient and preconditioner fields */	
						// store_field( shotfolder, shotid, GRADIENT      , &v, numberOfCells );
						// store_field( shotfolder, shotid, PRECONDITIONER, &v, numberOfCells );
        		
						break;
      		}
      		case( FM_KERNEL  ):
      		{
						start_t = dtime();

        		propagate_shot ( FWMODEL,
                        v, s, coeffs, rho,
                        forw_steps, back_steps -1,
                        dt,dz,dx,dy,
                        nz0, nzf, nx0, nxf, ny0, nyf,
                        stacki,
                        shotfolder,
                        NULL,
                        dimmz, dimmx, (nyf - ny0));

						end_t = dtime();
        
						print_stats("Forward Modelling finished in %lf seconds", end_t - start_t );
			 			break;
      		}
      		default:
      		{
        		print_error("Invalid propagation identifier");
        		abort();
      		}
    		}
				/* deallocate shot memory */
   			free_memory_shot  ( &coeffs, &s, &v, &rho);
				// __free( io_buffer );
			} /* end of ompss pragma running on workers*/
		} /* end of work scheduling loop */
		#pragma omp taskwait
    
		deep_booster_free(&workers.intercomm);
};


int main(int argc, char *argv[])
{
	/* inialize nanos runtime */
	nanos_mpi_init( &argc, &argv );
	MPI_Comm spawn_comm = MPI_COMM_WORLD;

	schedule_t S = load_schedule( argv[1] );

  for(int i=0; i<S.nfreqs; i++)
  {
		
		real waveletFreq = S.freq[i];
		integer stacki   = S.stacki[i];
		real dt          = S.dt[i];
		integer forws    = S.forws[i];
		integer backs    = S.backs[i];
		real dz          = S.dz[i];
		real dx          = S.dx[i];
		real dy          = S.dy[i];
		integer dimmz    = S.dimmz[i];
		integer dimmx    = S.dimmx[i];
		integer dimmy    = S.dimmy[i];
		integer MaxYPlanesPerWorker = S.ppd[i];
		integer nworkers = S.nworkers[i];

		print_info("At %.2Hz, we'll allocate %d slaves and %d workers", waveletFreq,  S.nshots, S.nworkers[i] );

		/* allocate slave nodes */
		booster_alloc_t slaves = allocate_slaves( S.nshots );

		for(int grad=0; grad<S.ngrads; grad++) /* inversion iterations */
		{
			print_info("Processing %d-th gradient iteration.", grad);

			for(int shot=0; shot<S.nshots; shot++)
			{
				#pragma omp task onto(slaves.intercomm, shot) in([200]S.outputfolder) label(rtm_kernel) copy_deps
				{

					print_debug("stacki %d forws %d backs %d", stacki, forws, backs );

					char shotfolder[200];
          sprintf(shotfolder, "%s/shot.%2.2fHz.%03d", S.outputfolder, waveletFreq, shot);
					create_folder( shotfolder );
					
					store_shot_parameters ( shot, &stacki, &dt, &forws, &backs, 
																	&dz, &dx, &dy, 
																	&dimmz, &dimmx, &dimmy,
																	&MaxYPlanesPerWorker,
																	S.outputfolder,
																	waveletFreq);

					kernel( RTM_KERNEL, waveletFreq, shot, S.outputfolder, nworkers);
					
          print_info("\tGradient loop processed for %d-th shot", shot);
				}
			}
			#pragma omp taskwait

			/* shot gathering */
			// gather_shots( outputfolder, nshots, numberOfCells );
	
			for(int test=0; test<S.ntests; test++)
			{
				print_info("Processing %d-th test iteration.", S.ntests);
				 
				for(int shot=0; shot<S.nshots; shot++)
				{
					#pragma omp task onto(slaves.intercomm, shot) in([200]S.outputfolder) label(test_iterations) copy_deps
					{
						char shotfolder[200];
            sprintf(shotfolder, "%s/test.%05d.shot.%2.2fHz.%03d", 
                    S.outputfolder, test, waveletFreq, shot);
						create_folder( shotfolder );
					
						store_shot_parameters ( shot, &stacki, &dt, &forws, &backs, 
																	&dz, &dx, &dy, 
																	&dimmz, &dimmx, &dimmy, 
																	&MaxYPlanesPerWorker,
																	S.outputfolder,
																	waveletFreq);

						kernel( FM_KERNEL , waveletFreq, shot, S.outputfolder, nworkers);
				
            print_info("\t\tTest loop processed for the %d-th shot", shot);
					}// end of ompss pragma
				}
				#pragma omp taskwait
			}
		} /* end of test loop */
    deep_booster_free(&slaves.intercomm);
  } /* end of frequency loop */
  print_info("-------- FWI propagator Finished ------------------- \n");

	/* clean UP */
	int rank, spawn_rank, size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_rank( spawn_comm, &spawn_rank );
	
	/* wait for all the processes to finish */
	print_info("Waiting for all processes to finish");
	MPI_Barrier( MPI_COMM_WORLD );

	print_info("Waiting for Nanos runtime to finish");
	/* explicitely finalize nanos runtime */
	nanos_mpi_finalize();

	print_info("End of the program");
	
  return 0;
}
