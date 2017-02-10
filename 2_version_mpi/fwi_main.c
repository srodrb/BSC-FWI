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
#include <fwi_sched.h>


/*
 * In order to generate a source for injection,
 * /system/support/bscgeo/src/wavelet.c
 * functions can be used.
 */
void kernel( propagator_t propagator, real waveletFreq, int shotid, char* outputfolder, char* shotfolder)
{
    /* find ourselves into the MPI space */
    int rank, ranksize;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &ranksize);

    /* aux variables, just to make it more readable */
    const int FIRSTRANK = 0;
    const int LASTRANK  = ranksize -1;  
    
    /* local variables */
    int stacki;
    double start_t, end_t;
    real dt,dz,dx,dy;
    integer dimmz;
    integer dimmx;
    integer dimmy;
    integer MaxYPlanesPerWorker;
    integer forw_steps;
    integer back_steps;

    /* load shot parameters */
    load_shot_parameters( shotid, &stacki, &dt, &forw_steps, &back_steps, 
				&dz, &dx, &dy,
				&dimmz, &dimmx, &dimmy,
				&MaxYPlanesPerWorker,
				outputfolder, waveletFreq );
    
    /* Max number of y-planes per worker does not differenciate between HALO or
     * maingrid planes! */

    /* Compute the integration limits in order to load the correct slice from the input
     * velocity model. These are not the limits for the wave propagator! (they are local,
     * i.e. starts at zero!) */
    const integer y0 = (rank == FIRSTRANK) ? 0     : (MaxYPlanesPerWorker * rank) - HALO;
    const integer yf = (rank == LASTRANK ) ? dimmy : y0 + MaxYPlanesPerWorker;

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


    /* Compute number of cell for each MPI rank and local integtion limits */
    real    *rho;
    v_t     v;
    s_t     s;
    coeff_t coeffs;


    /* allocate shot memory */
    alloc_memory_shot  ( dimmz, dimmx, (nyf - ny0), &coeffs, &s, &v, &rho);

    /* load initial model from a binary file */
    load_local_velocity_model ( waveletFreq, dimmz, dimmx, y0, yf, &coeffs, &s, &v, rho);

    /* Allocate memory for IO buffer */
    real* io_buffer = (real*) __malloc( ALIGN_REAL, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

    /* inspects every array positions for leaks. Enabled when DEBUG flag is defined */
    check_memory_shot  ( dimmz, dimmx, (nyf - ny0), &coeffs, &s, &v, rho);

    /* Perform forward, backward or test propagations */
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
   			    io_buffer,
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
   			    io_buffer,
   			    dimmz, dimmx, (nyf - ny0));
     
	    end_t = dtime();
    
	    print_stats("Backward propagation finished in %lf seconds", end_t - start_t );

#if defined(DO_NOT_PERFORM_IO)
    
	    print_info("Warning: we are not creating gradient nor preconditioner "
	        	    "fields, because IO is not enabled for this execution" );
#else
    	    if ( rank == 0 ) 
    	    {
	       char fnameGradient[300];
   	       char fnamePrecond[300];
   	       sprintf( fnameGradient, "%s/gradient_%05d.dat", shotfolder, shotid );
   	       sprintf( fnamePrecond , "%s/precond_%05d.dat" , shotfolder, shotid );
   
	       FILE* fgradient = safe_fopen( fnameGradient, "ab", __FILE__, __LINE__ );
   	       FILE* fprecond  = safe_fopen( fnamePrecond , "ab", __FILE__, __LINE__ );
   
	       print_info("Storing local gradient field in %s", fnameGradient );
   	       safe_fwrite( io_buffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, fgradient, __FILE__, __LINE__ );
   
	       print_info("Storing local preconditioner field in %s", fnamePrecond);
    	       safe_fwrite( io_buffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, fprecond , __FILE__, __LINE__ );
    
	       safe_fclose( fnameGradient, fgradient, __FILE__, __LINE__ );
    	       safe_fclose( fnamePrecond , fprecond , __FILE__, __LINE__ );
    	    }
#endif
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
   			    io_buffer,
   			    dimmz, dimmx, (nyf - ny0) );
    
	    end_t = dtime();
    
	    print_stats("Forward Modelling finished in %lf seconds", end_t - start_t );
     
	    break;
    	}
    	default:
     	{
            print_error("Invalid propagation identifier");
    	    abort();
    	}
    } /* end case */

    // liberamos la memoria alocatada en el shot
    free_memory_shot  ( &coeffs, &s, &v, &rho);
    __free( io_buffer );
};

int main(int argc, char* argv[])
{
    double tstart, tend;
    tstart = dtime();
		
    /* Initialize MPI environment */
    MPI_Init ( &argc, &argv );
    int mpi_rank, subdomains;
    MPI_Comm_size( MPI_COMM_WORLD, &subdomains);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);

    /* load parameters from schedule file */
    schedule_t S = load_schedule(argv[1]);

    for(int i=0; i<S.nfreqs; i++)
    {
        /* Process one frequency at a time */
        real waveletFreq = S.freq[i];
	integer stacki   = S.stacki[i];
	real dt          = S.dt[i];

	integer forw_steps = S.forws[i];
	integer back_steps = S.backs[i];

	real dx = S.dx[i];
        real dy = S.dy[i];
        real dz = S.dz[i];

	integer dimmz = S.dimmz[i];
	integer dimmx = S.dimmx[i];
	integer dimmy = S.dimmy[i];

	integer MaxYPlanesPerWorker   = S.ppd[i];

        print_info("\n------ Computing %d-th frequency (%.2fHz).  -----\n", i, waveletFreq); 

        const integer numberOfCells = dimmz * dimmx * dimmx;
        const size_t VolumeMemory  = numberOfCells * sizeof(real) * 58;

        print_stats("Local domain size for freq %f [%d][%d][%d] is %lu bytes (%lf GB)", 
                    waveletFreq, dimmz, dimmx, dimmy, VolumeMemory, TOGB(VolumeMemory) );

        for(int grad=0; grad<S.ngrads; grad++) /* iteracion de inversion */
        {
            print_info("Processing %d-gradient iteration", grad);

            for(int shot=0; shot<S.nshots; shot++)
            {
                char shotfolder[200];
                sprintf(shotfolder, "%s/shot.%2.2fHz.%03d", S.outputfolder, waveletFreq, shot);
                
                if ( mpi_rank == 0 ) 
                {
                    create_folder( shotfolder );

                    store_shot_parameters( shot, &stacki, &dt, &forw_steps, &back_steps,
                                           &dz, &dx, &dy, 
                                           &dimmz, &dimmx, &dimmy, 
					   &MaxYPlanesPerWorker,
                                           S.outputfolder, waveletFreq );
                }

                MPI_Barrier( MPI_COMM_WORLD );

                kernel( RTM_KERNEL, waveletFreq, shot, S.outputfolder, shotfolder);

                print_info("\tGradient loop processed for %d-th shot", shot);
                //update_shot()
            }

           MPI_Barrier( MPI_COMM_WORLD );

	   // gather_shots();
	   
	   for(int test=0; test<S.ntests; test++)
            {
                print_info("\tProcessing %d-th test iteration", test);
                
                for(int shot=0; shot<S.nshots; shot++)
                {
                    char shotfolder[200];
                    sprintf(shotfolder, "%s/test.%05d.shot.%2.2fHz.%03d", 
                            S.outputfolder, test, waveletFreq, shot);

                    if ( mpi_rank == 0)
                    {
                        create_folder( shotfolder );    

                        store_shot_parameters( shot, &stacki, &dt, &forw_steps, &back_steps, 
                                               &dz, &dx, &dy, 
                                               &dimmz, &dimmx, &dimmy,
					       &MaxYPlanesPerWorker,
                                               S.outputfolder, waveletFreq );
                    }

                    MPI_Barrier( MPI_COMM_WORLD );

                    kernel( FM_KERNEL , waveletFreq, shot, S.outputfolder, shotfolder);
                
                    print_info("\t\tTest loop processed for the %d-th shot", shot);
                }
            } /* end of test loop */
        } /* end of gradient loop */
    } /* end of frequency loop */

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    tend = dtime() - tstart;

    fprintf(stderr, "FWI Program finished in %lf seconds\n", tend);

    return 0;
}
