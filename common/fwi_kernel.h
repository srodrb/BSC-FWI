/*
 * =====================================================================================
 *
 *       Filename:  fwi_kernel.h
 *
 *    Description:  Kernel propagator implementation
 *
 *        Version:  1.0
 *        Created:  14/12/15 12:10:24
 *       Revision:  none
 *       Compiler:  icc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */
#pragma once

#include "fwi_propagator.h"
#include "fwi_comms.h"

/*
 * Ensures that the domain contains a minimum number of planes.
 * Just needed for debugging when running small cases.
 */
void check_domain_dimensions ( const integer dimmz,
                               const integer dimmx,
                               const integer dimmy);

void set_array_to_random_real(real* restrict array,
                              const integer length);

void set_array_to_constant(real* restrict array,
                           const real value,
                           const integer length);

void alloc_memory_shot( const integer dimmz,
												const integer dimmx,
												const integer dimmy,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
                        real    **rho);

void free_memory_shot( coeff_t *c,
                       s_t     *s,
                       v_t     *v,
                       real    **rho);

void check_memory_shot( const integer dimmz,
												const integer dimmx,
												const integer dimmy,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
                        real    *rho);

/* --------------- I/O RELATED FUNCTIONS -------------------------------------- */
void load_local_velocity_model ( const real    waveletFreq,
													const integer dimmz,
													const integer dimmx,
													const integer FirstYPlane,
													const integer LastYPlane,
                          coeff_t *c,
                          s_t     *s,
                          v_t     *v,
                          real    *rho);

void write_snapshot ( char          *folder,
                      const int     suffix,
                      v_t          *v,
                      const integer dimmz,
											const integer dimmx,
											const integer dimmy);

void read_snapshot ( char          *folder,
                     const int     suffix,
                     v_t          *v,
                     const integer dimmz,
										 const integer dimmx,
										 const integer dimmy);



/* --------------- WAVE PROPAGATOR FUNCTIONS --------------------------------- */

void propagate_shot ( time_d        direction,
                     v_t           v,
                     s_t           s,
                     coeff_t       coeffs,
                     real          *rho,
                     int           timesteps,
                     int           ntbwd,
                     real          dt,
                     real          dzi,
                     real          dxi,
                     real          dyi,
                     integer       nz0,
                     integer       nzf,
                     integer       nx0,
                     integer       nxf,
                     integer       ny0,
                     integer       nyf,
                     integer       stacki,
                     char          *folder,
                     real          *UNUSED(dataflush),
                     integer       dimmz,
                     integer       dimmx,
                     integer       dimmy);

/* --------------- BOUNDARY EXCHANGES ---------------------------------------- */
void EXCHANGE (const real*   sendbuf, 
	                          real*   recvbuf, 
	                    const integer dst, 
	                    const integer src, 
	                    const integer message_size,
	                    const char*   file,
	                    const integer line);


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

