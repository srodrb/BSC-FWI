#include "fwi_comms.h"


void exchange_velocity_boundaries ( v_t v, 
                                    const integer plane_size, 
                                    const integer rank,
                                    const integer nranks,
                                    const integer nyf, 
                                    const integer ny0 )
{
    const integer num_planes = HALO;
    const integer nelems     = num_planes * plane_size;

    const integer left_recv  = ny0;
    const integer left_send  = ny0+HALO;

    const integer right_recv = nyf-HALO;
    const integer right_send = nyf-2*HALO;
    
    if ( rank != 0 )
    {
        // [RANK-1] <---> [RANK] communication
        EXCHANGE( &v.tl.u[left_send], &v.tl.u[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.tl.v[left_send], &v.tl.v[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.tl.w[left_send], &v.tl.w[left_recv], rank-1, rank, nelems );

        EXCHANGE( &v.tr.u[left_send], &v.tr.u[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.tr.v[left_send], &v.tr.v[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.tr.w[left_send], &v.tr.w[left_recv], rank-1, rank, nelems );

        EXCHANGE( &v.bl.u[left_send], &v.bl.u[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.bl.v[left_send], &v.bl.v[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.bl.w[left_send], &v.bl.w[left_recv], rank-1, rank, nelems );

        EXCHANGE( &v.br.u[left_send], &v.br.u[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.br.v[left_send], &v.br.v[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.br.w[left_send], &v.br.w[left_recv], rank-1, rank, nelems );
    }

    if ( rank != nranks -1 )  //task to exchange stress boundaries
    {
        //                [RANK] <---> [RANK+1] communication
        EXCHANGE( &v.tl.u[right_send], &v.tl.u[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.tl.v[right_send], &v.tl.v[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.tl.w[right_send], &v.tl.w[right_recv], rank+1, rank, nelems );

        EXCHANGE( &v.tr.u[right_send], &v.tr.u[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.tr.v[right_send], &v.tr.v[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.tr.w[right_send], &v.tr.w[right_recv], rank+1, rank, nelems );

        EXCHANGE( &v.bl.u[right_send], &v.bl.u[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.bl.v[right_send], &v.bl.v[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.bl.w[right_send], &v.bl.w[right_recv], rank+1, rank, nelems );

        EXCHANGE( &v.br.u[right_send], &v.br.u[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.br.v[right_send], &v.br.v[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.br.w[right_send], &v.br.w[right_recv], rank+1, rank, nelems );
    }

		print_debug("Velocity boundaries exchanged correctly");
};

void exchange_stress_boundaries ( s_t s, 
                                  const integer plane_size, 
                                  const integer rank,
                                  const integer nranks,
                                  const integer nyf, 
                                  const integer ny0 )
{
    const integer num_planes = HALO;
    const integer nelems     = num_planes * plane_size;

    const integer left_recv  = ny0;
    const integer left_send  = ny0+HALO;

    const integer right_recv = nyf-HALO;
    const integer right_send = nyf-2*HALO;

    if ( rank != 0 )
    {
        // [RANK-1] <---> [RANK] communication
        EXCHANGE( &s.tl.zz[left_send], &s.tl.zz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tl.xz[left_send], &s.tl.xz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tl.yz[left_send], &s.tl.yz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tl.xx[left_send], &s.tl.xx[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tl.xy[left_send], &s.tl.xy[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tl.yy[left_send], &s.tl.yy[left_recv], rank-1, rank, nelems );

        EXCHANGE( &s.tr.zz[left_send], &s.tr.zz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tr.xz[left_send], &s.tr.xz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tr.yz[left_send], &s.tr.yz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tr.xx[left_send], &s.tr.xx[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tr.xy[left_send], &s.tr.xy[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tr.yy[left_send], &s.tr.yy[left_recv], rank-1, rank, nelems );

        EXCHANGE( &s.bl.zz[left_send], &s.bl.zz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.bl.xz[left_send], &s.bl.xz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.bl.yz[left_send], &s.bl.yz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.bl.xx[left_send], &s.bl.xx[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.bl.xy[left_send], &s.bl.xy[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.bl.yy[left_send], &s.bl.yy[left_recv], rank-1, rank, nelems );

        EXCHANGE( &s.br.zz[left_send], &s.br.zz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.br.xz[left_send], &s.br.xz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.br.yz[left_send], &s.br.yz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.br.xx[left_send], &s.br.xx[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.br.xy[left_send], &s.br.xy[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.br.yy[left_send], &s.br.yy[left_recv], rank-1, rank, nelems );
    }
    
    if ( rank != nranks-1 )
    {
        //                [RANK] <---> [RANK+1] communication
        EXCHANGE( &s.tl.zz[right_send], &s.tl.zz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tl.xz[right_send], &s.tl.xz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tl.yz[right_send], &s.tl.yz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tl.xx[right_send], &s.tl.xx[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tl.xy[right_send], &s.tl.xy[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tl.yy[right_send], &s.tl.yy[right_recv], rank+1, rank, nelems );

        EXCHANGE( &s.tr.zz[right_send], &s.tr.zz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tr.xz[right_send], &s.tr.xz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tr.yz[right_send], &s.tr.yz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tr.xx[right_send], &s.tr.xx[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tr.xy[right_send], &s.tr.xy[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tr.yy[right_send], &s.tr.yy[right_recv], rank+1, rank, nelems );

        EXCHANGE( &s.bl.zz[right_send], &s.bl.zz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.bl.xz[right_send], &s.bl.xz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.bl.yz[right_send], &s.bl.yz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.bl.xx[right_send], &s.bl.xx[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.bl.xy[right_send], &s.bl.xy[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.bl.yy[right_send], &s.bl.yy[right_recv], rank+1, rank, nelems );

        EXCHANGE( &s.br.zz[right_send], &s.br.zz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.br.xz[right_send], &s.br.xz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.br.yz[right_send], &s.br.yz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.br.xx[right_send], &s.br.xx[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.br.xy[right_send], &s.br.xy[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.br.yy[right_send], &s.br.yy[right_recv], rank+1, rank, nelems );
    }

		print_debug("Stress boundaries exchanged correctly");
};
