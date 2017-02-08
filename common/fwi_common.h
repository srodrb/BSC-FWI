/*
 * =====================================================================================
 *
 *       Filename:  fwi_common.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  10/12/15 10:34:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */
#pragma once

// When included before <stdlib.h>, solves implicit declaration of posix_memalign()
// http://stackoverflow.com/questions/32438554/warning-implicit-declaration-of-posix-memalign
#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>
#include "fwi_constants.h"

#define I "%d"     // integer printf symbol

typedef enum {RTM_KERNEL, FM_KERNEL} propagator_t;
typedef enum {FORWARD   , BACKWARD, FWMODEL}  time_d;

#if defined(DISTRIBUTED_MEMORY_IMPLEMENTATION)
	#include "mpi.h"
#endif

double TOGB(size_t bytes);

/*  Compiler compatiblity macros */
#ifdef __GNUC__
  /* http://stackoverflow.com/questions/25667901/assume-clause-in-gcc*/ \
    #define __assume(_cond) do { if (!(_cond)) __builtin_unreachable(); } while (0)
#endif

/*  Compiler macro to suppress unused variable warnings:
 *  This throws an ambiguity error in Mercurium!  */
// #ifdef UNUSED
// #elif defined(__GNUC__)
//   #define UNUSED(x) (x) __attribute__((unused))
//#else
    #define UNUSED(x) x
//#endif

#define CHECK(error) { checkErrors((error), __FILE__, __LINE__); }
static inline void checkErrors(const integer error, const char *filename, int line)
{
    if ( error < 0 ) {                     
        fprintf(stderr, "ERROR: %d in %s:%d\n", error, filename, line);
        exit(-1);
    }
};

char *read_env_variable(const char *varname);
FILE* safe_fopen  ( const char *filename, char *mode, char* srcfilename, int linenumber);
void  safe_fclose ( const char *filename, FILE* stream, char* srcfilename, int linenumber);
void  safe_fwrite ( void *ptr, size_t size, size_t nmemb, FILE *stream, char* srcfilename, int linenumber );
void  safe_fread  ( void *ptr, size_t size, size_t nmemb, FILE *stream, char* srcfilename, int linenumber );
integer roundup(integer number, integer multiple);


int max_int( int a, int b);

double dtime(void);

void read_fwi_parameters (const char *fname,
                          real *lenz,
                          real *lenx,
                          real *leny,
                          real *vmin,
                          real *srclen,
                          real *rcvlen,
													int  *nshots,
													int  *ngrads,
													int  *ntests,
													real  *workmem,
													real  *slavemem,
                          char *outputfolder);

void store_shot_parameters(int     shotid,
                           int     *stacki,
                           real    *dt,
                           int     *nt_fwd,
                           int     *nt_bwd,
                           real    *dz,
                           real    *dx,
                           real    *dy,
                           integer *dimmz,
                           integer *dimmx,
                           integer *dimmy,
													 integer *LocalYPlanes,
                           char    *outputfolder,
                           real    waveletFreq);

void load_shot_parameters(int     shotid,
                          int     *stacki,
                          real    *dt,
                          int     *nt_fwd,
                          int     *nt_bwd,
                          real    *dz,
                          real    *dx,
                          real    *dy,
                          integer *dimmz,
                          integer *dimmx,
                          integer *dimmy,
													integer *LocalYPlanes,
                          char    *outputfolder,
                          real    waveletFreq);

void load_freqlist (  const char*  filename,
                      int*   nfreqs,
                      real** freqlist );

void* __malloc ( const size_t alignment, const integer size);
void  __free   ( void *ptr );

void create_output_volumes(char* outputfolder, integer VolumeMemory);

int mkdir_p(const char *dir);

void create_folder(const char *folder);


#define print_error(M, ...)    fwi_writelog(__FILE__, __LINE__, __func__, "ERROR", M, ##__VA_ARGS__)
#define print_info(M, ...)     fwi_writelog(__FILE__, __LINE__, __func__, "INFO ", M, ##__VA_ARGS__)

#if defined(COLLECT_STATS)
	#define print_stats(M, ...)  fwi_writelog(__FILE__, __LINE__, __func__, "STATS", M, ##__VA_ARGS__)
#else
	#define print_stats(M, ...)
#endif

#if defined(DEBUG)
  #define print_debug(M, ...)  fwi_writelog(__FILE__, __LINE__, __func__, "DEBUG", M, ##__VA_ARGS__)
#else
  #define print_debug(M, ...)
#endif

void fwi_writelog(const char *SourceFileName, 
                  const int LineNumber,
                  const char *FunctionName,
                  const char* MessageHeader,
                  const char *fmt,
                  ...);
