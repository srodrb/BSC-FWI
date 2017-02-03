#include "fwi_common.h"

int max_int( int a, int b)
{
    return ((a >= b) ? a : b);
};

inline double dtime(void)
{
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday( &mytime, (struct timezone*) 0);
    tseconds = (double) (mytime.tv_sec + (double) mytime.tv_usec * 1.0e-6);
    return (tseconds);
};

inline double TOGB(size_t bytes)
{
    return (bytes / (1024.f * 1024.f * 1024.f));
};

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
                          char *outputfolder)
{
    FILE *fp = safe_fopen(fname, "r", __FILE__, __LINE__ );

    CHECK( fscanf( fp, "%f\n", (real*) lenz   ) );
    CHECK( fscanf( fp, "%f\n", (real*) lenx   ) );
    CHECK( fscanf( fp, "%f\n", (real*) leny   ) );
    CHECK( fscanf( fp, "%f\n", (real*) vmin   ) );
    CHECK( fscanf( fp, "%f\n", (real*) srclen ) );
    CHECK( fscanf( fp, "%f\n", (real*) rcvlen ) );
    CHECK( fscanf( fp, "%d\n", (int*)  nshots  ) );
    CHECK( fscanf( fp, "%d\n", (int*)  ngrads  ) );
    CHECK( fscanf( fp, "%d\n", (int*)  ntests  ) );
    CHECK( fscanf( fp, "%f\n", (real*)  workmem ) );
    CHECK( fscanf( fp, "%f\n", (real*)  slavemem) );
    CHECK( fscanf( fp, "%s\n",  outputfolder  ) );

		print_debug("Worker memory %f slave memory %f", *workmem, *slavemem);

    print_debug("Len (z,x,y) (%.2f,%.2f,%.2f)\n \
				vmin %.2f scrlen %.2f rcvlen %.2f outputfolder '%s'\n \
				worker memory %.5f GB slave memory %.5fGB",
      *lenz, *lenx, *leny, *vmin, *srclen, *rcvlen, outputfolder, *workmem, *slavemem );

    fclose(fp);
};

/*
  This function is intended to round up a number (number) to the nearest multiple of the register
  size. In this way, we assure that the dimensions of the domain are suited to the most aggressive
  compiler optimizations.
 */
integer roundup(integer number, integer multiple)
{
    if (multiple == 0)
        return number;

    int remainder = number % multiple;
    if (remainder == 0)
        return number;

    return number + multiple - remainder;
};

/*
 NAME:allocate_shot_memory
 PURPOSE: Create files to store final preconditioner and gradient results. Must be initialized with zeroes.

 outputfolder     (in) folder where snapshot data is store
 VolumeMemory     (in) memory needed to store the domain

 RETURN none
 */
void create_output_volumes(char *outputfolder, integer VolumeMemory)
{
    print_debug("Creating output files in %s", outputfolder);

#ifdef DO_NOT_PERFORM_IO
    print_info("Warning: we are not doing any IO here.");
#else
    char fnamePrecond[300], fnameGradient[300];

    sprintf( fnameGradient, "%s/resultGradient.res", outputfolder);
    sprintf( fnamePrecond , "%s/resultPrecond.res", outputfolder);

    FILE *fGradient = safe_fopen( fnameGradient, "wb", __FILE__, __LINE__ );
    FILE *fPrecond  = safe_fopen( fnamePrecond , "wb", __FILE__, __LINE__ );

    int numIts = ceil( VolumeMemory / IO_CHUNK_SIZE );

    /* create buffer array */
    real *tmparray = (real*) __malloc( ALIGN_REAL, IO_CHUNK_SIZE );

    /* perform the accumulation of the chunks */
    for (int i=0; i<numIts; i++) {
        safe_fwrite(tmparray, 1, IO_CHUNK_SIZE, fGradient, __FILE__, __LINE__ );
        safe_fwrite(tmparray, 1, IO_CHUNK_SIZE, fPrecond , __FILE__, __LINE__ );
    }

    __free(tmparray);

    // close files
    safe_fclose( fnameGradient, fGradient, __FILE__, __LINE__ );
    safe_fclose( fnamePrecond , fPrecond , __FILE__, __LINE__ );
#endif
}

/*
 NAME:create_folder
 PURPOSE:During execution creates temporal folders to organize necessary data for the execution

 folder      (in) name of the temporal folder created
 parent_rank (in) name of the rank related to the data archived in to the folder
 shotID      (in) identifier of the shot related to the data to be archived in to the folder

 RETURN none
 */
void create_folder(const char *folder)
{
    if (mkdir_p(folder) != 0) {
        print_error("cant create folder %s (Err code: %s)", folder, strerror(errno));
        exit(-1);
    }
    print_debug("Folder '%s' created",folder);
};

/*
 NAME: mkdir_p
 PURPOSE: creates the hierarchy of folders requested, if they do not exist.

 RETURN 0 if successful, !=0 otherwise
 */
int mkdir_p(const char *dir)
{
    char tmp[256];
    char *p = NULL;
    size_t len;

    snprintf(tmp, sizeof(tmp),"%s",dir);
    len = strlen(tmp);

    if(tmp[len - 1] == '/')
        tmp[len - 1] = 0;

    for(p = tmp + 1; *p; p++) {
        if(*p == '/') {
            *p = 0;
            int rc = mkdir(tmp, S_IRWXU);
            if (rc != 0 && errno != EEXIST) {
                print_error("Error creating folder %s (Err code %s)", tmp, strerror(errno));
                return -1;
            }

            *p = '/';
        }
    }

    int rc = mkdir(tmp, S_IRWXU);
    if (rc != 0 && errno != EEXIST) {
        print_error("Error creating folder %s (Err code %s)", tmp, strerror(errno));
        return -1;
    }

    return 0;
}

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
                           char    *outputfolder, 
                           real    waveletFreq)
{
    char name[200];

    sprintf(name, "%s/shotparams_%2.1f.%05d.dat", 
            outputfolder, waveletFreq, shotid);

    print_debug("Storing parameters for freq %.3fHz, shot %d into %s", waveletFreq, shotid, name);

    FILE *fp = safe_fopen(name, "w", __FILE__, __LINE__);

    fprintf(fp, "%f\n",  (real   ) *dz     );
    fprintf(fp, "%f\n",  (real   ) *dx     );
    fprintf(fp, "%f\n",  (real   ) *dy     );
    fprintf(fp,  I"\n", (integer) *dimmz  );
    fprintf(fp,  I"\n", (integer) *dimmx  );
    fprintf(fp,  I"\n", (integer) *dimmy  );
    fprintf(fp, "%d\n",  (int    ) *nt_fwd );
    fprintf(fp, "%d\n",  (int    ) *nt_bwd );
    fprintf(fp, "%f\n",  (real   ) *dt     );
    fprintf(fp, "%d\n",  (int    ) *stacki );

    fclose(fp);
};

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
                          char    *outputfolder,
                          real    waveletFreq)
{
    char name[200];

    sprintf(name, "%s/shotparams_%2.1f.%05d.dat", outputfolder, waveletFreq, shotid);
    print_debug("Loading parameters for freq %.3fHz, shot %d into %s", waveletFreq, shotid, name);

    FILE *fp = safe_fopen(name, "r", __FILE__, __LINE__);

    CHECK( fscanf(fp, "%f\n",  (real*   ) dz     ) );
    CHECK( fscanf(fp, "%f\n",  (real*   ) dx     ) );
    CHECK( fscanf(fp, "%f\n",  (real*   ) dy     ) );
    CHECK( fscanf(fp,  I"\n",  (integer*) dimmz  ) );
    CHECK( fscanf(fp,  I"\n",  (integer*) dimmx  ) );
    CHECK( fscanf(fp,  I"\n",  (integer*) dimmy  ) );
    CHECK( fscanf(fp, "%d\n",  (int*    ) nt_fwd ) );
    CHECK( fscanf(fp, "%d\n",  (int*    ) nt_bwd ) );
    CHECK( fscanf(fp, "%f\n",  (real*   ) dt     ) );
    CHECK( fscanf(fp, "%d\n",  (int*    ) stacki ) );

    safe_fclose( name, fp, __FILE__, __LINE__);
};

void load_freqlist( const char* filename, int *nfreqs, real **freqlist )
{
    int count  = 0;
    real freq;

    FILE *freqfile = safe_fopen( filename, "r", __FILE__, __LINE__);

    while( 1 )
    {
        int n = fscanf( freqfile, "%f", &freq);

        if ( n == 1 )
        {
            count += 1;
        }
        else if (errno != 0)
        {
            print_error("Error while reading freqlist file");
            break;
        }
        else if ( n == EOF )
        {
            break;
        }
    }


    /* Allocate memory for frequencies */
    *freqlist = (real*) __malloc( ALIGN_REAL, count * sizeof(real));

    /* return to initial position */
    fseek( freqfile, 0, SEEK_SET);
    count = 0;



    /* read again the file, storing the wavelet frequencies this time */
    while( 1 )
    {
        int n = fscanf( freqfile, "%f", &freq);

        if ( n == 1 )
        {
            (*freqlist)[count++] = freq;
        }
        else if (errno != 0)
        {
            print_error("Error while reading freqlist file");
            break;
        }
        else if ( n == EOF )
        {
            break;
        }
    }
    fclose( freqfile );

    *nfreqs = count;

    print_info("A total of %d frequencies were found...", *nfreqs );
    for( int i=0; i<count; i++)
        print_info("     %.2f Hz", (*freqlist)[i] );
};

void* __malloc( size_t alignment, const integer size)
{
    void *buffer;
    int error;
    
    if( (error=posix_memalign( &buffer, alignment, size)) != 0)
    {
        print_error("Cant allocate buffer correctly");
        abort();
    }
    
    return (buffer);
};

void __free ( void* ptr)
{
    free( ptr );
};

/*
	Reads an environmental variable.
 */
char *read_env_variable (const char* varname)
{	
	char* s = getenv(varname);
	
	if ( s == NULL )
	{
		fprintf(stderr, "%s: ERROR: unable to read  %s env. var\n", __FUNCTION__, varname);
		abort();
	}

#ifdef DEBUG
	printf("%s: %s variable value is :%s\n", __FUNCTION__, varname, s);
#endif

	return (s);
};


FILE* safe_fopen(const char *filename, char *mode, char* srcfilename, int linenumber)
{
    FILE* temp = fopen( filename, mode);
    
    if( temp == NULL){
        print_error("Cant open filename %s, openmode '%s' (called from %s - %d)", filename, mode, srcfilename, linenumber);
        exit(-1);
    }
    return temp;
};

void safe_fclose ( const char *filename, FILE* stream, char* srcfilename, int linenumber)
{
    if ( fclose( stream ) != 0)
    {
        print_error("Cant close filename %s (called from %s - %d)", filename, srcfilename, linenumber);
        abort();
    }

/*if ( unlink(filename)  != 0)
  {
    fprintf(stderr, "%s:%d: Cant unlink file %s correctly!\n", srcfilename, linenumber, filename );
    abort();
  }*/
};


inline void safe_fwrite (void *ptr, size_t size, size_t nmemb, FILE *stream, char* srcfilename, int linenumber)
{
#if defined(DO_NOT_PERFORM_IO)
    print_info("Warning: we are not doing any IO (called from %s).", __FUNCTION__);
#else
    if( stream == NULL ){
        print_error("Invalid stream\n");
        abort();
    }
    double start_t = dtime();
    size_t res = fwrite( ptr, size, nmemb, stream);
    double end_t = dtime() - start_t;
    
   
    if( res != nmemb )
    {
        print_error("Error while fwrite (called from %s - %d)", srcfilename, linenumber );
        abort();
    }
    
    double mbytes = (1.0 * size * nmemb) / (1024.0 * 1024.0);
		print_stats("WRITE Time %lf, elements %lu bytes %lu, MB %lf MB/s %lf", end_t, nmemb, size*nmemb, mbytes, mbytes / end_t);
#endif
};

inline void safe_fread (void *ptr, size_t size, size_t nmemb, FILE *stream, char* srcfilename, int linenumber)
{
#if defined(DO_NOT_PERFORM_IO)
    print_info("Warning: we are not doing any IO (called from %s).", __FUNCTION__);
#else
    if( stream == NULL ){
        print_error("Invalid\n");
        abort();
    }
    
		double start_t = dtime();
		size_t res = fread( ptr, size, nmemb, stream);
		double end_t = dtime() - start_t;

    
		if( res != nmemb )
    {
        print_error("Error while fread (called from %s - %d)", srcfilename, linenumber);
        print_error("Trying to read %lu elements, only %lu were recovered", nmemb, res);
        abort();
    }
    
    double mbytes = (1.0 * size * nmemb) / (1024.0 * 1024.0);
		print_stats("READ Time %lf, elements %lu bytes %lu, MB %lf MB/s %lf", end_t, nmemb, size*nmemb, mbytes, mbytes / end_t );
#endif
};



void fwi_writelog(const char *SourceFileName, 
                  const int LineNumber,
                  const char *FunctionName,
                  const char* MessageHeader,
                  const char *fmt,
                  ...)
{
    const char LogFileName[] = "./fwi.log";
    
    FILE *fp = safe_fopen ( LogFileName, "a", __FILE__, __LINE__ );
    
    va_list args;
    va_start(args, fmt);
    fprintf(fp, "%s :[%s:%d:%s] :: ", MessageHeader, SourceFileName, LineNumber, FunctionName );
    vfprintf(fp, fmt, args);
    fprintf(fp, "\n");
    va_end(args);
    
    safe_fclose ( LogFileName, fp, __FILE__, __LINE__);
};
