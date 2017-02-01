#include "fwi_constants.h"

/* extern variables declared in the header file */
const integer  WRITTEN_FIELDS =   12; /* >= 12.  */
const integer  HALO           =    4; /* >= 4    */ 
const integer  SIMD_LENGTH    =    8; /* # of real elements fitting into regs */
const real     IT_FACTOR      = 0.02;
const real     IO_CHUNK_SIZE  = 1024.f * 1024.f;

const size_t ALIGN_INT     = 16;
const size_t ALIGN_INTEGER = 16;
const size_t ALIGN_REAL    = 64;


