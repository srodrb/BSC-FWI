#pragma once

#include <stdlib.h>

typedef int integer;
typedef float real;

/* simulation parameters */
extern const integer WRITTEN_FIELDS;
extern const integer HALO;
extern const integer SIMD_LENGTH;
extern const real    IT_FACTOR;
extern const real    IO_CHUNK_SIZE;

extern const size_t ALIGN_INT;
extern const size_t ALIGN_INTEGER;
extern const size_t ALIGN_REAL;


