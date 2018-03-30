#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "sz.h"

unsigned char * SZ_compress_float_3D_MDQ_openmp(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size);

void decompressDataSeries_float_3D_openmp(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data);

