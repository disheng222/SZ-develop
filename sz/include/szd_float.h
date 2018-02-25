/**
 *  @file szd_float.h
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief Header file for the szd_float.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _SZD_Float_H
#define _SZD_Float_H

#ifdef __cplusplus
extern "C" {
#endif

#include "TightDataPointStorageF.h"

void decompressDataSeries_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps);
void decompressDataSeries_float_2D(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps);
void decompressDataSeries_float_3D(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps);
void decompressDataSeries_float_4D(float** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageF* tdps);
void decompressDataSeries_float_3D_RA(float** data, size_t r1, size_t r2, size_t r3, unsigned char * comp_data);
void decompressDataSeries_float_1D_RA(float** data, size_t r1, unsigned char * comp_data);
unsigned short decompressDataSeries_float_1D_RA_block_1D_pred(float * data, float mean, size_t dim_0, int block_dim_0, double realPrecision, int * type, float * unpredictable_data);
//unsigned short decompressDataSeries_float_3D_RA_block(float * data, float mean, size_t dim_0, size_t dim_1, size_t dim_2, int block_dim_0, int block_dim_1, int block_dim_2, double realPrecision, int * type, float * unpredictable_data);
unsigned short decompressDataSeries_float_3D_RA_block_3D_pred(float * data, float mean, size_t dim_0, size_t dim_1, size_t dim_2, int block_dim_0, int block_dim_1, int block_dim_2, double realPrecision, int * type, float * unpredictable_data);
//unsigned short decompressDataSeries_float_3D_RA_block_1D_pred(float * data, float mean, size_t dim_0, size_t dim_1, size_t dim_2, int block_dim_0, int block_dim_1, int block_dim_2, double realPrecision, int * type, float * unpredictable_data);
void getSnapshotData_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps, int errBoundMode);
void getSnapshotData_float_2D(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps, int errBoundMode);
void getSnapshotData_float_3D(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps, int errBoundMode);
void getSnapshotData_float_4D(float** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageF* tdps, int errBoundMode);

int SZ_decompress_args_float(float** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZD_Float_H  ----- */