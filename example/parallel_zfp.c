/**
 *  @file test_compress_parallel.c
 *  @author Dingwen Tao
 *  @date January, 2017
 *  @brief This is an example of using compression interface in parallel
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zfp.h"
#include "mpi.h"
#include "rw.h"

unsigned char * zfp_compress_3D(float * array, double tolerance, size_t r1, size_t r2, size_t r3, size_t *out_size){
	int status = 0;    /* return value: 0 = success */
	zfp_type type;     /* array scalar type */
	zfp_field* field;  /* array meta data */
	zfp_stream* zfp;   /* compressed stream */
	void* buffer;      /* storage for compressed stream */
	size_t bufsize;    /* byte size of compressed buffer */
	bitstream* stream; /* bit stream to write to or read from */
	size_t zfpsize;    /* byte size of compressed stream */

	/* allocate meta data for the 3D array a[nz][ny][nx] */
	type = zfp_type_float;
	field = zfp_field_3d(array, type, r3, r2, r1);

	/* allocate meta data for a compressed stream */
	zfp = zfp_stream_open(NULL);

	/* set compression mode and parameters via one of three functions */
	/*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
	/*  zfp_stream_set_precision(zfp, precision); */
	zfp_stream_set_accuracy(zfp, tolerance);

	/* allocate buffer for compressed data */
	bufsize = zfp_stream_maximum_size(zfp, field);
	buffer = malloc(bufsize);

	/* associate bit stream with allocated buffer */
	stream = stream_open(buffer, bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);

	zfpsize = zfp_compress(zfp, field);
    if (!zfpsize) {
      fprintf(stderr, "compression failed\n");
      status = 1;
    }	

	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
	*out_size = zfpsize;
	return (unsigned char *)buffer;
}

float * zfp_decompress_3D(unsigned char * comp_data, double tolerance, size_t buffer_size, size_t r1, size_t r2, size_t r3){
	int status = 0;    /* return value: 0 = success */
	zfp_type type;     /* array scalar type */
	zfp_field* field;  /* array meta data */
	zfp_stream* zfp;   /* compressed stream */
	void* buffer;      /* storage for compressed stream */
	size_t bufsize;    /* byte size of compressed buffer */
	bitstream* stream; /* bit stream to write to or read from */
	size_t zfpsize;    /* byte size of compressed stream */

	/* allocate meta data for the 3D array a[nz][ny][nx] */
	float * array = (float *) malloc(r1 * r2 * r3 * sizeof(float));
	type = zfp_type_float;
	field = zfp_field_3d(array, type, r3, r2, r1);

	/* allocate meta data for a compressed stream */
	zfp = zfp_stream_open(NULL);

	/* set compression mode and parameters via one of three functions */
	/*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
	/*  zfp_stream_set_precision(zfp, precision); */
	zfp_stream_set_accuracy(zfp, tolerance);

	/* allocate buffer for compressed data */
	bufsize = zfp_stream_maximum_size(zfp, field);
	// buffer = malloc(bufsize);
	buffer = (void *) comp_data;
	bufsize = buffer_size;

	/* associate bit stream with allocated buffer */
	stream = stream_open(buffer, bufsize);
	zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);

    if (!zfp_decompress(zfp, field)) {
      fprintf(stderr, "decompression failed\n");
      status = 1;
    }
	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
	return array;
}

// USAGE
// mpirun -np 16 parallel sz.config folder_num r3 r2 r1
int main(int argc, char * argv[])
{

	size_t r5=0,r4=0,r3=0,r2=0,r1=0;
	char *cfgFile;

	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);	
	
	if(argc < 3)
	{
		printf("Test case: testfloat_compress [config_file] [srcFilePath] [dimension sizes...]\n");
		printf("Example: testfloat_compress sz.config testfloat_8_8_128.dat 8 8 128\n");
		exit(0);
	}

	cfgFile=argv[1];
	
	if(argc>=4)
	  r1 = atoi(argv[3]); //8
	if(argc>=5)
	  r2 = atoi(argv[4]); //8
	if(argc>=6)
	  r3 = atoi(argv[5]); //128
	if(argc>=7)
	  r4 = atoi(argv[6]);
	if(argc>=8)
	  r5 = atoi(argv[7]);
	
	if (world_rank == 0) printf ("Start parallel compressing ... \n");
	if (world_rank == 0) printf("size: %d\n", world_size);
	double start, end;
	double costReadOri = 0.0, costReadZip = 0.0, costWriteZip = 0.0, costWriteOut = 0.0, costComp = 0.0, costDecomp = 0.0;

	MPI_Barrier(MPI_COMM_WORLD);

	int total_folder_num = atoi(argv[2]);
	int rank_folder_num = total_folder_num / world_size;
	int count = 0;
	char file[6][30] ={"dark_matter_density.log10.dat", "temperature.dat", "baryon_density.log10.dat", "velocity_x.dat", "velocity_y.dat", "velocity_z.dat"};
	double tolerance[6] = {20, 1150000, 7, 8500000, 8500000, 8500000};
	char folder[50] = "/lcrc/project/ECP-EZ/public/compression/datasets";
	char filename[100];
	char zip_filename[100];
	// char out_filename[100];
	size_t inSize, outSize; 
	size_t nbEle;
	int status;
	size_t offset = total_folder_num;
	if(offset > 2048) offset = 0; 
	int var_num = 2;
	while (count < rank_folder_num) 
	{
		int folder_index = world_rank * rank_folder_num + count + offset;
		for(int i=0; i<1; i++){
			sprintf(filename, "%s/%d/%s", folder, folder_index, file[var_num]);
			sprintf(zip_filename, "%s/%d/%s.zfp", folder, folder_index, file[var_num]);
			// sprintf(out_filename, "%s/%d/%s.sz.out", folder, i, file[i]);
			// printf("%s\n", filename);
			// printf("%s\n", zip_filename);
			// printf("%s\n", out_filename);

			// Read Input Data
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0) start = MPI_Wtime();
			float *dataIn = readFloatData(filename, &nbEle, &status);
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0){
				end = MPI_Wtime();
				costReadOri += end - start;
			}
			
			// Compress Input Data
			if (world_rank == 0) printf ("Compressing %s\n", filename);
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0) start = MPI_Wtime();
			// unsigned char *bytesOut = SZ_compress_args(SZ_FLOAT, dataIn, &outSize, REL, 0, rel_bound[0], 0, 0, r5, r4, r3, r2, r1);
			unsigned char * bytesOut = zfp_compress_3D(dataIn, tolerance[var_num], r1, r2, r3, &outSize);
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0){
				end = MPI_Wtime();
				costComp += end - start;
			}
			free (dataIn);

			// Write Compressed Data
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0) start = MPI_Wtime();
			writeByteData(bytesOut, outSize, zip_filename, &status);
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0){
				end = MPI_Wtime();
				costWriteZip += end - start;
			}
			free(bytesOut);
			// Read Compressed Data
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0) start = MPI_Wtime();
			unsigned char *bytesIn = readByteData(zip_filename, &inSize, &status);
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0){
				end = MPI_Wtime();
				costReadZip += end - start;
			}
			// Decompress Compressed Data
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0) start = MPI_Wtime();
			// float *dataOut = SZ_decompress(SZ_FLOAT, bytesIn, inSize, r5, r4, r3, r2, r1);
			float * dataOut = zfp_decompress_3D(bytesIn, tolerance[0], inSize, r1, r2, r3);
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0){
				end = MPI_Wtime();
				costDecomp += end - start; 
			}
			free(bytesIn);
			// Write Decompressed Data
			// start = MPI_Wtime();
			// writeFloatData_inBytes(dataOut, nbEle, out_filename, &status);
			// end = MPI_Wtime();
			// costWriteOut += end - start;
			free(dataOut);
			MPI_Barrier(MPI_COMM_WORLD);
		}
		count ++;
	}
	

	double globalcostReadOri, globalcostReadZip, globalcostWriteZip, globalcostWriteOut, globalcostComp, globalcostDecomp;

	MPI_Reduce(&costReadOri, &globalcostReadOri, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&costReadZip, &globalcostReadZip, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&costWriteZip, &globalcostWriteZip, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&costWriteOut, &globalcostWriteOut, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&costComp, &globalcostComp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&costDecomp, &globalcostDecomp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (world_rank == 0)
	{
		printf ("ZFP Finish parallel compressing.\n");
		printf ("Timecost of reading original files = %.2f seconds\n", costReadOri);
		printf ("Timecost of reading compressed files = %.2f seconds\n", costReadZip);
		printf ("Timecost of writing compressed files = %.2f seconds\n", costWriteZip);
		printf ("Timecost of writing decompressed files = %.2f seconds\n", costWriteOut);
		printf ("Timecost of compressing using %d processes = %.2f seconds\n", world_size, costComp);
		printf ("Timecost of decompressing using %d processes = %.2f seconds\n\n", world_size, costDecomp);
	}


	MPI_Finalize();

	return 0;
}
