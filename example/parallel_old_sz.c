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
#include "sz.h"
#include "rw.h"
#include "mpi.h"

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
	
	SZ_Init(NULL);

	if (world_rank == 0) printf ("Start parallel compressing ... \n");
	if (world_rank == 0) printf("size: %d\n", world_size);
	double start, end;
	double costReadOri = 0.0, costReadZip = 0.0, costWriteZip = 0.0, costWriteOut = 0.0, costComp = 0.0, costDecomp = 0.0;

	MPI_Barrier(MPI_COMM_WORLD);

	int total_folder_num = atoi(argv[2]);
	int rank_folder_num = total_folder_num / world_size;
	int count = 0;
	char file[6][30] ={"dark_matter_density.log10.dat", "temperature.dat", "baryon_density.log10.dat", "velocity_x.dat", "velocity_y.dat", "velocity_z.dat"};
	//double rel_bound[6] = {0.09, 0.103, 0.2, 0.006, 0.0105, 0.005};
	double rel_bound[6] = {0.055, 0.0023, 0.017, 0.0018, 0.0018, 0.0018};
	szRandomAccess = SZ_NO_RANDOM_ACCESS;
	char folder[50] = "/lcrc/project/ECP-EZ/public/compression/datasets";
	char filename[100];
	char zip_filename[100];
	// char out_filename[100];
	size_t inSize, outSize; 
	size_t nbEle;
	int status;
	size_t offset = total_folder_num;
	if(offset > 2048) offset = 0; 

	while (count < rank_folder_num) 
	{
		int folder_index = world_rank * rank_folder_num + count + offset;
		for(int i=0; i<6; i++){
			sprintf(filename, "%s/%d/%s", folder, folder_index, file[i]);
			sprintf(zip_filename, "%s/%d/%s.sz_old", folder, folder_index, file[i]);
			// sprintf(out_filename, "%s/%d/%s.sz.out", folder, i, file[i]);
			// printf("%s\n", filename);
			// printf("%s\n", zip_filename);
			// printf("%s\n", out_filename);

			// Read Input Data
			if(world_rank == 0){
				start = MPI_Wtime();
				dataIn = readFloatData(filename, &nbEle, &status);
				end = MPI_Wtime();
				printf("data read time: %.2f\n", end - start);
				start = MPI_Wtime();
				MPI_Bcast(dataIn, nbEle, MPI_FLOAT, 0, MPI_COMM_WORLD);
				end = MPI_Wtime();
				printf("broadcast time: %.2f\n", end - start);
			}
			else{
				nbEle = 512 * 512 * 512;
				dataIn = (float *) malloc(nbEle * sizeof(float));
				MPI_Bcast(dataIn, nbEle, MPI_FLOAT, 0, MPI_COMM_WORLD);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0){
				end = MPI_Wtime();
				costReadOri += end - start;
			}
			
			// Compress Input Data
			stateNum = 65536;
			allNodes = stateNum * 2;
			if (world_rank == 0) printf ("Compressing %s\n", filename);
			MPI_Barrier(MPI_COMM_WORLD);
			if(world_rank == 0) start = MPI_Wtime();
			// unsigned char *bytesOut = SZ_compress(SZ_FLOAT, dataIn, &outSize, r5, r4, r3, r2, r1);
			unsigned char *bytesOut = SZ_compress_args(SZ_FLOAT, dataIn, &outSize, REL, 0, rel_bound[var_num], 0, 0, r5, r4, r3, r2, r1);
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
			float *dataOut = SZ_decompress(SZ_FLOAT, bytesIn, inSize, r5, r4, r3, r2, r1);
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
		printf ("SZ_OLD Finish parallel compressing.\n");
		printf ("Timecost of reading original files = %.2f seconds\n", costReadOri);
		printf ("Timecost of reading compressed files = %.2f seconds\n", costReadZip);
		printf ("Timecost of writing compressed files = %.2f seconds\n", costWriteZip);
		printf ("Timecost of writing decompressed files = %.2f seconds\n", costWriteOut);
		printf ("Timecost of compressing using %d processes = %.2f seconds\n", world_size, costComp);
		printf ("Timecost of decompressing using %d processes = %.2f seconds\n\n", world_size, costDecomp);
	}


	SZ_Finalize();

	MPI_Finalize();

	return 0;
}
