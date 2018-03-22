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


int main(int argc, char * argv[])
{

	int r5=0,r4=0,r3=0,r2=0,r1=0;
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
	
	sprintf(inDir, "%s", argv[2]);
	sprintf(outDir, "%s/out", inDir);
	if (world_rank == 0) printf ("Input directory = %s\n", inDir);
	if (world_rank == 0) printf ("Output directory = %s\n", outDir);

	
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

	if (world_rank == 0) printf("cfgFile=%s\n", cfgFile); 
	
	SZ_Init(cfgFile);

	if (world_rank == 0) printf ("Start parallel compressing ... \n");

	double start, end;
	double costReadOri = 0.0, costReadZip = 0.0, costWriteZip = 0.0, costWriteOut = 0.0, costComp = 0.0, costDecomp = 0.0;

	MPI_Barrier(MPI_COMM_WORLD);

	int total_folder_num = 8192;
	int rank_folder_num = total_folder_num / world_size;
	int count = 0;
	char file[6][30] ={"dark_matter_density.dat", "temperature.dat", "baryon_density.dat", "velocity_x.dat", "velocity_y.dat", "velocity_z.dat"};

	char folder[30] = "/global/cscratch1/sd/shdi/xin";
	char filename[100];
	char zip_filename[100];
	char out_filename[100];
	size_t inSize, outSize; 
	size_t nbEle;
	int status;
	while (count < rank_folder_num) 
	{
		int folder_index = world_rank * rank_folder_num + count;
		for(int i=0; i<6; i++){
			sprintf(filename, "%s/%d/%s", folder, i, file[i]);
			sprintf(out_filename, "%s/%d/%s.sz", folder, i, file[i]);
			sprintf(out_filename, "%s/%d/%s.sz.out", folder, i, file[i]);
			// printf("%s\n", filename);

			// Read Input Data
			start = MPI_Wtime();
			float *dataIn = readFloatData(filename, &nbEle, &status);
			end = MPI_Wtime();
			costReadOri += end - start;
			MPI_Barrier(MPI_COMM_WORLD);
			
			// Compress Input Data
			if (world_rank == 0) printf ("Compressing %s\n", filename);
			start = MPI_Wtime();
			unsigned char *bytesOut = SZ_compress(SZ_FLOAT, dataIn, &outSize, r5, r4, r3, r2, r1);
			end = MPI_Wtime();
			costComp += end - start;
			free (dataIn);
			MPI_Barrier(MPI_COMM_WORLD);

			// Write Compressed Data
			start = MPI_Wtime();
			writeByteData(bytesOut, outSize, zip_filename, &status);
			end = MPI_Wtime();
			costWriteZip += end - start;
			free(bytesOut);
			MPI_Barrier(MPI_COMM_WORLD);

			// Read Compressed Data
			start = MPI_Wtime();
			unsigned char *bytesIn = readByteData(zip_filename, &inSize, &status);
			end = MPI_Wtime();
			costReadZip += end - start;
			MPI_Barrier(MPI_COMM_WORLD);

			// Decompress Compressed Data
			start = MPI_Wtime();
			float *dataOut = SZ_decompress(SZ_FLOAT, bytesIn, inSize, r5, r4, r3, r2, r1);
			end = MPI_Wtime();
			costDecomp += end - start; 
			free(bytesIn);
			MPI_Barrier(MPI_COMM_WORLD);

			// Write Decompressed Data
			start = MPI_Wtime();
			writeFloatData_inBytes(dataOut, nbEle, out_filename, &status);
			end = MPI_Wtime();
			costWriteOut += end - start;
			free(dataOut);
			MPI_Barrier(MPI_COMM_WORLD);
		}

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
		printf ("Finish parallel compressing.\n");
		printf ("Timecost of reading original files = %.2f seconds\n", globalcostReadOri/world_size);
		printf ("Timecost of reading compressed files = %.2f seconds\n", globalcostReadZip/world_size);
		printf ("Timecost of writing compressed files = %.2f seconds\n", globalcostWriteZip/world_size);
		printf ("Timecost of writing decompressed files = %.2f seconds\n", globalcostWriteOut/world_size);
		printf ("Timecost of compressing using %d processes = %.2f seconds\n", world_size, globalcostComp/world_size);
		printf ("Timecost of decompressing using %d processes = %.2f seconds\n\n", world_size, globalcostDecomp/world_size);
	}

	SZ_Finalize();

	MPI_Finalize();

	return 0;
}
