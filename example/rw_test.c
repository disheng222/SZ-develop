#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rw.h"
#include <time.h>
int main(){

	clock_t start, end;
	double elaplsed_time = 0.0;
	int status;
	size_t nbEle;
	char zfp_file[100] = "/lcrc/project/ECP-EZ/public/compression/datasets/0/dark_matter_density.log10.dat.zfp";
	char sz_file[100] = "/lcrc/project/ECP-EZ/public/compression/datasets/0/dark_matter_density.log10.dat.sz_old";
	printf("read zfp compressed data\n");
	start = clock();
	float * data = readByteData(zfp_file, &nbEle, &status);
	end = clock();
	elapsed_time = ((double)(end - start)) /CLOCKS_PER_SEC;
	printf("zfp read time: %.4f\n", elapsed_time);

	printf("write zfp compressed data\n");
	start = clock();
	writeByteData(data, nbEle, "tmp.zfp", &status);
	end = clock();
	elapsed_time = ((double)(end - start)) /CLOCKS_PER_SEC;
	printf("zfp write time: %.4f\n", elapsed_time);

	//
	printf("read sz compressed data\n");
	start = clock();
	float * sz_data = readByteData(sz_file, &nbEle, &status);
	end = clock();
	elapsed_time = ((double)(end - start)) /CLOCKS_PER_SEC;
	printf("sz read time: %.4f\n", elapsed_time);

	printf("write sz compressed data\n");
	start = clock();
	writeByteData(sz_data, nbEle, "tmp.sz", &status);
	end = clock();
	elapsed_time = ((double)(end - start)) /CLOCKS_PER_SEC;
	printf("sz write time: %.4f\n", elapsed_time);

	free(data);
	free(sz_data);

}