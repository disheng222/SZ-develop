/**
 *  @file convertBinToPFM.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief Converting floating-point Bin data to PFM format data.
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "sz.h"
#include "rw.h"

int main(int argc, char * argv[])
{
	int x = 1;
        char *y = (char*)&x;

        if(*y==1)
                sysEndianType = LITTLE_ENDIAN_SYSTEM;
        else //=0
                sysEndianType = BIG_ENDIAN_SYSTEM;

	int status = 0;
	size_t r5=0,r4=0,r3=0,r2=0,r1=0;
	char outDir[640], oriFilePath[640], outputFilePath[640];
	char *cfgFile;

	if(argc < 3)
	{
		printf("Test case: convertBinToPFM [srcFilePath] [dimension sizes...]\n");
		printf("Example: convertBinToPFM testfloat_8_8_128.dat 8 8 128\n");
		printf("output: [srcFilePath].pfm\n");
		exit(0);
	}

	sprintf(oriFilePath, "%s", argv[1]);
	if(argc>=3)
		r1 = atoi(argv[2]); //8
	if(argc>=4)
		r2 = atoi(argv[3]); //8
	if(argc>=5)
		r3 = atoi(argv[4]); //128
	if(argc>=6)
		r4 = atoi(argv[5]);
	if(argc>=7)
		r5 = atoi(argv[6]);

	sprintf(outputFilePath, "%s.pfm", oriFilePath);

	size_t nbEle;
	float *data = readFloatData(oriFilePath, &nbEle, &status);
	if(status != SZ_SCES)
	{
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
	}

	convertToPFM_float(data, r5, r4, r3, r2, r1, LITTLE_ENDIAN, outputFilePath, &status);

	if(status != SZ_SCES)
	{
		printf("Error: data file %s cannot be written!\n", outputFilePath);
		exit(0);
	}

	printf("done\n");
	free(data);

	return 0;
}
