/**
 *  @file test_compress.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "sz.h"
#include "rw.h"

int main(int argc, char * argv[])
{
    size_t r5=0,r4=0,r3=0,r2=0,r1=0;
    char outDir[640], oriFilePath[640], outputFilePath[640];
    char *cfgFile;
    
    if(argc < 3)
    {
		printf("Test case: changeDataSize [config_file] [srcFilePath] [dimension sizes...]\n");
		printf("Example: changeDataSize sz.config testfloat_8_8_128.dat 8 8 128\n");
		exit(0);
    }
   
    cfgFile=argv[1];
    sprintf(oriFilePath, "%s", argv[2]);
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
   
    printf("cfgFile=%s\n", cfgFile); 
    int status = SZ_Init(cfgFile);
    if(status == SZ_NSCS)
		exit(0);
    sprintf(outputFilePath, "%s.sz", oriFilePath);
   
    size_t nbEle;
    float *data = readFloatData(oriFilePath, &nbEle, &status);
    int i = 0;
    for(i=0;i<10;i++)
	    printf("data[$d]=%f\n", i, data[i]);
    if(status != SZ_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
   
    size_t outSize; 
    writeFloatData_inBytes(data, 1000000, "test.dat", &status);
    printf("done\n");
    free(data);
    SZ_Finalize();
    
    return 0;
}
