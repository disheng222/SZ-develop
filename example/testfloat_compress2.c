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

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;

static int
initialize(int offset, double tolerance)
{
  sz_params parameter;
  parameter.max_quant_intervals = 65536;
  parameter.quantization_intervals = 0;
  parameter.dataEndianType = LITTLE_ENDIAN_DATA;
  parameter.sol_ID = SZ;
  parameter.layers = 1;
  parameter.sampleDistance = 100;
  parameter.predThreshold = 0.99f;
  parameter.offset = offset;
  parameter.szMode = SZ_BEST_SPEED;
  parameter.gzipMode = 1; // DEFAULT_COMPRESSION = -1; BEST_SPEED = 1; BEST_COMPRESSION = 9
  parameter.errorBoundMode = ABS;
  parameter.absErrBound = tolerance;

  // unused, copied from config file
  parameter.relBoundRatio = 1e-4;
  parameter.psnr = 80;
  parameter.pw_relBoundRatio = 1e-2;
  parameter.segment_size = 25;
  parameter.pwr_type = SZ_PWR_MIN_TYPE;

  if (SZ_Init_Params(&parameter)) {
    fprintf(stderr, "initialization failed\n");
    return 0;
  }
  else
    return 1;
}


void cost_start()
{
        gettimeofday(&costStart, NULL);
}

void cost_end()
{
        double elapsed;
        struct timeval costEnd;
        gettimeofday(&costEnd, NULL);
        elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
        totalCost += elapsed;
}


int main(int argc, char * argv[])
{
    size_t r5=0,r4=0,r3=0,r2=0,r1=0;
    char outDir[640], oriFilePath[640], outputFilePath[640];
    char *cfgFile;
    
    if(argc < 3)
    {
		printf("Test case: testfloat_compress [config_file] [srcFilePath] [dimension sizes...]\n");
		printf("Example: testfloat_compress sz.config testfloat_8_8_128.dat 8 8 128\n");
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
    //int status = SZ_Init(cfgFile);
    int status = initialize(0, 1E-3);
    if(status == 0)
	exit(0);
    sprintf(outputFilePath, "%s.sz", oriFilePath);
   
    size_t nbEle;
    float *data = readFloatData(oriFilePath, &nbEle, &status);
    /*size_t i, j, index = 0;
    for(i=0;i<5;i++)
    {
	index = i*1000000;
    	for(j=0;j<1000;j++)
	{
		printf("%f\n", data[index+j]);
	}
	printf("===========%zu\n", index);
    }*/
    if(status != SZ_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
    //float *revValue = (float *)malloc(sizeof(float));
    //*revValue = 1.0E36;
   
    size_t outSize; 
    //char *bytes = (char *)malloc(nbEle*sizeof(float)); //
    //SZ_compress_args2(SZ_FLOAT, data, bytes, &outSize, ABS, 0.0001, 0.0001, r5, r4, r3, r2, r1);    
    //char *bytes = SZ_compress_rev(SZ_FLOAT, data, revValue, &outSize, r5, r4, r3, r2, r1);
    cost_start();
    unsigned char *bytes = SZ_compress(SZ_FLOAT, data, &outSize, r5, r4, r3, r2, r1);
    cost_end();
    printf("timecost=%f\n",totalCost); 
    writeByteData(bytes, outSize, outputFilePath, &status);
    if(status != SZ_SCES)
    {
        printf("Error: data file %s cannot be written!\n", outputFilePath);
        exit(0);
    }

   /* //check mem leakage of decompression
    int i;
    float *ddata;
    for(i=0;i<30;i++)
    {
	printf("start %d\n",i);
	ddata = SZ_decompress(SZ_FLOAT, bytes, outSize, r5, r4, r3, r2, r1);
	free(ddata);
	printf("end %d\n", i);
	sleep(1);
    }
*/
    printf("done\n");
    free(bytes); 
    free(data);
    SZ_Finalize();
    
    return 0;
}
