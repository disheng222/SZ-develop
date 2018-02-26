#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sz.h"
#include "rw.h"

int main(){

	// size_t nbEle = 134217728;
	// nbEle = 25000000;
	size_t nbEle;
	int status = 0;
	char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/dark_matter_density.raw";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/temperature.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/baryon_density.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/velocity_y.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/QCLOUDf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/Uf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/CLOUDf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/Pf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/TCf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/Wf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/PRECIPf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/QGRAUPf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/QICEf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/qmcpack.dat";
	float *data = readFloatData(path, &nbEle, &status);
	SZ_Init("sz.config");
	printf("read data done, size: %ld\n", nbEle);
	// nbEle = 288*115*69*69;
	size_t comped_size;
	// unsigned char *comp_data = SZ_compress_args(SZ_FLOAT, data, &comped_size, ABS, 1, 1, 1, 0, 0, 0, 0, 0, 134217728);
	// unsigned char *comp_data = SZ_compress_args(SZ_FLOAT, data, &comped_size, ABS, 4780.302, 1, 1, 0, 0, 0, 512, 512, 512);
	// unsigned char *comp_data = SZ_compress_args(SZ_FLOAT, data, &comped_size, ABS, 0.01056972, 1, 1, 0, 0, 0, 100, 500, 500);
	// printf("quantization bins: %d\n", intvCapacity);
	// unsigned char * comp_data = SZ_compress_float_1D_MDQ_RA(data, 25000000, 1e-2, &comped_size);
	// unsigned char * comp_data = SZ_compress_float_3D_MDQ_RA(data, 100, 500, 500, 0.001056972, &comped_size);
	// unsigned char * comp_data = SZ_compress_float_3D_MDQ_RA(data, 512, 512, 512, 0.55, &comped_size);
	// unsigned char * comp_data = SZ_compress_float_3D_MDQ_nonblocked(data, 100, 500, 500, 0.001056972, &comped_size);
	unsigned char * comp_data = SZ_compress_float_3D_MDQ_nonblocked(data, 512, 512, 512, 0.48, &comped_size);

	// printf("comped size: %ld\n", comped_size);
	unsigned char * out;
	unsigned long size = zlib_compress5(comp_data, comped_size, &out, 1);
	printf("%ld\n", size);
	free(comp_data);

	unsigned long tmpSize = zlib_uncompress5(out, size, &comp_data, comped_size);
	printf("%ld\n", tmpSize);
	// exit(0);


	float * result;
	// result = SZ_decompress(SZ_FLOAT, comp_data, comped_size, 0, 0, 0, 0, 134217728);
	// result = SZ_decompress(SZ_FLOAT, comp_data, comped_size, 0, 0, 512, 512, 512);
	// result = SZ_decompress(SZ_FLOAT, comp_data, comped_size, 0, 0, 100, 500, 500);
	// decompressDataSeries_float_1D_RA(&result, 25000000, comp_data + 24);
	// decompressDataSeries_float_3D_RA(&result, 100, 500, 500, comp_data + 24);
	// decompressDataSeries_float_3D_RA(&result, 512, 512, 512, comp_data + 24);
	// decompressDataSeries_float_3D_nonblocked(&result, 100, 500, 500, comp_data + 24);
	decompressDataSeries_float_3D_nonblocked(&result, 512, 512, 512, comp_data + 24);

	printf("decompression done, nbEle: %ld\n", nbEle);
	fflush(stdout);

	{
		size_t byteLength = size;
		float * ori_data = data;
		float * data = result;
		size_t i = 0;
	    float Max = 0, Min = 0, diffMax = 0;
	    Max = ori_data[0];
	    Min = ori_data[0];
	    diffMax = fabs(data[0] - ori_data[0]);
	    size_t k = 0;
	    double sum1 = 0, sum2 = 0;
	    for (i = 0; i < nbEle; i++)
	    {
	        sum1 += ori_data[i];
			sum2 += data[i];
	    }
	    double mean1 = sum1/nbEle;
	    double mean2 = sum2/nbEle;

	    double sum3 = 0, sum4 = 0;
	    double sum = 0, prodSum = 0, relerr = 0;
	   
	    double maxpw_relerr = 0; 
	    for (i = 0; i < nbEle; i++)
	    {
	        if (Max < ori_data[i])
	        	Max = ori_data[i];
	        if (Min > ori_data[i])
	        	Min = ori_data[i];
	        
	        float err = fabs(data[i] - ori_data[i]);

			if(ori_data[i]!=0)
			{
				if(fabs(ori_data[i])>1)
					relerr = err/ori_data[i];
				else
					relerr = err;
				if(maxpw_relerr<relerr)
					maxpw_relerr = relerr;
	        }

			if (diffMax < err)
				diffMax = err;
	        prodSum += (ori_data[i]-mean1)*(data[i]-mean2);
	        sum3 += (ori_data[i] - mean1)*(ori_data[i]-mean1);
	        sum4 += (data[i] - mean2)*(data[i]-mean2);
			sum += err*err;	
	    }
	    double std1 = sqrt(sum3/nbEle);
	    double std2 = sqrt(sum4/nbEle);
	    double ee = prodSum/nbEle;
	    double acEff = ee/std1/std2;
	 
	    double mse = sum/nbEle;
	    double range = Max - Min;
	    double psnr = 20*log10(range)-10*log10(mse);
	    double nrmse = sqrt(mse)/range;
	     
	    double compressionRatio = 1.0*nbEle*sizeof(float)/byteLength;

	    printf ("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
	    printf ("Max absolute error = %.10f\n", diffMax);
	    printf ("Max relative error = %f\n", diffMax/(Max-Min));
	    printf ("Max pw relative error = %f\n", maxpw_relerr);
	    printf ("PSNR = %f, NRMSE= %.20G\n", psnr,nrmse);
	    printf ("acEff=%f\n", acEff);
	    printf ("compressionRatio = %f\n", compressionRatio);
	}


	printf("write decompressed data to file\n");
	fflush(stdout);
	// writeFloatData_inBytes(result, nbEle, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/decompressed.dat", &status);
	free(out);
	free(data);
	free(result);
	free(comp_data);
	SZ_Finalize();
	
}
