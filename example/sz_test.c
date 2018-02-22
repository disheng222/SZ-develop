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
	//char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/dark_matter_density.raw";
	char path[100] = "/home/sdi/Data/NYX/dark_matter_density.raw";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/QCLOUDf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/qmcpack.dat";
	float *data = readFloatData(path, &nbEle, &status);
	SZ_Init("sz.config");
	SZ_Reset();

	// TightDataPointStorageRA * tdps_ra = SZ_compress_float_3D_MDQ_RA(data, 100, 500, 500, 4e-6);
	// nbEle = 288*115*69*69;
	// TightDataPointStorageRA * tdps_ra = SZ_compress_float_3D_MDQ_RA(data, 33120, 69, 69, 0.03301001489);
	// TightDataPointStorageRA * tdps_ra = SZ_compress_float_3D_MDQ(data, 33120, 69, 69, 0.01);
	size_t comped_size;
	 unsigned char * comp_data = SZ_compress_float_1D_MDQ_RA(data, 134217728, 1, &comped_size);
	//unsigned char * comp_data = SZ_compress_float_3D_MDQ_RA(data, 512, 512, 512, 1, &comped_size);
	//unsigned char * comp_data = SZ_compress_float_3D_MDQ_RA(data, 100, 500, 500, 1e-4, &comped_size);

	unsigned char * out;
	unsigned long size = zlib_compress5(comp_data, comped_size, &out, gzipMode);
	printf("%ld\n", size);

	SZ_Reset();
	free(comp_data);
	unsigned long tmpSize = zlib_uncompress5(out, size, &comp_data, comped_size);
	printf("%ld\n", tmpSize);


	float * result;
	 decompressDataSeries_float_1D_RA(&result, 134217728, comp_data);
	//decompressDataSeries_float_3D_RA(&result, 512, 512, 512, comp_data);
	//decompressDataSeries_float_3D_RA(&result, 100, 500, 500, comp_data);
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
	        if (Max < ori_data[i]) Max = ori_data[i];
	        if (Min > ori_data[i]) Min = ori_data[i];
	        
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
