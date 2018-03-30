#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sz.h"
#include "rw.h"
#include "sz_omp.h"
#include <time.h>

#define SZ_OMP -6
#define SZ_REG_NONBLOCKED -5
#define SZ_REG -4
#define REG -3
#define ADAPTIVE -2
#define ORI -1
#define BLOCKED 0
#define NONBLOCKED 1
#define BLOCKED_MULTIMEAN 2
#define NONBLOCKED_MULTIMEAN 3
#define HC 0
#define NYX 1

int main(int argc, char ** argv){

	// size_t nbEle = 134217728;
	// nbEle = 25000000;
	size_t nbEle;
	int status = 0;
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/dark_matter_density.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/dark_matter_density.dat.log10";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/temperature.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/temperaturelog10.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/baryon_density.dat.log10";
	char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/velocity_x.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/velocity_y.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/velocity_z.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/gokhan/u-130x640x641.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/gokhan/v-130x641x640.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/Hurricane/step48/CLOUDf48.bin.dat.log10";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/Hurricane/step48/QCLOUDf48.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/Hurricane/step48/QSNOWf48.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/Hurricane/QCLOUDf48log10.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/Hurricane/step48/TCf48.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/CLOUDf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/Pf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/TCf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/Wf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/PRECIPf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/QGRAUPf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/QICEf01.bin.dat";
	// char path[100] = "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/qmcpack.dat";
	float *data = readFloatData(path, &nbEle, &status);
	SZ_Init("../sz.config");
	printf("read data done, size: %ld\n", nbEle);
	// nbEle = 288*115*69*69;
	size_t comped_size;
	double HC_errorBound = 0.000001056972;
	double NYX_errorBound = 0.02;

	int verbose = 1;
	// int dim = 1;
	// int dim = 2;
	// int dim = 3;
	// int mode = BLOCKED;
	// int mode = NONBLOCKED;
	// int mode = BLOCKED_MULTIMEAN;
	// int mode = NONBLOCKED_MULTIMEAN;
	// int mode = ORI;
	// int mode = ADAPTIVE;
	// int mode = REG;
	// int mode = SZ_REG;
	// int mode = SZ_REG_NONBLOCKED;
	int mode = SZ_OMP;
	// int dataset = HC;
	int dataset = NYX;

	unsigned char * comp_data;
	size_t n1, n2, n3;
	// n1 = 130;
	// n2 = 640;
	// n3 = 641;
	// n1 = 100;
	// n2 = 500;
	// n3 = 500;
	n1 = 512;
	n2 = 512;
	n3 = 512;

	// double GH_errorBound = 1.4216423988342285156*1.1;
	// double GH_errorBound = 0.22826410293579101562;
	// double GH_errorBound = 0.106201141357421875 * 0.1; // TCf48
	// double GH_errorBound = 0.137789345703125 * 5; // dark matter
	// double GH_errorBound = 478.03025*5; // temperature
	// double GH_errorBound = 11.5862234375; // baryon
	double GH_errorBound = 82283.640; // velocity_x
	// double GH_errorBound = 100438.536; // velocity_y
	// double GH_errorBound = 723242.40 * 0.01; // velocity_z
	// double GH_errorBound = 0.00000078183296136558055878;
	// double GH_errorBound = 0.0000020479534287005662918; // QCLOUDf48

	double elapsed_time = 0.0;
	elapsed_time = -omp_get_wtime();
	switch(mode){
		case BLOCKED: 	comp_data = SZ_compress_float_3D_MDQ_RA(data, n1, n2, n3, GH_errorBound, &comped_size);
						break;
		case NONBLOCKED:	comp_data = SZ_compress_float_3D_MDQ_nonblocked(data, n1, n2, n3, GH_errorBound, &comped_size);
							break;
		case NONBLOCKED_MULTIMEAN:	comp_data = SZ_compress_float_3D_MDQ_nonblocked_multi_means(data, n1, n2, n3, GH_errorBound, &comped_size);
									break;
		case ORI:	comp_data = SZ_compress_float_3D_MDQ_nonblocked_ori(data, n1, n2, n3, GH_errorBound, &comped_size);
					break;
		case ADAPTIVE:	comp_data = SZ_compress_float_3D_MDQ_nonblocked_adaptive(data, n1, n2, n3, GH_errorBound, &comped_size);
						break;
		case REG:	comp_data = SZ_compress_float_3D_MDQ_RA_all_by_regression(data, n1, n2, n3, GH_errorBound, &comped_size);
					break;
		case SZ_REG:	comp_data = SZ_compress_float_3D_MDQ_RA_blocked_with_regression(data, n1, n2, n3, GH_errorBound, &comped_size);
						break;
		case SZ_REG_NONBLOCKED:	comp_data = SZ_compress_float_3D_MDQ_nonblocked_with_blocked_regression(data, n1, n2, n3, GH_errorBound, &comped_size);
								break;
		case SZ_OMP:	comp_data = SZ_compress_float_3D_MDQ_openmp(data, n1, n2, n3, GH_errorBound, &comped_size);
						break;
	}
	elapsed_time += omp_get_wtime();
	printf("compression elapsed time: %.4f\n", elapsed_time);
	// exit(0);

	elapsed_time = -omp_get_wtime();
	unsigned char * out;
	unsigned long size = zlib_compress5(comp_data, comped_size, &out, 1);
	elapsed_time += omp_get_wtime();
	printf("gzip elapsed_time: %.4f\n", elapsed_time);
	free(comp_data);

	elapsed_time = -omp_get_wtime();
	unsigned long tmpSize = zlib_uncompress5(out, size, &comp_data, comped_size);
	elapsed_time += omp_get_wtime();
	printf("gzip decompress elapsed_time: %.4f\n", elapsed_time);
	// exit(0);
	float * result;
	elapsed_time = -omp_get_wtime();
	switch(mode){
		case BLOCKED: 	decompressDataSeries_float_3D_RA(&result, n1, n2, n3, comp_data + 24);
						break;
		case NONBLOCKED:	decompressDataSeries_float_3D_nonblocked(&result, n1, n2, n3, comp_data + 24);
							break;
		case NONBLOCKED_MULTIMEAN:	decompressDataSeries_float_3D_nonblocked_multi_means(&result, n1, n2, n3, comp_data + 24);
									break;
		case ORI:	decompressDataSeries_float_3D_nonblocked_ori(&result, n1, n2, n3, comp_data + 24);
					break;
		case ADAPTIVE:	decompressDataSeries_float_3D_nonblocked_adaptive(&result, n1, n2, n3, comp_data + 24);
						break;
		case REG:	decompressDataSeries_float_3D_RA_all_by_regression(&result, n1, n2, n3, comp_data + 24);
					break;
		case SZ_REG:	decompressDataSeries_float_3D_RA_blocked_with_regression(&result, n1, n2, n3, comp_data + 24);
						break;
		case SZ_REG_NONBLOCKED:	decompressDataSeries_float_3D_nonblocked_with_blocked_regression(&result, n1, n2, n3, comp_data + 24);
								break;
		case SZ_OMP:	decompressDataSeries_float_3D_openmp(&result, n1, n2, n3, comp_data + 24);
						break;
	}
	elapsed_time += omp_get_wtime();
	printf("decompression elapsed time: %.4f\n", elapsed_time);
	printf("decompression done, nbEle: %ld\n", nbEle);
	fflush(stdout);

	if(verbose)
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
        unsigned int tmp_count = 0;
	    for (i = 0; i < nbEle; i++)
	    {
	        if (Max < ori_data[i])
	        	Max = ori_data[i];
	        if (Min > ori_data[i])
	        	Min = ori_data[i];
	        
	        float err = fabs(data[i] - ori_data[i]);
	        if(err > 0.83){
	        	tmp_count ++;
	        }

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
        printf("tmp count %d\n", tmp_count);
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


	// printf("write decompressed data to file\n");
	// fflush(stdout);
	// {
	// 	float min = 1000000;
	// 	for(size_t i=0; i<nbEle; i++){
	// 		if(result[i] > 0 && result[i] < min){
	// 			min = result[i];
	// 		}
	// 	}
	// 	printf("min: %.4f\n", min);
	// 	for(size_t i=0; i<nbEle; i++){
	// 		if(result[i] > 0){
	// 			result[i] = log10(result[i]);
	// 		}
	// 		else result[i] = log10(min);
	// 	}

	// writeFloatData_inBytes(result, nbEle, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/dataimage-generator/dark_matter_density_log10_decompressed.dat", &status);
	// }
	// free(out);
	free(data);
	free(result);
	free(comp_data);
	SZ_Finalize();
	
}
