/**
 *  @file sz_float.c
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
 *  @brief SZ_Init, Compression and Decompression functions
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "sz.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "TightDataPointStorageF.h"
#include "sz_float.h"
#include "sz_float_pwr.h"
#include "szd_float.h"
#include "szd_float_pwr.h"
#include "zlib.h"
#include "rw.h"

// descend sorting
int cmp(const void * e1, const void * e2) 
{
    size_t f = *((size_t*)e1);
    size_t s = *((size_t*)e2);
    if (f > s) return -1;
    if (f < s) return 1;
    return 0;
}

unsigned int optimize_intervals_float_1D(float *oriData, size_t dataLength, double realPrecision)
{	
	size_t i = 0, radiusIndex;
	float pred_value = 0, pred_err;
	size_t *intervals = (size_t*)malloc(maxRangeRadius*sizeof(size_t));
	memset(intervals, 0, maxRangeRadius*sizeof(size_t));
	size_t totalSampleSize = dataLength/sampleDistance;
	for(i=2;i<dataLength;i++)
	{
		if(i%sampleDistance==0)
		{
			//pred_value = 2*oriData[i-1] - oriData[i-2];
			pred_value = oriData[i-1];
			pred_err = fabs(pred_value - oriData[i]);
			radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
			if(radiusIndex>=maxRangeRadius)
				radiusIndex = maxRangeRadius - 1;			
			intervals[radiusIndex]++;
		}
	}
	//compute the appropriate number
	size_t targetCount = totalSampleSize*predThreshold;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
		
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);
	
	if(powerOf2<32)
		powerOf2 = 32;
	
	free(intervals);
	//printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);
	return powerOf2;
}

unsigned int optimize_intervals_and_compute_dense_position_float_1D(float *oriData, size_t dataLength, double realPrecision, float * dense_pos){

	// compute mean
	float mean = 0.0;
	size_t mean_distance = (int) (sqrt(dataLength));
	size_t mean_sample_size = (int) (dataLength / mean_distance - 1);
	float * data_pos = oriData;
	for(size_t i=0; i<mean_sample_size; i++){
		mean += *data_pos;
		data_pos += mean_distance;
	}
	if(mean_sample_size > 0) mean /= mean_sample_size;
	size_t range = 256;
	size_t radius = 128;
	size_t * freq_intervals = (size_t *) malloc(range*sizeof(size_t));
	memset(freq_intervals, 0, range*sizeof(size_t));
	// compute optimized intervals
	size_t i = 0, radiusIndex;
	float pred_value = 0, pred_err;
	size_t *intervals = (size_t*)malloc(maxRangeRadius*sizeof(size_t));
	memset(intervals, 0, maxRangeRadius*sizeof(size_t));
	size_t totalSampleSize = (dataLength - 1)/sampleDistance;

	float mean_diff;
	ptrdiff_t freq_index;
	data_pos = oriData + 1;
	for(size_t i=0; i<totalSampleSize; i++){
		// optimize interval
		pred_value = data_pos[-1];
		pred_err = fabs(pred_value - oriData[i]);
		radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
		if(radiusIndex>=maxRangeRadius)
			radiusIndex = maxRangeRadius - 1;			
		intervals[radiusIndex]++;
		// collect frequency
		mean_diff = data_pos[0] - mean;
		freq_index = (ptrdiff_t)(mean_diff/realPrecision) + radius;
		if(freq_index <= 0){
			freq_intervals[0] ++;
		}
		else if(freq_index >= range){
			freq_intervals[range - 1] ++;
		}
		else{
			freq_intervals[freq_index] ++;
		}
		data_pos += sampleDistance;
	}

	//compute the appropriate number
	size_t targetCount = totalSampleSize*predThreshold*1.009;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
		
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);
	
	if(powerOf2<32)
		powerOf2 = 32;
	
	free(intervals);
	//printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);

	// compute estimated dense position
	size_t max_sum = 0;
	size_t max_index = 0;
	size_t tmp_sum;
	size_t * freq_pos = freq_intervals;
	for(size_t i=1; i<range; i++){
		tmp_sum = freq_pos[0] + freq_pos[1];
		if(tmp_sum > max_sum){
			max_sum = tmp_sum;
			max_index = i;
		}
		freq_pos ++;
	}
	*dense_pos = mean + realPrecision * (ptrdiff_t)(max_index - 1 - radius);
	// printf("real precision: %.4f dense_pos: %.4f\n", realPrecision, dense_pos[0]);
	free(freq_intervals);
	return powerOf2;

}

unsigned int optimize_intervals_float_2D(float *oriData, size_t r1, size_t r2, double realPrecision)
{	
	size_t i,j, index;
	size_t radiusIndex;
	float pred_value = 0, pred_err;
	size_t *intervals = (size_t*)malloc(maxRangeRadius*sizeof(size_t));
	memset(intervals, 0, maxRangeRadius*sizeof(size_t));
	size_t totalSampleSize = (r1-1)*(r2-1)/sampleDistance;

	//float max = oriData[0];
	//float min = oriData[0];

	for(i=1;i<r1;i++)
	{
		for(j=1;j<r2;j++)
		{
			if((i+j)%sampleDistance==0)
			{
				index = i*r2+j;
				pred_value = oriData[index-1] + oriData[index-r2] - oriData[index-r2-1];
				pred_err = fabs(pred_value - oriData[index]);
				radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
				if(radiusIndex>=maxRangeRadius)
					radiusIndex = maxRangeRadius - 1;
				intervals[radiusIndex]++;

			//	if (max < oriData[index]) max = oriData[index];
			//	if (min > oriData[index]) min = oriData[index];
			}			
		}
	}
	//compute the appropriate number
	size_t targetCount = totalSampleSize*predThreshold;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	//	struct timeval costStart, costEnd;
	//	double cost_est = 0;
	//
	//	gettimeofday(&costStart, NULL);
	//
	//	//compute estimate of bit-rate and distortion
	//	double est_br = 0;
	//	double est_psnr = 0;
	//	double c1 = log2(targetCount)+1;
	//	double c2 = -20.0*log10(realPrecision) + 20.0*log10(max-min) + 10.0*log10(3);
	//
	//	for (i = 0; i < powerOf2/2; i++)
	//	{
	//		int count = intervals[i];
	//		if (count != 0)
	//			est_br += count*log2(count);
	//		est_psnr += count;
	//	}
	//
	//	//compute estimate of bit-rate
	//	est_br -= c1*est_psnr;
	//	est_br /= totalSampleSize;
	//	est_br = -est_br;
	//
	//	//compute estimate of psnr
	//	est_psnr /= totalSampleSize;
	//	printf ("sum of P(i) = %lf\n", est_psnr);
	//	est_psnr = -10.0*log10(est_psnr);
	//	est_psnr += c2;
	//
	//	printf ("estimate bitrate = %.2f\n", est_br);
	//	printf ("estimate psnr = %.2f\n",est_psnr);
	//
	//	gettimeofday(&costEnd, NULL);
	//	cost_est = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
	//
	//	printf ("analysis time = %f\n", cost_est);

	free(intervals);
	//printf("maxRangeRadius = %d, accIntervals=%d, powerOf2=%d\n", maxRangeRadius, accIntervals, powerOf2);
	return powerOf2;
}

unsigned int optimize_intervals_float_3D(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision)
{	
	size_t i,j,k, index;
	size_t radiusIndex;
	size_t r23=r2*r3;
	float pred_value = 0, pred_err;
	size_t *intervals = (size_t*)malloc(maxRangeRadius*sizeof(size_t));
	memset(intervals, 0, maxRangeRadius*sizeof(size_t));
	size_t totalSampleSize = (r1-1)*(r2-1)*(r3-1)/sampleDistance;

	//float max = oriData[0];
	//float min = oriData[0];

	for(i=1;i<r1;i++)
	{
		for(j=1;j<r2;j++)
		{
			for(k=1;k<r3;k++)
			{			
				if((i+j+k)%sampleDistance==0)
				{
					index = i*r23+j*r3+k;
					pred_value = oriData[index-1] + oriData[index-r3] + oriData[index-r23] 
					- oriData[index-1-r23] - oriData[index-r3-1] - oriData[index-r3-r23] + oriData[index-r3-r23-1];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (pred_err/realPrecision+1)/2;
					if(radiusIndex>=maxRangeRadius)
					{
						radiusIndex = maxRangeRadius - 1;
						//printf("radiusIndex=%d\n", radiusIndex);
					}
					intervals[radiusIndex]++;

					//	if (max < oriData[index]) max = oriData[index];
					//	if (min > oriData[index]) min = oriData[index];
				}
			}
		}
	}
	//compute the appropriate number
	size_t targetCount = totalSampleSize*predThreshold;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;
	
	//	struct timeval costStart, costEnd;
	//	double cost_est = 0;
	//
	//	gettimeofday(&costStart, NULL);
	//
	//	//compute estimate of bit-rate and distortion
	//	double est_br = 0;
	//	double est_psnr = 0;
	//	double c1 = log2(targetCount)+1;
	//	double c2 = -20.0*log10(realPrecision) + 20.0*log10(max-min) + 10.0*log10(3);
	//
	//	for (i = 0; i < powerOf2/2; i++)
	//	{
	//		int count = intervals[i];
	//		if (count != 0)
	//			est_br += count*log2(count);
	//		est_psnr += count;
	//	}
	//
	//	//compute estimate of bit-rate
	//	est_br -= c1*est_psnr;
	//	est_br /= totalSampleSize;
	//	est_br = -est_br;
	//
	//	//compute estimate of psnr
	//	est_psnr /= totalSampleSize;
	//	printf ("sum of P(i) = %lf\n", est_psnr);
	//	est_psnr = -10.0*log10(est_psnr);
	//	est_psnr += c2;
	//
	//	printf ("estimate bitrate = %.2f\n", est_br);
	//	printf ("estimate psnr = %.2f\n",est_psnr);
	//
	//	gettimeofday(&costEnd, NULL);
	//	cost_est = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
	//
	//	printf ("analysis time = %f\n", cost_est);

	free(intervals);
	//printf("targetCount=%d, sum=%d, totalSampleSize=%d, ratio=%f, accIntervals=%d, powerOf2=%d\n", targetCount, sum, totalSampleSize, (double)sum/(double)totalSampleSize, accIntervals, powerOf2);
	return powerOf2;
}

unsigned int optimize_intervals_and_compute_dense_position_float_3D(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, float * dense_pos){	
	// compute mean
	float mean = 0.0;
	size_t len = r1 * r2 * r3;
	size_t mean_distance = (int) (sqrt(len));
	// printf("mean_distance %ld, mean_sample_size %ld\n", mean_distance, mean_sample_size);
	float * data_pos = oriData;
	size_t offset_count = 0;
	size_t offset_count_2 = 0;
	size_t mean_count = 0;
	while(data_pos - oriData < len){
		mean += *data_pos;
		mean_count ++;
		data_pos += mean_distance;
		offset_count += mean_distance;
		offset_count_2 += mean_distance;
		if(offset_count >= r3){
			offset_count = 0;
			data_pos -= 1;
		}
		if(offset_count_2 >= r2 * r3){
			offset_count_2 = 0;
			data_pos -= 1;
		}
	}
	if(mean_count > 0) mean /= mean_count;
	size_t range = 8192;
	size_t radius = 4096;
	size_t * freq_intervals = (size_t *) malloc(range*sizeof(size_t));
	memset(freq_intervals, 0, range*sizeof(size_t));
	float * sum_intervals = (float *) malloc(range*sizeof(float));
	memset(sum_intervals, 0, range*sizeof(float));
	float value_sum = 0;

	size_t i,j,k, index;
	size_t radiusIndex;
	size_t r23=r2*r3;
	float pred_value = 0, pred_err;
	size_t *intervals = (size_t*)malloc(maxRangeRadius*sizeof(size_t));
	memset(intervals, 0, maxRangeRadius*sizeof(size_t));
	size_t *intervals_2D = (size_t*)malloc(maxRangeRadius*sizeof(size_t));
	memset(intervals_2D, 0, maxRangeRadius*sizeof(size_t));

	float mean_diff;
	ptrdiff_t freq_index;
	size_t totalSampleSize = 0;
	for(i=1;i<r1;i++)
	{
		for(j=1;j<r2;j++)
		{
			for(k=1;k<r3;k++)
			{			
				if((i+j+k)%sampleDistance==0)
				{
					index = i*r23+j*r3+k;
					mean_diff = oriData[index] - mean;
					if(mean_diff > 0) freq_index = (ptrdiff_t)(mean_diff/realPrecision) + radius;
					else freq_index = (ptrdiff_t)(mean_diff/realPrecision) - 1 + radius;
					if(freq_index <= 0){
						freq_intervals[0] ++;
					}
					else if(freq_index >= range){
						freq_intervals[range - 1] ++;
					}
					else{
						freq_intervals[freq_index] ++;
						sum_intervals[freq_index] += oriData[index];
						value_sum += oriData[index];
					}
					totalSampleSize ++;
				}
			}
		}
	}

	// collect frequency
	size_t max_sum = 0;
	size_t max_index = 0;
	size_t tmp_sum;
	size_t * freq_pos = freq_intervals + 1;
	for(size_t i=1; i<range-2; i++){
		tmp_sum = freq_pos[0] + freq_pos[1];
		if(tmp_sum > max_sum){
			max_sum = tmp_sum;
			max_index = i;
		}
		freq_pos ++;
	}
	*dense_pos = mean + realPrecision * (ptrdiff_t)(max_index + 1 - radius);
	printf("Max frequency: %.6f index: %ld dense_pos: %.6f\n", max_sum * 1.0 / totalSampleSize, max_index, *dense_pos);

	// compute intervals
	float dense_point = *dense_pos;
	size_t unflushed_count = 0;

	// debug
	int * typeArray = (int *) malloc(totalSampleSize*sizeof(int));
	size_t typeArray_count = 0;
	int * typeArray_flag = (int *) malloc(totalSampleSize*sizeof(int));
	memset(typeArray, 0, maxRangeRadius*sizeof(int));

	for(i=1;i<r1;i++)
	{
		for(j=1;j<r2;j++)
		{
			for(k=1;k<r3;k++)
			{			
				if((i+j+k)%sampleDistance==0)
				{
					// 3D prediction
					index = i*r23+j*r3+k;
					if(fabs(oriData[index] - dense_point) <= realPrecision){
						typeArray[typeArray_count ++] = 513;
						continue;
					}
					pred_value = oriData[index-1] + oriData[index-r3] + oriData[index-r23] 
					- oriData[index-1-r23] - oriData[index-r3-1] - oriData[index-r3-r23] + oriData[index-r3-r23-1];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (pred_err/realPrecision+1)/2;
					if(radiusIndex>=maxRangeRadius){
						radiusIndex = maxRangeRadius - 1;
						typeArray[typeArray_count ++] = 514;
					}
					else{
						if(typeArray_flag[radiusIndex] > 0){
							typeArray[typeArray_count ++] = 2 * radiusIndex;
							typeArray_flag[radiusIndex] = 0;
						}
						else{
							typeArray[typeArray_count ++] = 2 * radiusIndex + 1;
							typeArray_flag[radiusIndex] = 1;
						}
					}
					intervals[radiusIndex]++;

					// 2D prediction
					pred_value = oriData[index-1] + oriData[index-r3] - oriData[index-1-r3];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (pred_err/realPrecision+1)/2;
					if(radiusIndex>=maxRangeRadius){
						radiusIndex = maxRangeRadius - 1;
					}
					intervals_2D[radiusIndex]++;

					// 2D another direction
					pred_value = oriData[index-1] + oriData[index-r23] - oriData[index-r23-1];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (pred_err/realPrecision+1)/2;
					if(radiusIndex>=maxRangeRadius){
						radiusIndex = maxRangeRadius - 1;
					}
					intervals_2D[radiusIndex]++;

					// 2D 3rd direction
					pred_value = oriData[index-r3] + oriData[index-r23] - oriData[index-r23-r3];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (pred_err/realPrecision+1)/2;
					if(radiusIndex>=maxRangeRadius){
						radiusIndex = maxRangeRadius - 1;
					}
					intervals_2D[radiusIndex]++;

					unflushed_count ++;
				}
			}
		}
	}
	printf("unflushed_count: %ld with respect to total sample count %ld\n", unflushed_count, totalSampleSize);
	//compute the appropriate number
	size_t targetCount = unflushed_count*predThreshold;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	printf("Estimated powerOf2: %ld, sum %ld percent %.4f\n", powerOf2, sum, sum * 1.0 / unflushed_count);
	// use block or not
	// qsort(intervals, maxRangeRadius, sizeof(size_t), cmp);
	// qsort(intervals_2D, maxRangeRadius, sizeof(size_t), cmp);

	unsigned int block_opt[5] = {4, 8, 12, 16, 20};
	float block_3D_est_prob[5];
	for(i=0; i<5; i++){
		float tmp = 1 - 1.0 / block_opt[i]; 
		block_3D_est_prob[i] = tmp * tmp * tmp;
	}
	float block_br[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	float inf_br = 0.0;
	float p2D_br = 0.0;
	double block_prob;
	double block_accum[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	double probs_3D, probs_2D;
	size_t accum_2D = 0;
	size_t accum_3D = 0;
	for(i=0; i<powerOf2/2; i++){
		probs_3D = intervals[i] * 1.0 / 2 / totalSampleSize;
		// printf("%d: %.6f\n", i, probs_3D);
		accum_3D += intervals[i];
		// INTERVALS_2D already contains number of 3 directions
		probs_2D = intervals_2D[i] * 1.0 / 2 / (3 * totalSampleSize);
		accum_2D += intervals_2D[i];
		for(int j=0; j<5; j++){
			block_prob = (block_3D_est_prob[j] * probs_3D + (1 - block_3D_est_prob[j]) * probs_2D);
			if(block_prob > 0){
				block_br[j] += - block_prob * log2(block_prob);
				block_accum[j] += block_prob;
			}
		}
		if(probs_3D > 0){
			inf_br += - probs_3D * log2(probs_3D);
			accum_3D += probs_3D;
		}
		if(probs_2D > 0){
			p2D_br += - probs_2D * log2(probs_2D);
			accum_2D += probs_2D;
		}
	}
	printf("accum count: %ld %ld\n", accum_3D, accum_2D);
	{
		//huffman encoding result
		SZ_Reset(allNodes, stateNum);
		unsigned char * tmp = (unsigned char *)malloc(totalSampleSize*sizeof(int));
		size_t enCodeSize = 0;
		encode_withTree(typeArray, totalSampleSize, &tmp, &enCodeSize);
		printf("Huffman bitrate: %.4f\n", enCodeSize * 32.0 / (totalSampleSize * sizeof(int)));
		free(tmp);
		// int status;
		// writeIntData_inBytes(typeArray, totalSampleSize, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/type_array.dat", &status);
		// writeIntData_inBytes(intervals, maxRangeRadius, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/intervals.dat", &status);
		free(typeArray);
	}

	double flushed_prob = 1 - unflushed_count * 1.0 / totalSampleSize; 
	float flushed_br = - flushed_prob * log2(flushed_prob);
	printf("Flushed prob: %.4f\n", flushed_prob);
	double unpred_prob;
	float unpred_br;
	double unpred_3D = (unflushed_count - accum_3D) * 1.0 / totalSampleSize;
	double unpred_2D = (unflushed_count * 3 - accum_2D) * 1.0 / (3 * totalSampleSize);
	printf("Unpredicted prob: ");
	for(int j=0; j<5; j++){
		unpred_prob = (block_3D_est_prob[j] * unpred_3D + (1 - block_3D_est_prob[j]) * unpred_2D);
		if(unpred_prob > 0) unpred_br = unpred_prob * log2(unpred_prob);
		else{
			unpred_prob = 0;
			unpred_br = 0;
		}
		printf("%.6f ", unpred_prob);		
		//				2 p log p  ---------------		unpred 									block overhead							
		block_br[j] = 2 * block_br[j] + flushed_br + unpred_br;// + unpred_prob * 32;// + 40.0 / (block_opt[j] * block_opt[j] * block_opt[j]);
	}
	printf("%.6f %.6f\n", unpred_3D, unpred_2D);
	unpred_prob = unpred_3D;
	unpred_br = unpred_prob * log2(unpred_prob);
	if(unpred_prob > 0) unpred_br = unpred_prob * log2(unpred_prob);
	else{
		unpred_prob = 0;
		unpred_br = 0;
	}
	inf_br = 2 * inf_br + flushed_br + unpred_br;// + unpred_prob * 32;
	unpred_prob = unpred_2D;
	unpred_br = unpred_prob * log2(unpred_prob);
	if(unpred_prob > 0) unpred_br = unpred_prob * log2(unpred_prob);
	else{
		unpred_prob = 0;
		unpred_br = 0;
	}
	p2D_br = 2 * p2D_br + flushed_br + unpred_br;// + unpred_prob * 32;

	printf("Blocked bit cost: ");
	for(int j=0; j<5; j++){
		printf("%.6f ", block_br[j]);
	}
	printf("\nNonblocked bit cost: 3D %.6f\n2D %.6f\n", inf_br, p2D_br);

	printf("Overhead over 3D:\n");
	for(int j=0; j<5; j++){
		printf("%.4f%% ", (block_br[j] / inf_br - 1)*100);
	}
	printf("\n");
	// TODO: estimate mse and PSNR???
	// exit(0);

	// float block_cost[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	// for(int j=0; j<5; j++){
	// 	block_cost[j] = block_br[j] + 0;
	// }
	// printf("Overall cost: ");
	// for(int j=0; j<5; j++){
	// 	printf("%.6f ", block_cost[j]);
	// }
	// printf("%.6f %.6f\n", inf_br, p2D_br);


	size_t less_than_mean_freq_count = 0;
	for(i=0; i<accIntervals/2; i++){
		if(max_sum * 2 > intervals[i]){
			less_than_mean_freq_count += 2;
		}
	}
	// printf("Mean freq greater than %.2f of the quantization_intervals\n", less_than_mean_freq_count * 1.0 / accIntervals);
	free(sum_intervals);
	free(freq_intervals);
	free(intervals_2D);
	free(intervals);

	return powerOf2;
}

unsigned int optimize_intervals_and_compute_mean_intervals_float_3D(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, float * dense_pos, unsigned int * out_mean_count, float ** means){	
	// compute mean
	float mean = 0.0;
	size_t len = r1 * r2 * r3;
	size_t mean_distance = (int) (sqrt(len));
	// printf("mean_distance %ld, mean_sample_size %ld\n", mean_distance, mean_sample_size);
	float * data_pos = oriData;
	size_t offset_count = 0;
	size_t offset_count_2 = 0;
	size_t mean_count = 0;
	while(data_pos - oriData < len){
		mean += *data_pos;
		mean_count ++;
		data_pos += mean_distance;
		offset_count += mean_distance;
		offset_count_2 += mean_distance;
		if(offset_count >= r3){
			offset_count = 0;
			data_pos -= 1;
		}
		if(offset_count_2 >= r2 * r3){
			offset_count_2 = 0;
			data_pos -= 1;
		}
	}
	if(mean_count > 0) mean /= mean_count;
	size_t range = 8192;
	size_t radius = 4096;
	size_t * freq_intervals = (size_t *) malloc(range*sizeof(size_t));
	memset(freq_intervals, 0, range*sizeof(size_t));
	float * sum_intervals = (float *) malloc(range*sizeof(float));
	memset(sum_intervals, 0, range*sizeof(float));
	float value_sum = 0;

	size_t i,j,k, index;
	size_t radiusIndex;
	size_t r23=r2*r3;
	float pred_value = 0, pred_err;
	size_t *intervals = (size_t*)malloc(maxRangeRadius*sizeof(size_t));
	memset(intervals, 0, maxRangeRadius*sizeof(size_t));
	size_t *intervals_2D = (size_t*)malloc(maxRangeRadius*sizeof(size_t));
	memset(intervals_2D, 0, maxRangeRadius*sizeof(size_t));

	float mean_diff;
	ptrdiff_t freq_index;
	size_t totalSampleSize = 0;
	for(i=1;i<r1;i++)
	{
		for(j=1;j<r2;j++)
		{
			for(k=1;k<r3;k++)
			{			
				if((i+j+k)%sampleDistance==0)
				{
					index = i*r23+j*r3+k;
					mean_diff = oriData[index] - mean;
					if(mean_diff > 0) freq_index = (ptrdiff_t)(mean_diff/realPrecision) + radius;
					else freq_index = (ptrdiff_t)(mean_diff/realPrecision) - 1 + radius;
					if(freq_index <= 0){
						freq_intervals[0] ++;
					}
					else if(freq_index >= range){
						freq_intervals[range - 1] ++;
					}
					else{
						freq_intervals[freq_index] ++;
						sum_intervals[freq_index] += oriData[index];
						value_sum += oriData[index];
					}
					totalSampleSize ++;
				}
			}
		}
	}

	// collect frequency
	size_t max_sum = 0;
	size_t max_index = 0;
	size_t tmp_sum;
	size_t * freq_pos = freq_intervals + 1;
	for(size_t i=1; i<range-2; i++){
		tmp_sum = freq_pos[0] + freq_pos[1];
		if(tmp_sum > max_sum){
			max_sum = tmp_sum;
			max_index = i;
		}
		freq_pos ++;
	}
	*dense_pos = mean + realPrecision * (ptrdiff_t)(max_index + 1 - radius);
	printf("Max frequency: %.6f dense_pos: %.6f\n", max_sum * 1.0 / totalSampleSize, *dense_pos);

	// compute intervals
	float dense_point = *dense_pos;
	size_t unflushed_count = 0;

	// debug
	int * typeArray = (int *) malloc(totalSampleSize*sizeof(int));
	size_t typeArray_count = 0;
	int * typeArray_flag = (int *) malloc(totalSampleSize*sizeof(int));
	memset(typeArray, 0, maxRangeRadius*sizeof(int));

	for(i=1;i<r1;i++)
	{
		for(j=1;j<r2;j++)
		{
			for(k=1;k<r3;k++)
			{			
				if((i+j+k)%sampleDistance==0)
				{
					// 3D prediction
					index = i*r23+j*r3+k;
					if(fabs(oriData[index] - dense_point) <= realPrecision){
						typeArray[typeArray_count ++] = 513;
						continue;
					}
					pred_value = oriData[index-1] + oriData[index-r3] + oriData[index-r23] 
					- oriData[index-1-r23] - oriData[index-r3-1] - oriData[index-r3-r23] + oriData[index-r3-r23-1];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (pred_err/realPrecision+1)/2;
					if(radiusIndex>=maxRangeRadius){
						radiusIndex = maxRangeRadius - 1;
						typeArray[typeArray_count ++] = 514;
					}
					else{
						if(typeArray_flag[radiusIndex] > 0){
							typeArray[typeArray_count ++] = 2 * radiusIndex;
							typeArray_flag[radiusIndex] = 0;
						}
						else{
							typeArray[typeArray_count ++] = 2 * radiusIndex + 1;
							typeArray_flag[radiusIndex] = 1;
						}
					}
					intervals[radiusIndex]++;

					// 2D prediction
					pred_value = oriData[index-1] + oriData[index-r3] - oriData[index-1-r3];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (pred_err/realPrecision+1)/2;
					if(radiusIndex>=maxRangeRadius){
						radiusIndex = maxRangeRadius - 1;
					}
					intervals_2D[radiusIndex]++;

					// 2D another direction
					pred_value = oriData[index-1] + oriData[index-r23] - oriData[index-r23-1];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (pred_err/realPrecision+1)/2;
					if(radiusIndex>=maxRangeRadius){
						radiusIndex = maxRangeRadius - 1;
					}
					intervals_2D[radiusIndex]++;

					// 2D 3rd direction
					pred_value = oriData[index-r3] + oriData[index-r23] - oriData[index-r23-r3];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (pred_err/realPrecision+1)/2;
					if(radiusIndex>=maxRangeRadius){
						radiusIndex = maxRangeRadius - 1;
					}
					intervals_2D[radiusIndex]++;

					unflushed_count ++;
				}
			}
		}
	}
	printf("unflushed_count: %ld with respect to total sample count %ld\n", unflushed_count, totalSampleSize);
	//compute the appropriate number
	size_t targetCount = unflushed_count*predThreshold;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	printf("Estimated powerOf2: %ld, sum %ld percent %.4f\n", powerOf2, sum, sum * 1.0 / unflushed_count);
	// use block or not
	// qsort(intervals, maxRangeRadius, sizeof(size_t), cmp);
	// qsort(intervals_2D, maxRangeRadius, sizeof(size_t), cmp);

	unsigned int block_opt[5] = {4, 8, 12, 16, 20};
	float block_3D_est_prob[5];
	for(i=0; i<5; i++){
		float tmp = 1 - 1.0 / block_opt[i]; 
		block_3D_est_prob[i] = tmp * tmp * tmp;
	}
	float block_br[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	float inf_br = 0.0;
	float p2D_br = 0.0;
	double block_prob;
	double block_accum[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	double probs_3D, probs_2D;
	size_t accum_2D = 0;
	size_t accum_3D = 0;
	for(i=0; i<powerOf2/2; i++){
		probs_3D = intervals[i] * 1.0 / 2 / totalSampleSize;
		// printf("%d: %.6f\n", i, probs_3D);
		accum_3D += intervals[i];
		// INTERVALS_2D already contains number of 3 directions
		probs_2D = intervals_2D[i] * 1.0 / 2 / (3 * totalSampleSize);
		accum_2D += intervals_2D[i];
		for(int j=0; j<5; j++){
			block_prob = (block_3D_est_prob[j] * probs_3D + (1 - block_3D_est_prob[j]) * probs_2D);
			if(block_prob > 0){
				block_br[j] += - block_prob * log2(block_prob);
				block_accum[j] += block_prob;
			}
		}
		if(probs_3D > 0){
			inf_br += - probs_3D * log2(probs_3D);
			accum_3D += probs_3D;
		}
		if(probs_2D > 0){
			p2D_br += - probs_2D * log2(probs_2D);
			accum_2D += probs_2D;
		}
	}
	printf("accum count: %ld %ld\n", accum_3D, accum_2D);
	{
		//huffman encoding result
		SZ_Reset(allNodes, stateNum);
		unsigned char * tmp = (unsigned char *)malloc(totalSampleSize*sizeof(int));
		size_t enCodeSize = 0;
		encode_withTree(typeArray, totalSampleSize, &tmp, &enCodeSize);
		printf("Huffman bitrate: %.4f\n", enCodeSize * 32.0 / (totalSampleSize * sizeof(int)));
		free(tmp);
		// int status;
		// writeIntData_inBytes(typeArray, totalSampleSize, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/type_array.dat", &status);
		// writeIntData_inBytes(intervals, maxRangeRadius, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/intervals.dat", &status);
		free(typeArray);
	}

	double flushed_prob = 1 - unflushed_count * 1.0 / totalSampleSize; 
	float flushed_br = - flushed_prob * log2(flushed_prob);
	printf("Flushed prob: %.4f\n", flushed_prob);
	double unpred_prob;
	float unpred_br;
	double unpred_3D = (unflushed_count - accum_3D) * 1.0 / totalSampleSize;
	double unpred_2D = (unflushed_count * 3 - accum_2D) * 1.0 / (3 * totalSampleSize);
	printf("Unpredicted prob: ");
	for(int j=0; j<5; j++){
		unpred_prob = (block_3D_est_prob[j] * unpred_3D + (1 - block_3D_est_prob[j]) * unpred_2D);
		if(unpred_prob > 0) unpred_br = unpred_prob * log2(unpred_prob);
		else{
			unpred_prob = 0;
			unpred_br = 0;
		}
		printf("%.6f ", unpred_prob);		
		//				2 p log p  ---------------		unpred 									block overhead							
		block_br[j] = 2 * block_br[j] + flushed_br + unpred_br;// + unpred_prob * 32;// + 40.0 / (block_opt[j] * block_opt[j] * block_opt[j]);
	}
	printf("%.6f %.6f\n", unpred_3D, unpred_2D);
	unpred_prob = unpred_3D;
	unpred_br = unpred_prob * log2(unpred_prob);
	if(unpred_prob > 0) unpred_br = unpred_prob * log2(unpred_prob);
	else{
		unpred_prob = 0;
		unpred_br = 0;
	}
	inf_br = 2 * inf_br + flushed_br + unpred_br;// + unpred_prob * 32;
	unpred_prob = unpred_2D;
	unpred_br = unpred_prob * log2(unpred_prob);
	if(unpred_prob > 0) unpred_br = unpred_prob * log2(unpred_prob);
	else{
		unpred_prob = 0;
		unpred_br = 0;
	}
	p2D_br = 2 * p2D_br + flushed_br + unpred_br;// + unpred_prob * 32;

	printf("Blocked bit cost: ");
	for(int j=0; j<5; j++){
		printf("%.6f ", block_br[j]);
	}
	printf("\nNonblocked bit cost: 3D %.6f\n2D %.6f\n", inf_br, p2D_br);

	printf("Overhead over 3D:\n");
	for(int j=0; j<5; j++){
		printf("%.4f%% ", (block_br[j] / inf_br - 1)*100);
	}
	printf("\n");
	// TODO: estimate mse and PSNR???
	// exit(0);

	// float block_cost[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
	// for(int j=0; j<5; j++){
	// 	block_cost[j] = block_br[j] + 0;
	// }
	// printf("Overall cost: ");
	// for(int j=0; j<5; j++){
	// 	printf("%.6f ", block_cost[j]);
	// }
	// printf("%.6f %.6f\n", inf_br, p2D_br);

	// limit max_mean_count
	size_t max_mean_count = powerOf2 - accIntervals;
	if(max_mean_count == 0) max_mean_count = 1;
	// do not limit max_mean_count
	// size_t max_mean_count = accIntervals;
	// estimate on how many mean to use
	mean_count = 1;
	qsort(intervals, maxRangeRadius, sizeof(size_t), cmp);
	size_t threshold = intervals[(int)(accIntervals/4)];
	if(threshold < 0.05 * totalSampleSize) threshold = 0.05 * totalSampleSize;
	printf("threshold: %ld\n", threshold);
	size_t * left_freq_pos = freq_intervals + max_index - 2;
	size_t * right_freq_pos = freq_intervals + max_index + 2;
	size_t tmp_sum_left, tmp_sum_right;
	while(mean_count < max_mean_count){
		tmp_sum_left = left_freq_pos[0] + left_freq_pos[1];
		tmp_sum_right = right_freq_pos[0] + right_freq_pos[1];
		if(tmp_sum_left > tmp_sum_right){
			if(tmp_sum_left > threshold){
				mean_count ++;
				left_freq_pos -= 2;
			}
			else break;
		}
		else{
			if(tmp_sum_right > threshold){
				mean_count ++;
				right_freq_pos += 2;
			}
			else break;
		}
		if(left_freq_pos <= freq_intervals + 1) break;
		if(right_freq_pos >= freq_intervals + range - 2) break;
	}
	// split 1 mean to 2
	if(mean_count == 1 && max_sum > 0.9 * totalSampleSize){
		printf("Split one mean to 2\n");
		mean_count = 2;
		left_freq_pos = left_freq_pos - 1;
	}
	*means = (float *) malloc(mean_count * sizeof(float));
	freq_pos = left_freq_pos + 2;
	float * sum_pos = sum_intervals + (freq_pos - freq_intervals);
	for(unsigned int i=0; i<mean_count; i++){
		(*means)[i] = (sum_pos[0] + sum_pos[1]) / (freq_pos[0] + freq_pos[1]);
		if((*means)[i] > 10000){
			printf("mean: %.6f\n", (*means)[i]);
		}
		sum_pos += 2, freq_pos += 2;
	}
	*out_mean_count = (unsigned int)mean_count;
	for(int i=0; i<mean_count; i++){
		printf("%.6f ", (*means)[i]);
	}
	printf("\nQI count: %d Mean count: %d\n", accIntervals, mean_count);
	
	free(sum_intervals);
	free(freq_intervals);
	free(intervals_2D);
	free(intervals);

	// do not limit max_mean_count
	// powerOf2 = roundUpToPowerOf2(accIntervals + mean_count);
	return powerOf2;
}

unsigned int optimize_intervals_float_4D(float *oriData, size_t r1, size_t r2, size_t r3, size_t r4, double realPrecision)
{
	size_t i,j,k,l, index;
	size_t radiusIndex;
	size_t r234=r2*r3*r4;
	size_t r34=r3*r4;
	float pred_value = 0, pred_err;
	size_t *intervals = (size_t*)malloc(maxRangeRadius*sizeof(size_t));
	memset(intervals, 0, maxRangeRadius*sizeof(size_t));
	size_t totalSampleSize = (r1-1)*(r2-1)*(r3-1)*(r4-1)/sampleDistance;
	for(i=1;i<r1;i++)
	{
		for(j=1;j<r2;j++)
		{
			for(k=1;k<r3;k++)
			{
				for (l=1;l<r4;l++)
				{
					if((i+j+k+l)%sampleDistance==0)
					{
						index = i*r234+j*r34+k*r4+l;
						pred_value = oriData[index-1] + oriData[index-r3] + oriData[index-r34]
								- oriData[index-1-r34] - oriData[index-r4-1] - oriData[index-r4-r34] + oriData[index-r4-r34-1];
						pred_err = fabs(pred_value - oriData[index]);
						radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
						if(radiusIndex>=maxRangeRadius)
							radiusIndex = maxRangeRadius - 1;
						intervals[radiusIndex]++;
					}
				}
			}
		}
	}
	//compute the appropriate number
	size_t targetCount = totalSampleSize*predThreshold;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;

	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	return powerOf2;
}

TightDataPointStorageF* SZ_compress_float_1D_MDQ(float *oriData, 
size_t dataLength, double realPrecision, float valueRangeSize, float medianValue_f)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
		quantization_intervals = optimize_intervals_float_1D(oriData, dataLength, realPrecision);
	else
		quantization_intervals = intvCapacity;
	updateQuantizationInfo(quantization_intervals);	
	//clearHuffmanMem();
	size_t i;
	int reqLength;
	float medianValue = medianValue_f;
	short reqExpo = getPrecisionReqLength_float((float)realPrecision);
	short radExpo = getExponent_float(valueRangeSize/2);
	
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);	

	int* type = (int*) malloc(dataLength*sizeof(int));
		
	float* spaceFillingValue = oriData; //
	
	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);
	
	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);
	
	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);
	
	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);
	
	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);
	
	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;
	float last3CmprsData[3] = {0};

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
				
	//add the first data	
	type[0] = 0;
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_float(last3CmprsData, vce->data);
	//printf("%.30G\n",last3CmprsData[0]);	
		
	//add the second data
	type[1] = 0;
	compressSingleFloatValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_float(last3CmprsData, vce->data);
	//printf("%.30G\n",last3CmprsData[0]);	
	
	int state;
	float lcf, qcf;	
	double checkRadius;
	float curData;
	float pred;
	float predAbsErr;
	float min_pred, minErr, minIndex;
	int a = 0;		
	checkRadius = (intvCapacity-1)*realPrecision;
	double interval = 2*realPrecision;
	
	for(i=2;i<dataLength;i++)
	{
//		if(i==2869438)
//			printf("i=%d\n", i);
		curData = spaceFillingValue[i];
		//pred = 2*last3CmprsData[0] - last3CmprsData[1];
		pred = last3CmprsData[0];
		predAbsErr = fabs(curData - pred);	
		if(predAbsErr<=checkRadius)
		{
			state = (predAbsErr/realPrecision+1)/2;
			if(curData>=pred)
			{
				type[i] = intvRadius+state;
				pred = pred + state*interval;
			}
			else //curData<pred
			{
				type[i] = intvRadius-state;
				pred = pred - state*interval;
			}
/*			if(type[i]==0)
				printf("err:type[%d]=0\n", i);*/
				
			//double-check the prediction error in case of machine-epsilon impact	
			if(fabs(curData-pred)>realPrecision)
			{	
				type[i] = 0;				
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);		
				
				listAdd_float(last3CmprsData, vce->data);		
			}
			else
				listAdd_float(last3CmprsData, pred);
				
			continue;
		}
		
		//unpredictable data processing		
		type[i] = 0;		
		compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);

		listAdd_float(last3CmprsData, vce->data);	
	}//end of for
		
//	char* expSegmentsInBytes;
//	int expSegmentsInBytes_size = convertESCToBytes(esc, &expSegmentsInBytes);
	size_t exactDataNum = exactLeadNumArray->size;
	
	TightDataPointStorageF* tdps;
			
	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum, 
			type, exactMidByteArray->array, exactMidByteArray->size,  
			exactLeadNumArray->array,  
			resiBitArray->array, resiBitArray->size, 
			resiBitsLength,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

//sdi:Debug
/*	int sum =0;
	for(i=0;i<dataLength;i++)
		if(type[i]==0) sum++;
	printf("opt_quantizations=%d, exactDataNum=%d, sum=%d\n",quantization_intervals, exactDataNum, sum);*/

//	writeUShortData(type, dataLength, "compressStateBytes.sb");
//	unsigned short type_[dataLength];
//	SZ_Reset();
//	decode_withTree(tdps->typeArray, tdps->typeArray_size, type_);	
//	printf("tdps->typeArray_size=%d\n", tdps->typeArray_size);
	
//	printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d\n", 
//			exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size);
	
//	for(i = 3800;i<3844;i++)
//		printf("exactLeadNumArray->array[%d]=%d\n",i,exactLeadNumArray->array[i]);
	
	//free memory
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);	
	free(vce);
	free(lce);	
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);
	
	return tdps;
}

void SZ_compress_args_float_StoreOriData(float* oriData, size_t dataLength, TightDataPointStorageF* tdps, 
unsigned char** newByteData, size_t *outSize)
{
	int floatSize=sizeof(float);	
	size_t k = 0, i;
	tdps->isLossless = 1;
	size_t totalByteLength = 3 + MetaDataByteLength + SZ_SIZE_TYPE + 1 + floatSize*dataLength;
	*newByteData = (unsigned char*)malloc(totalByteLength);
	
	unsigned char dsLengthBytes[8];
	for (i = 0; i < 3; i++)//3
		(*newByteData)[k++] = versionNumber[i];

	if(SZ_SIZE_TYPE==4)//1
		(*newByteData)[k++] = 16; //00010000
	else
		(*newByteData)[k++] = 80;	//01010000: 01000000 indicates the SZ_SIZE_TYPE=8
	
	convertSZParamsToBytes(conf_params, &((*newByteData)[k]));
	k = k + MetaDataByteLength;	
	
	sizeToBytes(dsLengthBytes,dataLength); //SZ_SIZE_TYPE: 4 or 8	
	for (i = 0; i < SZ_SIZE_TYPE; i++)
		(*newByteData)[k++] = dsLengthBytes[i];
		
	if(sysEndianType==BIG_ENDIAN_SYSTEM)
		memcpy((*newByteData)+4+MetaDataByteLength+SZ_SIZE_TYPE, oriData, dataLength*floatSize);
	else
	{
		unsigned char* p = (*newByteData)+4+MetaDataByteLength+SZ_SIZE_TYPE;
		for(i=0;i<dataLength;i++,p+=floatSize)
			floatToBytes(p, oriData[i]);
	}	
	*outSize = totalByteLength;
}

void SZ_compress_args_float_NoCkRngeNoGzip_1D(unsigned char** newByteData, float *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f)
{
	SZ_Reset(allNodes, stateNum);	
	TightDataPointStorageF* tdps = SZ_compress_float_1D_MDQ(oriData, dataLength, realPrecision, valueRangeSize, medianValue_f);
			
	//TODO: return bytes....
	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);
	
	if(*outSize>dataLength*sizeof(float))
		SZ_compress_args_float_StoreOriData(oriData, dataLength+2, tdps, newByteData, outSize);
	
	free_TightDataPointStorageF(tdps);
}

TightDataPointStorageF* SZ_compress_float_2D_MDQ(float *oriData, size_t r1, size_t r2, double realPrecision, float valueRangeSize, float medianValue_f)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_2D(oriData, r1, r2, realPrecision);
		updateQuantizationInfo(quantization_intervals);
	}	
	else
		quantization_intervals = intvCapacity;
	size_t i,j; 
	int reqLength;
	float pred1D, pred2D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;
		
	size_t dataLength = r1*r2;	
	
	P0 = (float*)malloc(r2*sizeof(float));
	memset(P0, 0, r2*sizeof(float));
	P1 = (float*)malloc(r2*sizeof(float));
	memset(P1, 0, r2*sizeof(float));
		
	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);	

	int* type = (int*) malloc(dataLength*sizeof(int));
	//type[dataLength]=0;
		
	float* spaceFillingValue = oriData; //

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);
	
	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);
	
	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);
	
	type[0] = 0;
	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);
	
	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
			
	/* Process Row-0 data 0*/
	type[0] = 0;
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	float curData;

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	curData = spaceFillingValue[1];
	diff = curData - pred1D;

	itvNum =  fabs(diff)/realPrecision + 1;

	if (itvNum < intvCapacity)
	{
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;

		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(spaceFillingValue[1]-P1[1])>realPrecision)
		{	
			type[1] = 0;			
			compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
			
			P1[1] = vce->data;	
		}
	}
	else
	{
		type[1] = 0;
		compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-0 data 2 --> data r2-1 */
	for (j = 2; j < r2; j++)
	{
		pred1D = 2*P1[j-1] - P1[j-2];
		curData = spaceFillingValue[j];
		diff = curData - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
		
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[j])>realPrecision)
			{	
				type[j] = 0;				
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
				
				P1[j] = vce->data;	
			}
		}
		else
		{
			type[j] = 0;
			compressSingleFloatValue(vce,curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[j] = vce->data;
		}
	}

	/* Process Row-1 --> Row-r1-1 */
	size_t index;
	for (i = 1; i < r1; i++)
	{	
		/* Process row-i data 0 */
		index = i*r2;
		pred1D = P1[0];
		curData = spaceFillingValue[index];
		diff = curData - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;

			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P0[0])>realPrecision)
			{	
				type[index] = 0;				
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
				
				P0[0] = vce->data;	
			}
		}
		else
		{
			type[index] = 0;
			compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}
									
		/* Process row-i data 1 --> r2-1*/
		for (j = 1; j < r2; j++)
		{
			index = i*r2+j;
			pred2D = P0[j-1] + P1[j] - P1[j-1];

			curData = spaceFillingValue[index];
			diff = curData - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
			
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[j])>realPrecision)
				{	
					type[index] = 0;					
					compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
					
					P0[j] = vce->data;	
				}			
			}
			else
			{
				type[index] = 0;
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}
	
	if(r2!=1)
		free(P0);
	free(P1);			
	size_t exactDataNum = exactLeadNumArray->size;
	
	TightDataPointStorageF* tdps;
			
	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum, 
			type, exactMidByteArray->array, exactMidByteArray->size,  
			exactLeadNumArray->array,  
			resiBitArray->array, resiBitArray->size, 
			resiBitsLength, 
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

//	printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d\n", 
//			exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size);
	
//	for(i = 3800;i<3844;i++)
//		printf("exactLeadNumArray->array[%d]=%d\n",i,exactLeadNumArray->array[i]);
	
	//free memory
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);
	
	return tdps;	
}

/**
 * 
 * Note: @r1 is high dimension
 * 		 @r2 is low dimension 
 * */
void SZ_compress_args_float_NoCkRngeNoGzip_2D(unsigned char** newByteData, float *oriData, size_t r1, size_t r2, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f)
{
	SZ_Reset(allNodes, stateNum);	
		
	TightDataPointStorageF* tdps = SZ_compress_float_2D_MDQ(oriData, r1, r2, realPrecision, valueRangeSize, medianValue_f);

	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

	size_t dataLength = r1*r2;
	if(*outSize>dataLength*sizeof(float))
		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);
	
	free_TightDataPointStorageF(tdps);	
}

TightDataPointStorageF* SZ_compress_float_3D_MDQ(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, float valueRangeSize, float medianValue_f)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_3D(oriData, r1, r2, r3, realPrecision);
		updateQuantizationInfo(quantization_intervals);
	}	
	else
		quantization_intervals = intvCapacity;
	//clearHuffmanMem();
	size_t i,j,k; 
	int reqLength;
	float pred1D, pred2D, pred3D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	size_t dataLength = r1*r2*r3;
	size_t r23 = r2*r3;
	P0 = (float*)malloc(r23*sizeof(float));
	P1 = (float*)malloc(r23*sizeof(float));

	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);	

	int* type = (int*) malloc(dataLength*sizeof(int));

	float* spaceFillingValue = oriData; //

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);
	
	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));


	///////////////////////////	Process layer-0 ///////////////////////////
	/* Process Row-0 data 0*/
	type[0] = 0;
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	float curData;

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	curData = spaceFillingValue[1];
	diff = curData - pred1D;

	itvNum = fabs(diff)/realPrecision + 1;

	if (itvNum < intvCapacity)
	{
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
		
		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(curData-P1[1])>realPrecision)
		{	
			type[1] = 0;			
			compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
			
			P1[1] = vce->data;	
		}		
	}
	else
	{
		type[1] = 0;
		compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++)
	{
		pred1D = 2*P1[j-1] - P1[j-2];
		curData = spaceFillingValue[j];
		diff = curData - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
			
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[j])>realPrecision)
			{	
				type[j] = 0;				
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
				
				P1[j] = vce->data;	
			}			
		}
		else
		{
			type[j] = 0;
			compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[j] = vce->data;
		}
	}

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = P1[index-r3];
		curData = spaceFillingValue[index];
		diff = curData - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P1[index] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[index])>realPrecision)
			{	
				type[index] = 0;				
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
				
				P1[index] = vce->data;	
			}			
		}
		else
		{
			type[index] = 0;
			compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index] = vce->data;
		}

		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];

			curData = spaceFillingValue[index];
			diff = curData - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[index])>realPrecision)
				{	
					type[index] = 0;					
					compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
					
					P1[index] = vce->data;	
				}				
			}
			else
			{
				type[index] = 0;
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index] = vce->data;
			}
		}
	}


	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = P1[0];
		curData = spaceFillingValue[index];
		diff = curData - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P0[0])>realPrecision)
			{	
				type[index] = 0;				
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
				
				P0[0] = vce->data;	
			}			
		}
		else
		{
			type[index] = 0;
			compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}


	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = P0[j-1] + P1[j] - P1[j-1];
			curData = spaceFillingValue[index];
			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[j])>realPrecision)
				{	
					type[index] = 0;					
					compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
					
					P0[j] = vce->data;	
				}
			}
			else
			{
				type[index] = 0;
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
			curData = spaceFillingValue[index];
			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[index2D])>realPrecision)
				{	
					type[index] = 0;					
					compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
					
					P0[index2D] = vce->data;	
				}				
			}
			else
			{
				type[index] = 0;
				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
//				if(k==63&&i==43&&j==27)
//					printf("i=%d\n", i);
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
				curData = spaceFillingValue[index];
				diff = curData - pred3D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
					
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[index2D])>realPrecision)
					{	
						type[index] = 0;						
						compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
						updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
						memcpy(preDataBytes,vce->curBytes,4);
						addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);	
						
						P0[index2D] = vce->data;	
					}					
				}
				else
				{
					type[index] = 0;
					compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}
		}

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}
	if(r23!=1)
		free(P0);
	free(P1);
	size_t exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitsLength, 
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

//sdi:Debug
/*	int sum =0;
	for(i=0;i<dataLength;i++)
		if(type[i]==0) sum++;
	printf("opt_quantizations=%d, exactDataNum=%d, sum=%d\n",quantization_intervals, exactDataNum, sum);
*/

//	printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d\n",
//			exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size);

	//free memory
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);	
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);
	
	return tdps;	
}


void SZ_compress_args_float_NoCkRngeNoGzip_3D(unsigned char** newByteData, float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f)
{
	SZ_Reset(allNodes, stateNum);

	TightDataPointStorageF* tdps = SZ_compress_float_3D_MDQ(oriData, r1, r2, r3, realPrecision, valueRangeSize, medianValue_f);

	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

	int dataLength = r1*r2*r3;
	if(*outSize>dataLength*sizeof(float))
		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);
}


TightDataPointStorageF* SZ_compress_float_4D_MDQ(float *oriData, size_t r1, size_t r2, size_t r3, size_t r4, double realPrecision, float valueRangeSize, float medianValue_f)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_4D(oriData, r1, r2, r3, r4, realPrecision);
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;

	size_t i,j,k; 
	int reqLength;
	float pred1D, pred2D, pred3D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	size_t dataLength = r1*r2*r3*r4;

	size_t r234 = r2*r3*r4;
	size_t r34 = r3*r4;

	P0 = (float*)malloc(r34*sizeof(float));
	P1 = (float*)malloc(r34*sizeof(float));

	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));

	float* spaceFillingValue = oriData; //

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));


	size_t l;
	for (l = 0; l < r1; l++)
	{

		///////////////////////////	Process layer-0 ///////////////////////////
		/* Process Row-0 data 0*/
		size_t index = l*r234;
		size_t index2D = 0;

		type[index] = 0;
		compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[index2D] = vce->data;

		/* Process Row-0 data 1*/
		index = l*r234+1;
		index2D = 1;

		pred1D = P1[index2D-1];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P1[index2D] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
		}
		else
		{
			type[index] = 0;
			compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index2D] = vce->data;
		}

		/* Process Row-0 data 2 --> data r4-1 */
		for (j = 2; j < r4; j++)
		{
			index = l*r234+j;
			index2D = j;

			pred1D = 2*P1[index2D-1] - P1[index2D-2];
			diff = spaceFillingValue[index] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index2D] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				type[index] = 0;
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index2D] = vce->data;
			}
		}

		/* Process Row-1 --> Row-r3-1 */
		for (i = 1; i < r3; i++)
		{
			/* Process row-i data 0 */
			index = l*r234+i*r4;
			index2D = i*r4;

			pred1D = P1[index2D-r4];
			diff = spaceFillingValue[index] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index2D] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				type[index] = 0;
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index2D] = vce->data;
			}

			/* Process row-i data 1 --> data r4-1*/
			for (j = 1; j < r4; j++)
			{
				index = l*r234+i*r4+j;
				index2D = i*r4+j;

				pred2D = P1[index2D-1] + P1[index2D-r4] - P1[index2D-r4-1];

				diff = spaceFillingValue[index] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P1[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				}
				else
				{
					type[index] = 0;
					compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P1[index2D] = vce->data;
				}
			}
		}


		///////////////////////////	Process layer-1 --> layer-r2-1 ///////////////////////////

		for (k = 1; k < r2; k++)
		{
			/* Process Row-0 data 0*/
			index = l*r234+k*r34;
			index2D = 0;

			pred1D = P1[index2D];
			diff = spaceFillingValue[index] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				type[index] = 0;
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}

			/* Process Row-0 data 1 --> data r4-1 */
			for (j = 1; j < r4; j++)
			{
				index = l*r234+k*r34+j;
				index2D = j;

				pred2D = P0[index2D-1] + P1[index2D] - P1[index2D-1];
				diff = spaceFillingValue[index] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				}
				else
				{
					type[index] = 0;
					compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}

			/* Process Row-1 --> Row-r3-1 */
			for (i = 1; i < r3; i++)
			{
				/* Process Row-i data 0 */
				index = l*r234+k*r34+i*r4;
				index2D = i*r4;

				pred2D = P0[index2D-r4] + P1[index2D] - P1[index2D-r4];
				diff = spaceFillingValue[index] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				}
				else
				{
					type[index] = 0;
					compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}

				/* Process Row-i data 1 --> data r4-1 */
				for (j = 1; j < r4; j++)
				{
					index = l*r234+k*r34+i*r4+j;
					index2D = i*r4+j;

					pred3D = P0[index2D-1] + P0[index2D-r4]+ P1[index2D] - P0[index2D-r4-1] - P1[index2D-r4] - P1[index2D-1] + P1[index2D-r4-1];
					diff = spaceFillingValue[index] - pred3D;


					itvNum = fabs(diff)/realPrecision + 1;

					if (itvNum < intvCapacity)
					{
						if (diff < 0) itvNum = -itvNum;
						type[index] = (int) (itvNum/2) + intvRadius;
						P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
					}
					else
					{
						type[index] = 0;
						compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
						updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
						memcpy(preDataBytes,vce->curBytes,4);
						addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
						P0[index2D] = vce->data;
					}
				}
			}

			float *Pt;
			Pt = P1;
			P1 = P0;
			P0 = Pt;
		}
	}

	free(P0);
	free(P1);
	size_t exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitsLength,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

	//free memory
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);

	return tdps;
}

void SZ_compress_args_float_NoCkRngeNoGzip_4D(unsigned char** newByteData, float *oriData, size_t r1, size_t r2, size_t r3, size_t r4, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f)
{
	SZ_Reset(allNodes, stateNum);

	TightDataPointStorageF* tdps = SZ_compress_float_4D_MDQ(oriData, r1, r2, r3, r4, realPrecision, valueRangeSize, medianValue_f);

	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

	int dataLength = r1*r2*r3*r4;
	if(*outSize>dataLength*sizeof(float))
		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);
}

void SZ_compress_args_float_withinRange(unsigned char** newByteData, float *oriData, size_t dataLength, size_t *outSize)
{
	TightDataPointStorageF* tdps = (TightDataPointStorageF*) malloc(sizeof(TightDataPointStorageF));
	tdps->rtypeArray = NULL;
	tdps->typeArray = NULL;	
	tdps->leadNumArray = NULL;
	tdps->residualMidBits = NULL;
	
	tdps->allSameData = 1;
	tdps->dataSeriesLength = dataLength;
	tdps->exactMidBytes = (unsigned char*)malloc(sizeof(unsigned char)*4);
	tdps->pwrErrBoundBytes = NULL;
	tdps->isLossless = 0;
	float value = oriData[0];
	floatToBytes(tdps->exactMidBytes, value);
	tdps->exactMidBytes_size = 4;
	
	size_t tmpOutSize;
	//unsigned char *tmpByteData;
	convertTDPStoFlatBytes_float(tdps, newByteData, &tmpOutSize);

	//*newByteData = (unsigned char*)malloc(sizeof(unsigned char)*12); //for floating-point data (1+3+4+4)
	//memcpy(*newByteData, tmpByteData, 12);
	*outSize = tmpOutSize; //8+SZ_SIZE_TYPE; //8==3+1+4(float_size)
	free_TightDataPointStorageF(tdps);	
}

int SZ_compress_args_float_wRngeNoGzip(unsigned char** newByteData, float *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio, double pwrErrRatio)
{
	int status = SZ_SCES;
	size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
	float valueRangeSize = 0, medianValue = 0;
	
	float min = computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	float max = min+valueRangeSize;
	double realPrecision = getRealPrecision_float(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);
		
	if(valueRangeSize <= realPrecision)
	{
		SZ_compress_args_float_withinRange(newByteData, oriData, dataLength, outSize);
	}
	else
	{
//		SZ_compress_args_float_NoCkRngeNoGzip_2D(newByteData, oriData, r2, r1, realPrecision, outSize);
		if(r5==0&&r4==0&&r3==0&&r2==0)
		{
			if(errBoundMode>=PW_REL)
			{	
				//SZ_compress_args_float_NoCkRngeNoGzip_1D_pwr(newByteData, oriData, realPrecision, r1, outSize, min, max);
				SZ_compress_args_float_NoCkRngeNoGzip_1D_pwrgroup(newByteData, oriData, r1, absErr_Bound, relBoundRatio, pwrErrRatio, valueRangeSize, medianValue, outSize);
			}
			else
				SZ_compress_args_float_NoCkRngeNoGzip_1D(newByteData, oriData, r1, realPrecision, outSize, valueRangeSize, medianValue);
		}
		else if(r5==0&&r4==0&&r3==0)
		{
			if(errBoundMode>=PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_2D_pwr(newByteData, oriData, realPrecision, r2, r1, outSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_2D(newByteData, oriData, r2, r1, realPrecision, outSize, valueRangeSize, medianValue);
		}
		else if(r5==0&&r4==0)
		{
			if(errBoundMode>=PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_3D_pwr(newByteData, oriData, realPrecision, r3, r2, r1, outSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_3D(newByteData, oriData, r3, r2, r1, realPrecision, outSize, valueRangeSize, medianValue);
		}
		else if(r5==0)
		{
			if(errBoundMode>=PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_3D_pwr(newByteData, oriData, realPrecision, r4*r3, r2, r1, outSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_3D(newByteData, oriData, r4*r3, r2, r1, realPrecision, outSize, valueRangeSize, medianValue);
		}
	}
	return status;
}

int SZ_compress_args_float(unsigned char** newByteData, float *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio, double pwRelBoundRatio, int pwrType)
{
	errorBoundMode = errBoundMode;
	if(errBoundMode==PW_REL)
	{
		pw_relBoundRatio = pwRelBoundRatio;	
		pwr_type = pwrType;
		if(pwrType==SZ_PWR_AVG_TYPE && r3 != 0 )
		{
			printf("Error: Current version doesn't support 3D data compression with point-wise relative error bound being based on pwrType=AVG\n");
			exit(0);
			return SZ_NSCS;
		}
	}		
	int status = SZ_SCES;
	size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
	float valueRangeSize = 0, medianValue = 0;
	
	double realPrecision = absErrBound; 
	float min = computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	float max = min+valueRangeSize;
	if(errorBoundMode==PSNR)
	{
		errorBoundMode = conf_params->errorBoundMode = ABS;
		realPrecision = conf_params->absErrBound = computeABSErrBoundFromPSNR(psnr, (double)predThreshold, (double)valueRangeSize);
		//printf("realPrecision=%lf\n", realPrecision);
	}
	else
		realPrecision = getRealPrecision_float(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);		
		
	if(valueRangeSize <= realPrecision)
	{
		SZ_compress_args_float_withinRange(newByteData, oriData, dataLength, outSize);
	}
	else
	{
		size_t tmpOutSize = 0;
		unsigned char* tmpByteData;
		if (r2==0)
		{
			if(errorBoundMode>=PW_REL)
			{
				//SZ_compress_args_float_NoCkRngeNoGzip_1D_pwr(&tmpByteData, oriData, realPrecision, r1, &tmpOutSize, min, max);
				SZ_compress_args_float_NoCkRngeNoGzip_1D_pwrgroup(&tmpByteData, oriData, r1, absErr_Bound, relBoundRatio, pwRelBoundRatio, 
				valueRangeSize, medianValue, &tmpOutSize);
			}
			else
			{
				if(szRandomAccess == SZ_NO_RANDOM_ACCESS)
					SZ_compress_args_float_NoCkRngeNoGzip_1D(&tmpByteData, oriData, r1, realPrecision, &tmpOutSize, valueRangeSize, medianValue);
				else //if(szRandomAccess == SZ_YES_RANDOM_ACCESS)
					tmpByteData = SZ_compress_float_1D_MDQ_RA(oriData, r1, realPrecision, &tmpOutSize); 
			}
		}
		else
		if (r3==0) //2d data
		{
			if(errorBoundMode>=PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_2D_pwr(&tmpByteData, oriData, realPrecision, r2, r1, &tmpOutSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_2D(&tmpByteData, oriData, r2, r1, realPrecision, &tmpOutSize, valueRangeSize, medianValue);
		}
		else
		if (r4==0) //3d data
		{
			if(errorBoundMode>=PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_3D_pwr(&tmpByteData, oriData, realPrecision, r3, r2, r1, &tmpOutSize, min, max);
			else
			{
				if(szRandomAccess == SZ_NO_RANDOM_ACCESS)
					SZ_compress_args_float_NoCkRngeNoGzip_3D(&tmpByteData, oriData, r3, r2, r1, realPrecision, &tmpOutSize, valueRangeSize, medianValue);
				else 
					//tmpByteData = SZ_compress_float_3D_MDQ_RA(oriData, r3, r2, r1, realPrecision, &tmpOutSize);
					//tmpByteData = SZ_compress_float_3D_MDQ_nonblocked(oriData, r3, r2, r1, realPrecision, &tmpOutSize);
					tmpByteData = SZ_compress_float_3D_MDQ_nonblocked_with_blocked_regression(oriData, r3, r2, r1, realPrecision, &tmpOutSize);
			}
		}
		else
		if (r5==0)
		{
			if(errorBoundMode>=PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_3D_pwr(&tmpByteData, oriData, realPrecision, r4*r3, r2, r1, &tmpOutSize, min, max);
				//ToDO
				//SZ_compress_args_float_NoCkRngeNoGzip_4D_pwr(&tmpByteData, oriData, r4, r3, r2, r1, &tmpOutSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_4D(&tmpByteData, oriData, r4, r3, r2, r1, realPrecision, &tmpOutSize, valueRangeSize, medianValue);
		}
		else
		{
			printf("Error: doesn't support 5 dimensions for now.\n");
			status = SZ_DERR; //dimension error
		}
		//Call Gzip to do the further compression.
		if(szMode==SZ_BEST_SPEED)
		{
			*outSize = tmpOutSize;
			*newByteData = tmpByteData;
		}
		else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
		{
			int i = 0;
			*outSize = zlib_compress5(tmpByteData, tmpOutSize, newByteData, gzipMode);
			free(tmpByteData);
		}
		else
		{
			printf("Error: Wrong setting of szMode in the float compression.\n");
			status = SZ_MERR; //mode error			
		}
		
		SZ_ReleaseHuffman();
	}
	return status;
}


void computeReqLength_float(double realPrecision, short radExpo, int* reqLength, float* medianValue)
{
	short reqExpo = getPrecisionReqLength_double(realPrecision);
	*reqLength = 9+radExpo - reqExpo; //radExpo-reqExpo == reqMantiLength
	if(*reqLength<9)
		*reqLength = 9;
	if(*reqLength>32)
	{	
		*reqLength = 32;
		*medianValue = 0;
	}			
}

//TODO
int SZ_compress_args_float_subblock(unsigned char* compressedBytes, float *oriData,
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
size_t s5, size_t s4, size_t s3, size_t s2, size_t s1,
size_t e5, size_t e4, size_t e3, size_t e2, size_t e1,
size_t *outSize, int errBoundMode, double absErr_Bound, double relBoundRatio)
{
	int status = SZ_SCES;
	float valueRangeSize = 0, medianValue = 0;
	float min = computeRangeSize_float_subblock(oriData, &valueRangeSize, &medianValue, r5, r4, r3, r2, r1, s5, s4, s3, s2, s1, e5, e4, e3, e2, e1);
	float max = min+valueRangeSize;

	double realPrecision = getRealPrecision_float(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);

	if(valueRangeSize <= realPrecision)
	{
		//TODO
		//SZ_compress_args_float_withinRange_subblock();
	}
	else
	{
		if (r2==0)
		{
			if(errBoundMode>=PW_REL)
			{
				//TODO
				//SZ_compress_args_float_NoCkRngeNoGzip_1D_pwr_subblock();
				printf ("Current subblock version does not support point-wise relative error bound.\n");
			}
			else
				SZ_compress_args_float_NoCkRnge_1D_subblock(compressedBytes, oriData, realPrecision, outSize, valueRangeSize, medianValue, r1, s1, e1);
		}
		else
		if (r3==0)
		{
			//TODO
			if(errBoundMode>=PW_REL)
			{
				//TODO
				//SZ_compress_args_float_NoCkRngeNoGzip_2D_pwr_subblock();
				printf ("Current subblock version does not support point-wise relative error bound.\n");
			}
			else
				SZ_compress_args_float_NoCkRnge_2D_subblock(compressedBytes, oriData, realPrecision, outSize, valueRangeSize, medianValue, r2, r1, s2, s1, e2, e1);
		}
		else
		if (r4==0)
		{
			if(errBoundMode>=PW_REL)
			{
				//TODO
				//SZ_compress_args_float_NoCkRngeNoGzip_3D_pwr_subblock();
				printf ("Current subblock version does not support point-wise relative error bound.\n");
			}
			else
				SZ_compress_args_float_NoCkRnge_3D_subblock(compressedBytes, oriData, realPrecision, outSize, valueRangeSize, medianValue, r3, r2, r1, s3, s2, s1, e3, e2, e1);
		}
		else
		if (r5==0)
		{
			if(errBoundMode>=PW_REL)
			{
				//TODO
				//SZ_compress_args_float_NoCkRngeNoGzip_4D_pwr_subblock();
				printf ("Current subblock version does not support point-wise relative error bound.\n");
			}
			else
				SZ_compress_args_float_NoCkRnge_4D_subblock(compressedBytes, oriData, realPrecision, outSize, valueRangeSize, medianValue, r4, r3, r2, r1, s4, s3, s2, s1, e4, e3, e2, e1);
		}
		else
		{
			printf("Error: doesn't support 5 dimensions for now.\n");
			status = SZ_DERR; //dimension error
		}
	}
	SZ_ReleaseHuffman();
	return status;
}

void SZ_compress_args_float_NoCkRnge_1D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f,
size_t r1, size_t s1, size_t e1)
{
	SZ_Reset(allNodes, stateNum);
	TightDataPointStorageF* tdps = SZ_compress_float_1D_MDQ_subblock(oriData, realPrecision, valueRangeSize, medianValue_f, r1, s1, e1);

	if (szMode==SZ_BEST_SPEED)
		convertTDPStoFlatBytes_float_args(tdps, compressedBytes, outSize);
	else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
	{
		unsigned char *tmpCompBytes;
		size_t tmpOutSize;
		convertTDPStoFlatBytes_float(tdps, &tmpCompBytes, &tmpOutSize);
		*outSize = zlib_compress3(tmpCompBytes, tmpOutSize, compressedBytes, gzipMode);
		free(tmpCompBytes);
	}
	else
	{
		printf ("Error: Wrong setting of szMode in the double compression.\n");
	}

	//TODO
//	if(*outSize>dataLength*sizeof(float))
//		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);
}

void SZ_compress_args_float_NoCkRnge_2D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f,
size_t r2, size_t r1, size_t s2, size_t s1, size_t e2, size_t e1)
{
	SZ_Reset(allNodes, stateNum);
	TightDataPointStorageF* tdps = SZ_compress_float_2D_MDQ_subblock(oriData, realPrecision, valueRangeSize, medianValue_f, r2, r1, s2, s1, e2, e1);

	if (szMode==SZ_BEST_SPEED)
		convertTDPStoFlatBytes_float_args(tdps, compressedBytes, outSize);
	else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
	{
		unsigned char *tmpCompBytes;
		size_t tmpOutSize;
		convertTDPStoFlatBytes_float(tdps, &tmpCompBytes, &tmpOutSize);
		*outSize = zlib_compress3(tmpCompBytes, tmpOutSize, compressedBytes, gzipMode);
		free(tmpCompBytes);
	}
	else
	{
		printf ("Error: Wrong setting of szMode in the double compression.\n");
	}

	//TODO
//	if(*outSize>dataLength*sizeof(float))
//		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);
}

void SZ_compress_args_float_NoCkRnge_3D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f,
size_t r3, size_t r2, size_t r1, size_t s3, size_t s2, size_t s1, size_t e3, size_t e2, size_t e1)
{
	SZ_Reset(allNodes, stateNum);
	TightDataPointStorageF* tdps = SZ_compress_float_3D_MDQ_subblock(oriData, realPrecision, valueRangeSize, medianValue_f, r3, r2, r1, s3, s2, s1, e3, e2, e1);

	if (szMode==SZ_BEST_SPEED)
		convertTDPStoFlatBytes_float_args(tdps, compressedBytes, outSize);
	else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
	{
		unsigned char *tmpCompBytes;
		size_t tmpOutSize;
		convertTDPStoFlatBytes_float(tdps, &tmpCompBytes, &tmpOutSize);
		*outSize = zlib_compress3(tmpCompBytes, tmpOutSize, compressedBytes, gzipMode);
		free(tmpCompBytes);
	}
	else
	{
		printf ("Error: Wrong setting of szMode in the double compression.\n");
	}

	//TODO
//	if(*outSize>dataLength*sizeof(float))
//		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);
}

void SZ_compress_args_float_NoCkRnge_4D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f,
size_t r4, size_t r3, size_t r2, size_t r1, size_t s4, size_t s3, size_t s2, size_t s1, size_t e4, size_t e3, size_t e2, size_t e1)
{
	SZ_Reset(allNodes, stateNum);
	TightDataPointStorageF* tdps = SZ_compress_float_4D_MDQ_subblock(oriData, realPrecision, valueRangeSize, medianValue_f, r4, r3, r2, r1, s4, s3, s2, s1, e4, e3, e2, e1);

	if (szMode==SZ_BEST_SPEED)
		convertTDPStoFlatBytes_float_args(tdps, compressedBytes, outSize);
	else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
	{
		unsigned char *tmpCompBytes;
		size_t tmpOutSize;
		convertTDPStoFlatBytes_float(tdps, &tmpCompBytes, &tmpOutSize);
		*outSize = zlib_compress3(tmpCompBytes, tmpOutSize, compressedBytes, gzipMode);
		free(tmpCompBytes);
	}
	else
	{
		printf ("Error: Wrong setting of szMode in the double compression.\n");
	}

	//TODO
//	if(*outSize>dataLength*sizeof(float))
//		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);

}

unsigned int optimize_intervals_float_1D_subblock(float *oriData, double realPrecision, size_t r1, size_t s1, size_t e1)
{
	size_t dataLength = e1 - s1 + 1;
	oriData = oriData + s1;

	size_t i = 0;
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	size_t totalSampleSize = dataLength/sampleDistance;
	for(i=2;i<dataLength;i++)
	{
		if(i%sampleDistance==0)
		{
			pred_value = 2*oriData[i-1] - oriData[i-2];
			//pred_value = oriData[i-1];
			pred_err = fabs(pred_value - oriData[i]);
			radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
			if(radiusIndex>=maxRangeRadius)
				radiusIndex = maxRangeRadius - 1;
			intervals[radiusIndex]++;
		}
	}
	//compute the appropriate number
	size_t targetCount = totalSampleSize*predThreshold;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;

	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	//printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);
	return powerOf2;
}

unsigned int optimize_intervals_float_2D_subblock(float *oriData, double realPrecision, size_t r1, size_t r2, size_t s1, size_t s2, size_t e1, size_t e2)
{
	size_t R1 = e1 - s1 + 1;
	size_t R2 = e2 - s2 + 1;

	size_t i,j, index;
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	size_t totalSampleSize = R1*R2/sampleDistance;
	for(i=s1+1;i<=e1;i++)
	{
		for(j=s2+1;j<=e2;j++)
		{
			if((i+j)%sampleDistance==0)
			{
				index = i*r2+j;
				pred_value = oriData[index-1] + oriData[index-r2] - oriData[index-r2-1];
				pred_err = fabs(pred_value - oriData[index]);
				radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
				if(radiusIndex>=maxRangeRadius)
					radiusIndex = maxRangeRadius - 1;
				intervals[radiusIndex]++;
			}
		}
	}
	//compute the appropriate number
	size_t targetCount = totalSampleSize*predThreshold;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	//printf("maxRangeRadius = %d, accIntervals=%d, powerOf2=%d\n", maxRangeRadius, accIntervals, powerOf2);
	return powerOf2;
}

unsigned int optimize_intervals_float_3D_subblock(float *oriData, double realPrecision, size_t r1, size_t r2, size_t r3, size_t s1, size_t s2, size_t s3, size_t e1, size_t e2, size_t e3)
{
	size_t R1 = e1 - s1 + 1;
	size_t R2 = e2 - s2 + 1;
	size_t R3 = e3 - s3 + 1;

	size_t r23 = r2*r3;

	size_t i,j,k, index;
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	size_t totalSampleSize = R1*R2*R3/sampleDistance;
	for(i=s1+1;i<=e1;i++)
	{
		for(j=s2+1;j<=e2;j++)
		{
			for(k=s3+1;k<=e3;k++)
			{
				if((i+j+k)%sampleDistance==0)
				{
					index = i*r23+j*r3+k;
					pred_value = oriData[index-1] + oriData[index-r3] + oriData[index-r23]
					- oriData[index-1-r23] - oriData[index-r3-1] - oriData[index-r3-r23] + oriData[index-r3-r23-1];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
					if(radiusIndex>=maxRangeRadius)
						radiusIndex = maxRangeRadius - 1;
					intervals[radiusIndex]++;
				}
			}
		}
	}
	//compute the appropriate number
	size_t targetCount = totalSampleSize*predThreshold;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	return powerOf2;
}

unsigned int optimize_intervals_float_4D_subblock(float *oriData, double realPrecision,
size_t r1, size_t r2, size_t r3, size_t r4, size_t s1, size_t s2, size_t s3, size_t s4, size_t e1, size_t e2, size_t e3, size_t e4)
{
	size_t R1 = e1 - s1 + 1;
	size_t R2 = e2 - s2 + 1;
	size_t R3 = e3 - s3 + 1;
	size_t R4 = e4 - s4 + 1;

	size_t r34 = r3*r4;
	size_t r234 = r2*r3*r4;

	size_t i,j,k,l, index;
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	size_t totalSampleSize = R1*R2*R3*R4/sampleDistance;
	for(i=s1+1;i<=e1;i++)
	{
		for(j=s2+1;j<=e2;j++)
		{
			for(k=s3+1;k<=e3;k++)
			{
				for (l=s4+1;l<=e4;l++)
				{
					if((i+j+k+l)%sampleDistance==0)
					{
						index = i*r234+j*r34+k*r4+l;
						pred_value = oriData[index-1] + oriData[index-r4] + oriData[index-r34]
									- oriData[index-1-r34] - oriData[index-r4-1] - oriData[index-r4-r34] + oriData[index-r4-r34-1];
						pred_err = fabs(pred_value - oriData[index]);
						radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
						if(radiusIndex>=maxRangeRadius)
							radiusIndex = maxRangeRadius - 1;
						intervals[radiusIndex]++;
					}
				}
			}
		}
	}
	//compute the appropriate number
	size_t targetCount = totalSampleSize*predThreshold;
	size_t sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;

	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	return powerOf2;
}

TightDataPointStorageF* SZ_compress_float_1D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
size_t r1, size_t s1, size_t e1)
{
	size_t dataLength = e1 - s1 + 1;
	unsigned int quantization_intervals;
	if(optQuantMode==1)
		quantization_intervals = optimize_intervals_float_1D_subblock(oriData, realPrecision, r1, s1, e1);
	else
		quantization_intervals = intvCapacity;
	updateQuantizationInfo(quantization_intervals);

	size_t i; 
	int reqLength;
	float medianValue = medianValue_f;
	short reqExpo = getPrecisionReqLength_float((float)realPrecision);
	short radExpo = getExponent_float(valueRangeSize/2);

	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));

	float* spaceFillingValue = oriData + s1;

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	type[0] = 0;

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;
	float last3CmprsData[3] = {0};

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));

	//add the first data
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_float(last3CmprsData, vce->data);

	//add the second data
	type[1] = 0;
	compressSingleFloatValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_float(last3CmprsData, vce->data);

	int state;
	float lcf, qcf;
	double checkRadius;
	float curData;
	float pred;
	float predAbsErr;
	float min_pred, minErr, minIndex;
	int a = 0;
	checkRadius = (intvCapacity-1)*realPrecision;
	double interval = 2*realPrecision;

	for(i=2;i<dataLength;i++)
	{
		curData = spaceFillingValue[i];
		pred = 2*last3CmprsData[0] - last3CmprsData[1];
		predAbsErr = fabs(curData - pred);
		if(predAbsErr<=checkRadius)
		{
			state = (predAbsErr/realPrecision+1)/2;
			if(curData>=pred)
			{
				type[i] = intvRadius+state;
				pred = pred + state*interval;
			}
			else
			{
				type[i] = intvRadius-state;
				pred = pred - state*interval;
			}

			listAdd_float(last3CmprsData, pred);
			continue;
		}

		//unpredictable data processing
		type[i] = 0;
		compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);

		listAdd_float(last3CmprsData, vce->data);
	}

	size_t exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitsLength,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

	//free memory
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);

	return tdps;
}

TightDataPointStorageF* SZ_compress_float_2D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
size_t r1, size_t r2, size_t s1, size_t s2, size_t e1, size_t e2)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_2D_subblock(oriData, realPrecision, r1, r2, s1, s2, e1, e2);
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;

	size_t i,j; 
	int reqLength;
	float pred1D, pred2D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	size_t R1 = e1 - s1 + 1;
	size_t R2 = e2 - s2 + 1;
	size_t dataLength = R1*R2;

	P0 = (float*)malloc(R2*sizeof(float));
	memset(P0, 0, R2*sizeof(float));
	P1 = (float*)malloc(R2*sizeof(float));
	memset(P1, 0, R2*sizeof(float));

	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));

	float* spaceFillingValue = oriData; //

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));

	/* Process Row-s1 data s2*/
	size_t gIndex;
	size_t lIndex;

	gIndex = s1*r2+s2;
	lIndex = 0;

	type[lIndex] = 0;
	compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	/* Process Row-s1 data s2+1*/
	gIndex = s1*r2+(s2+1);
	lIndex = 1;

	pred1D = P1[0];
	diff = spaceFillingValue[gIndex] - pred1D;

	itvNum =  fabs(diff)/realPrecision + 1;

	if (itvNum < intvCapacity)
	{
		if (diff < 0) itvNum = -itvNum;
		type[lIndex] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
	}
	else
	{
		type[lIndex] = 0;
		compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-s1 data s2+2 --> data e2 */
	for (j = 2; j < R2; j++)
	{
		gIndex = s1*r2+(s2+j);
		lIndex = j;

		pred1D = 2*P1[j-1] - P1[j-2];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[j] = vce->data;
		}
	}

	/* Process Row-s1+1 --> Row-e1 */
	for (i = 1; i < R1; i++)
	{
		/* Process row-s1+i data s2 */
		gIndex = (s1+i)*r2+s2;
		lIndex = i*R2;

		pred1D = P1[0];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}

		/* Process row-s1+i data s2+1 --> e2 */
		for (j = 1; j < R2; j++)
		{
			gIndex = (s1+i)*r2+(s2+j);
			lIndex = i*R2+j;

//			printf ("global index = %d, local index = %d\n", gIndex, lIndex);

			pred2D = P0[j-1] + P1[j] - P1[j-1];

			diff = spaceFillingValue[gIndex] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	free(P0);
	free(P1);
	size_t exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitsLength,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

	//free memory
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);

	return tdps;
}

TightDataPointStorageF* SZ_compress_float_3D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
size_t r1, size_t r2, size_t r3, size_t s1, size_t s2, size_t s3, size_t e1, size_t e2, size_t e3)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_3D_subblock(oriData, realPrecision, r1, r2, r3, s1, s2, s3, e1, e2, e3);
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;

	size_t i,j,k; 
	int reqLength;
	float pred1D, pred2D, pred3D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	size_t R1 = e1 - s1 + 1;
	size_t R2 = e2 - s2 + 1;
	size_t R3 = e3 - s3 + 1;
	size_t dataLength = R1*R2*R3;

	size_t r23 = r2*r3;
	size_t R23 = R2*R3;

	P0 = (float*)malloc(R23*sizeof(float));
	P1 = (float*)malloc(R23*sizeof(float));

	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));
	//type[dataLength]=0;

	float* spaceFillingValue = oriData; //

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));


	///////////////////////////	Process layer-s1 ///////////////////////////
	/* Process Row-s2 data s3*/
	size_t gIndex; 	//global index
	size_t lIndex; 	//local index
	size_t index2D; 	//local 2D index

	gIndex = s1*r23+s2*r3+s3;
	lIndex = 0;
	index2D = 0;

	type[lIndex] = 0;
	compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[index2D] = vce->data;

	/* Process Row-s2 data s3+1*/
	gIndex = s1*r23+s2*r3+s3+1;
	lIndex = 1;
	index2D = 1;

	pred1D = P1[index2D-1];
	diff = spaceFillingValue[gIndex] - pred1D;

	itvNum = fabs(diff)/realPrecision + 1;

	if (itvNum < intvCapacity)
	{
		if (diff < 0) itvNum = -itvNum;
		type[lIndex] = (int) (itvNum/2) + intvRadius;
		P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
	}
	else
	{
		type[lIndex] = 0;
		compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[index2D] = vce->data;
	}

    /* Process Row-s2 data s3+2 --> data e3 */
	for (j = 2; j < R3; j++)
	{
		gIndex = s1*r23+s2*r3+s3+j;
		lIndex = j;
		index2D = j;

		pred1D = 2*P1[index2D-1] - P1[index2D-2];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index2D] = vce->data;
		}
	}

	/* Process Row-s2+1 --> Row-e2 */
	for (i = 1; i < R2; i++)
	{
		/* Process row-s2+i data s3 */
		gIndex = s1*r23+(s2+i)*r3+s3;
		lIndex = i*R3;
		index2D = i*R3;

		pred1D  = P1[index2D-R3];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index2D] = vce->data;
		}

		/* Process row-s2+i data s3+1 --> data e3*/
		for (j = 1; j < R3; j++)
		{
			gIndex = s1*r23+(s2+i)*r3+s3+j;
			lIndex = i*R3+j;
			index2D = i*R3+j;

			pred2D  = P1[index2D-1] + P1[index2D-R3] - P1[index2D-R3-1];
			diff = spaceFillingValue[gIndex] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P1[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index2D] = vce->data;
			}
		}
	}


	///////////////////////////	Process layer-s1+1 --> layer-e1 ///////////////////////////

	for (k = 1; k < R1; k++)
	{
		/* Process Row-s2 data s3*/
		gIndex = (s1+k)*r23+s2*r3+s3;
		lIndex = k*R23;
		index2D = 0;

		pred1D = P1[index2D];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P0[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[index2D] = vce->data;
		}

	    /* Process Row-s2 data s3+1 --> data e3 */
		for (j = 1; j < R3; j++)
		{
			gIndex = (s1+k)*r23+s2*r3+s3+j;
			lIndex = k*R23+j;
			index2D = j;

			pred2D = P0[index2D-1] + P1[index2D] - P1[index2D-1];
			diff = spaceFillingValue[gIndex] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}
		}

	    /* Process Row-s2+1 --> Row-e2 */
		for (i = 1; i < R2; i++)
		{
			/* Process Row-s2+i data s3 */
			gIndex = (s1+k)*r23+(s2+i)*r3+s3;
			lIndex = k*R23+i*R3;
			index2D = i*R3;

			pred2D = P0[index2D-R3] + P1[index2D] - P1[index2D-R3];
			diff = spaceFillingValue[gIndex] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}

			/* Process Row-s2+i data s3+1 --> data e3 */
			for (j = 1; j < R3; j++)
			{
				gIndex = (s1+k)*r23+(s2+i)*r3+s3+j;
				lIndex = k*R23+i*R3+j;
				index2D = i*R3+j;

//				printf ("global index = %d, local index = %d\n", gIndex, lIndex);

				pred3D = P0[index2D-1] + P0[index2D-R3]+ P1[index2D] - P0[index2D-R3-1] - P1[index2D-R3] - P1[index2D-1] + P1[index2D-R3-1];
				diff = spaceFillingValue[gIndex] - pred3D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[lIndex] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred3D + 2 * (type[lIndex] - intvRadius) * realPrecision;
				}
				else
				{
					type[lIndex] = 0;
					compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}
		}

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	free(P0);
	free(P1);
	size_t exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitsLength,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

	//free memory
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);

	return tdps;
}

TightDataPointStorageF* SZ_compress_float_4D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
size_t r1, size_t r2, size_t r3, size_t r4, size_t s1, size_t s2, size_t s3, size_t s4, size_t e1, size_t e2, size_t e3, size_t e4)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_4D_subblock(oriData, realPrecision, r1, r2, r3, r4, s1, s2, s3, s4, e1, e2, e3, e4);
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;

	size_t i,j,k; 
	int reqLength;
	float pred1D, pred2D, pred3D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	size_t R1 = e1 - s1 + 1;
	size_t R2 = e2 - s2 + 1;
	size_t R3 = e3 - s3 + 1;
	size_t R4 = e4 - s4 + 1;

	size_t dataLength = R1*R2*R3*R4;

	size_t r34 = r3*r4;
	size_t r234 = r2*r3*r4;
	size_t R34 = R3*R4;
	size_t R234 = R2*R3*R4;

	P0 = (float*)malloc(R34*sizeof(float));
	P1 = (float*)malloc(R34*sizeof(float));

	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));

	float* spaceFillingValue = oriData; //

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));


	size_t l;
	for (l = 0; l < R1; l++)
	{

		///////////////////////////	Process layer-s2 ///////////////////////////
		/* Process Row-s3 data s4*/
		size_t gIndex; 	//global index
		size_t lIndex; 	//local index
		size_t index2D; 	//local 2D index

		gIndex = (s1+l)*r234+s2*r34+s3*r4+s4;
		lIndex = l*R234;
		index2D = 0;

		type[lIndex] = 0;
		compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[index2D] = vce->data;

		/* Process Row-s3 data s4+1*/
		gIndex = (s1+l)*r234+s2*r34+s3*r4+s4+1;
		lIndex = l*R234+1;
		index2D = 1;

		pred1D = P1[index2D-1];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index2D] = vce->data;
		}

		/* Process Row-s3 data s4+2 --> data e4 */
		for (j = 2; j < R4; j++)
		{
			gIndex = (s1+l)*r234+s2*r34+s3*r4+s4+j;
			lIndex = l*R234+j;
			index2D = j;

			pred1D = 2*P1[index2D-1] - P1[index2D-2];
			diff = spaceFillingValue[gIndex] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index2D] = vce->data;
			}
		}

		/* Process Row-s3+1 --> Row-e3 */
		for (i = 1; i < R3; i++)
		{
			/* Process row-s2+i data s3 */
			gIndex = (s1+l)*r234+s2*r34+(s3+i)*r4+s4;
			lIndex = l*R234+i*R4;
			index2D = i*R4;

			pred1D  = P1[index2D-R4];
			diff = spaceFillingValue[gIndex] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index2D] = vce->data;
			}

			/* Process row-s3+i data s4+1 --> data e4*/
			for (j = 1; j < R4; j++)
			{
				gIndex = (s1+l)*r234+s2*r34+(s3+i)*r4+s4+j;
				lIndex = l*R234+i*R4+j;
				index2D = i*R4+j;

				pred2D  = P1[index2D-1] + P1[index2D-R4] - P1[index2D-R4-1];
				diff = spaceFillingValue[gIndex] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[lIndex] = (int) (itvNum/2) + intvRadius;
					P1[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
				}
				else
				{
					type[lIndex] = 0;
					compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P1[index2D] = vce->data;
				}
			}
		}


		///////////////////////////	Process layer-s2+1 --> layer-e2 ///////////////////////////

		for (k = 1; k < R2; k++)
		{
			/* Process Row-s3 data s4*/
			gIndex = (s1+l)*r234+(s2+k)*r34+s3*r4+s4;
			lIndex = l*R234+k*R34;
			index2D = 0;

			pred1D = P1[index2D];
			diff = spaceFillingValue[gIndex] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}

			/* Process Row-s3 data s4+1 --> data e4 */
			for (j = 1; j < R4; j++)
			{
				gIndex = (s1+l)*r234+(s2+k)*r34+s3*r4+s4+j;
				lIndex = l*R234+k*R34+j;
				index2D = j;

				pred2D = P0[index2D-1] + P1[index2D] - P1[index2D-1];
				diff = spaceFillingValue[gIndex] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[lIndex] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
				}
				else
				{
					type[lIndex] = 0;
					compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}

			/* Process Row-s3+1 --> Row-e3 */
			for (i = 1; i < R3; i++)
			{
				/* Process Row-s3+i data s4 */
				gIndex = (s1+l)*r234+(s2+k)*r34+(s3+i)*r4+s4;
				lIndex = l*R234+k*R34+i*R4;
				index2D = i*R4;

				pred2D = P0[index2D-R4] + P1[index2D] - P1[index2D-R4];
				diff = spaceFillingValue[gIndex] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[lIndex] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
				}
				else
				{
					type[lIndex] = 0;
					compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}

				/* Process Row-s3+i data s4+1 --> data e4 */
				for (j = 1; j < R4; j++)
				{
					gIndex = (s1+l)*r234+(s2+k)*r34+(s3+i)*r4+s4+j;
					lIndex = l*R234+k*R34+i*R4+j;
					index2D = i*R4+j;

//					printf ("global index = %d, local index = %d\n", gIndex, lIndex);

					pred3D = P0[index2D-1] + P0[index2D-R4]+ P1[index2D] - P0[index2D-R4-1] - P1[index2D-R4] - P1[index2D-1] + P1[index2D-R4-1];
					diff = spaceFillingValue[gIndex] - pred3D;

					itvNum = fabs(diff)/realPrecision + 1;

					if (itvNum < intvCapacity)
					{
						if (diff < 0) itvNum = -itvNum;
						type[lIndex] = (int) (itvNum/2) + intvRadius;
						P0[index2D] = pred3D + 2 * (type[lIndex] - intvRadius) * realPrecision;
					}
					else
					{
						type[lIndex] = 0;
						compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
						updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
						memcpy(preDataBytes,vce->curBytes,4);
						addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
						P0[index2D] = vce->data;
					}
				}
			}

			float *Pt;
			Pt = P1;
			P1 = P0;
			P0 = Pt;
		}

	}

	free(P0);
	free(P1);
	size_t exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitsLength,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

	//free memory
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);

	return tdps;
}

size_t SZ_compress_float_1D_MDQ_RA_block_1D_pred(float * block_ori_data, float * mean, float dense_pos, size_t dim_0, size_t block_dim_0, double realPrecision, int * type, DynamicFloatArray * unpredictable_data){


	float sum = 0.0;
	float * data_pos;
	size_t mean_count = 0;
	data_pos = block_ori_data;
	for(size_t i=0; i<block_dim_0; i++){
		if(fabs(*data_pos - dense_pos) <= realPrecision){
			sum += *data_pos;
			mean_count ++;
		}
		data_pos ++;
	}
	if(mean_count > 0) mean[0] = sum / mean_count;
	else mean[0] = 0;
	// printf("SZ_compress_float_1D_MDQ_RA_block mean computation done: %.2f %.2f %d\n", mean[0], sum, mean_count);
	// fflush(stdout);

	unsigned short unpredictable_count = 0;

	float * cur_data_pos = block_ori_data;
	float curData;
	double itvNum;
	double diff;
	float last_over_thres = mean[0];
	float pred1D;
	size_t type_index = 0;
	data_pos = block_ori_data;
	for(size_t i=0; i<block_dim_0; i++){
		curData = *data_pos;
		if(fabs(curData - mean[0]) <= realPrecision){
			type[type_index] = 1;
		}
		else{
			pred1D = last_over_thres;
			diff = curData - pred1D;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity){
				if (diff < 0) itvNum = -itvNum;
				type[type_index] = (int) (itvNum/2) + intvRadius;	
				last_over_thres = pred1D + 2 * (type[type_index] - intvRadius) * realPrecision;
			}
			else{
				type[type_index] = 0;
				// unpredictable_data[unpredictable_count ++] = curData;
				addDFA_Data(unpredictable_data, curData);
				unpredictable_count ++;
				last_over_thres = curData;
			}
		}
		type_index ++;
		data_pos ++;
	}
	return unpredictable_count;

}

size_t SZ_compress_float_2D_MDQ_RA_block_2D_pred(float * block_ori_data, float * mean, float dense_pos, size_t dim_0, size_t dim_1, size_t block_dim_0, size_t block_dim_1, double realPrecision, float * P0, float * P1, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1;

	data_pos = block_ori_data;
	size_t mean_count = 0;
	for(size_t i=0; i<block_dim_0; i++){
		for(size_t j=0; j<block_dim_1; j++){
			if(fabs(*data_pos - dense_pos) <= realPrecision){
				sum += *data_pos;
				mean_count ++;
			}
			data_pos ++;
		}
		data_pos += dim0_offset - block_dim_1;
	}
	if(mean_count > 0) mean[0] = sum / mean_count;
	else mean[0] = 0;

	size_t unpredictable_count = 0;
	size_t r1, r2;
	r1 = block_dim_0;
	r2 = block_dim_1;

	float * cur_data_pos = block_ori_data;
	float curData;
	float pred1D, pred2D;
	double itvNum;
	double diff;
	size_t i, j;
	/* Process Row-0 data 0*/
	curData = cur_data_pos[0];
	if(fabs(curData - mean[0]) <= realPrecision){
		type[0] = 1;
		P1[0] = mean[0];
	}
	else{
		pred1D = mean[0];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type[0] = (int) (itvNum/2) + intvRadius;
			P1[0] = pred1D + 2 * (type[0] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[0])>realPrecision){	
				type[0] = 0;
				P1[0] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}		
		}
		else{
			type[0] = 0;
			P1[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}

	/* Process Row-0 data 1*/
	curData = cur_data_pos[1];
	if(fabs(curData - mean[0]) <= realPrecision){
		type[1] = 1;
		P1[1] = mean[0];
	}
	else{
		pred1D = P1[0];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type[1] = (int) (itvNum/2) + intvRadius;
			P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[1])>realPrecision){	
				type[1] = 0;
				P1[1] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}		
		}
		else{
			type[1] = 0;
			P1[1] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}

    /* Process Row-0 data 2 --> data r2-1 */
	for (j = 2; j < r2; j++)
	{
		curData = cur_data_pos[j];
		if(fabs(curData - mean[0]) <= realPrecision){
			type[j] = 1;
			P1[j] = mean[0];
		}
		else{
			pred1D = 2*P1[j-1] - P1[j-2];
			diff = curData - pred1D;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity){
				if (diff < 0) itvNum = -itvNum;
				type[j] = (int) (itvNum/2) + intvRadius;
				P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[j])>realPrecision){	
					type[j] = 0;
					P1[j] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}			
			}
			else{
				type[j] = 0;
				P1[j] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
	}
	cur_data_pos += dim0_offset;
	/* Process Row-1 --> Row-r1-1 */
	size_t index;
	for (i = 1; i < r1; i++)
	{	
		/* Process row-i data 0 */
		index = i*r2;
		curData = cur_data_pos[0];
		if(fabs(curData - mean[0]) <= realPrecision){
			type[index] = 1;
			P0[0] = mean[0];
		}
		else{
			pred1D = P1[0];
			diff = curData - pred1D;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity){
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[0] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[0])>realPrecision){	
					type[index] = 0;
					P0[0] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}			
			}
			else{
				type[index] = 0;
				P0[0] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}				
									
		/* Process row-i data 1 --> r2-1*/
		for (j = 1; j < r2; j++)
		{
			index = i*r2+j;
			curData = cur_data_pos[j];
			if(fabs(curData - mean[0]) <= realPrecision){
				type[index] = 1;
				P0[j] = mean[0];
			}
			else{
				pred2D = P0[j-1] + P1[j] - P1[j-1];
				diff = curData - pred2D;
				itvNum = fabs(diff)/realPrecision + 1;
				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
					
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[j])>realPrecision)
					{	
						type[index] = 0;
						P0[j] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}				
				}
				else
				{
					type[index] = 0;
					P0[j] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
		}
		cur_data_pos += dim0_offset;

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}
	return unpredictable_count;
}

unsigned char * SZ_compress_float_2D_MDQ_nonblocked(float *oriData, size_t r1, size_t r2, double realPrecision, size_t * comp_size){
	unsigned int quantization_intervals;
	float dense_pos =0.532267;
	if(optQuantMode==1)
	{
		quantization_intervals = 256;
		// dense_pos = 1;
		// quantization_intervals = optimize_intervals_and_compute_dense_position_float_2D(oriData, r1, r2, realPrecision, &dense_pos);
		printf("number of bins: %d\nerror bound %.4f dense position %.4f\n", quantization_intervals, realPrecision, dense_pos);
		//dense_pos = realPrecision;
		//dense_pos = 2.5867;
		//quantization_intervals = 512;
		// if(quantization_intervals < 128) quantization_intervals = 128;
		updateQuantizationInfo(quantization_intervals);
		intvCapacity = quantization_intervals - 2;
	}	
	else{
		quantization_intervals = intvCapacity;
		intvCapacity = quantization_intervals - 2;
	}

	// calculate block dims
	size_t num_blocks = 1;
	size_t num_elements = r1 * r2;

	size_t dim0_offset = r2;
	
	float *P0, *P1; // buffer
	size_t buffer_size = r2 * sizeof(float);
	P0 = (float *) malloc(buffer_size);
	P1 = (float *) malloc(buffer_size);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	// int unpred_data_max_size = ((int)(num_block_elements * 0.2) + 1) ;
	size_t unpred_data_max_size = num_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);

	// int unpredictable_count = 0;
	size_t total_unpred = 0;
	size_t index = 0;
	// NOTE: Currently max unpred count cannot exceed unsigned int limit
	unsigned int max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	float * unpredictable_data = result_unpredictable_data;
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	float mean;
	unsigned int unpredictable_count = SZ_compress_float_2D_MDQ_RA_block_2D_pred(data_pos, &mean, dense_pos, r1, r2, r1, r2, realPrecision, P0, P1, type, unpredictable_data);
	if(unpredictable_count > max_unpred_count){
		max_unpred_count = unpredictable_count;
	}
	total_unpred += unpredictable_count;
	
	//debug
	size_t flushed_count = 0;
	for(size_t i = 0;i<num_elements;i++){
		if(result_type[i] == 1){
			flushed_count ++;
		}
	}
	printf("flushed count: %d\n", flushed_count);
	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	// fflush(stdout);
	free(P0);
	free(P1);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	size_t totalEncodeSize = 0;
	memcpy(result_pos, &unpredictable_count, sizeof(unsigned int));
	result_pos += num_blocks * sizeof(unsigned int);
	size_t unpredictableEncodeSize = unpredictable_count * sizeof(float);
	memcpy(result_pos, unpredictable_data, unpredictableEncodeSize);
	result_pos += unpredictableEncodeSize;
	memcpy(result_pos, &mean, sizeof(float));
	result_pos += num_blocks * sizeof(float);
	
	encode(type, num_elements, result_pos, &enCodeSize);
	result_pos += enCodeSize;

	printf("type array size: %ld\n", enCodeSize);
	totalEncodeSize = result_pos - result;
	printf("Total size %ld\n", totalEncodeSize);
	free(result_unpredictable_data);
	free(result_type);
	SZ_ReleaseHuffman();

	*comp_size = totalEncodeSize;
	return result;
}

unsigned char * SZ_compress_float_2D_MDQ_RA(float *oriData, size_t r1, size_t r2, double realPrecision, size_t * comp_size){

	unsigned int quantization_intervals;
	float dense_pos =0.532267;
	if(optQuantMode==1)
	{
		quantization_intervals = 256;
		// dense_pos = 1;
		// quantization_intervals = optimize_intervals_and_compute_dense_position_float_2D(oriData, r1, r2, realPrecision, &dense_pos);
		printf("number of bins: %d\nerror bound %.4f dense position %.4f\n", quantization_intervals, realPrecision, dense_pos);
		updateQuantizationInfo(quantization_intervals);
		intvCapacity = quantization_intervals - 2;
	}	
	else{
		quantization_intervals = intvCapacity;
		intvCapacity = quantization_intervals - 2;
	}

	// calculate block dims
	size_t num_x, num_y;
	size_t block_size = 128;
	COMPUTE_2D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	COMPUTE_2D_NUMBER_OF_BLOCKS(r2, num_y, block_size);

	size_t split_index_x, split_index_y;
	size_t early_blockcount_x, early_blockcount_y;
	size_t late_blockcount_x, late_blockcount_y;
	COLL_BASE_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y;
	size_t num_blocks = num_x * num_y;
	size_t num_elements = r1 * r2;
	size_t dim0_offset = r2;

	float *P0, *P1; // buffer
	size_t buffer_size = early_blockcount_y * sizeof(float);
	P0 = (float *) malloc(buffer_size);
	P1 = (float *) malloc(buffer_size);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	// int unpred_data_max_size = ((int)(num_block_elements * 0.2) + 1) ;
	// NOTE: change unpred_data_max_size to DFA or not?
	size_t unpred_data_max_size = max_num_block_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);

	// NOTE: unpredictable count is limited to unsigned short in each block
	unsigned short * unpredictable_count = (unsigned short *) malloc(num_blocks * sizeof(int));
	float * mean = malloc(num_blocks * sizeof(float));
	size_t total_unpred = 0;
	size_t index = 0;
	unsigned short max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	float * unpredictable_data = result_unpredictable_data;
	size_t offset_x, offset_y;
	size_t current_blockcount_x, current_blockcount_y;
	size_t type_offset = 0;
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
			offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
			data_pos = oriData + offset_x * dim0_offset + offset_y;

			current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
			current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
			type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x;
			type = result_type + type_offset;
			// printf("i j k: %d %d %d, offset %ld %ld %ld type offset %ld\n", i, j, k, offset_x, offset_y, offset_z, type_offset);
			index = i * num_y + j;
			unpredictable_data = result_unpredictable_data + index * unpred_data_max_size;
			unpredictable_count[index] = SZ_compress_float_2D_MDQ_RA_block_2D_pred(data_pos, mean + index, dense_pos, r1, r2, current_blockcount_x, current_blockcount_y, realPrecision, P0, P1, type, unpredictable_data);
			if(unpredictable_count[index] > max_unpred_count){
				max_unpred_count = unpredictable_count[index];
			}
			total_unpred += unpredictable_count[index];
		}
	}
	//debug
	size_t flushed_count = 0;
	for(size_t i = 0;i<num_elements;i++){
		if(result_type[i] == 1){
			flushed_count ++;
		}
	}
	printf("flushed count: %d\n", flushed_count);	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	fflush(stdout);
	free(P0);
	free(P1);
	// size_t typeArray_size;
	// unsigned char * typeArray;// = (unsigned char *) malloc(num_elements * sizeof(int));
	// encode_withTree(result_type, num_elements, &typeArray, &typeArray_size);
	// free(typeArray);
	// printf("typeArray_size: %ld\n", typeArray_size);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;
	intToBytes_bigEndian(result_pos, block_size);
	result_pos += 4;
	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	size_t unpredictableEncodeSize;
	size_t totalEncodeSize = 0;
	unsigned short * block_pos = (unsigned short *) result_pos;
	unsigned char * block_start_pos = NULL;
	result_pos += num_blocks * sizeof(unsigned short); // skip block size
	memcpy(result_pos, unpredictable_count, num_blocks * sizeof(unsigned short));
	result_pos += num_blocks * sizeof(unsigned short);
	memcpy(result_pos, mean, num_blocks * sizeof(float));
	result_pos += num_blocks * sizeof(float);

	// printf("compress offset to start: %ld\n", result_pos - result);
	// fflush(stdout);
	size_t current_block_elements;
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			// printf("i j k: %d %d %d\n", i, j, k);
			index = i * num_y + j;
			block_start_pos = result_pos;
			// memcpy(result_pos, mean+index, 4);
			// result_pos += 4;

			if(unpredictable_count[index] > 0){
				unpredictable_data = result_unpredictable_data + index * unpred_data_max_size;
				unpredictableEncodeSize = unpredictable_count[index] * sizeof(float);
				memcpy(result_pos, unpredictable_data, unpredictableEncodeSize);
				result_pos += unpredictableEncodeSize;
			}
			offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
			offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;

			current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
			current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
			type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x;

			type = result_type + type_offset;
			current_block_elements = current_blockcount_x * current_blockcount_y;
			enCodeSize = 0;
			encode(type, current_block_elements, result_pos, &enCodeSize);
			// if(i == 1){
			// 	printf("dist to start %ld\n", result_pos - result);
			// }

			result_pos += enCodeSize;
			*block_pos = result_pos - block_start_pos;
			block_pos ++;
		}
	}
	totalEncodeSize = result_pos - result;
	// printf("Total size %ld\n", totalEncodeSize);
	free(mean);
	free(result_unpredictable_data);
	free(unpredictable_count);
	free(result_type);
	SZ_ReleaseHuffman();

	*comp_size = totalEncodeSize;
	return result;
}

size_t SZ_compress_float_3D_MDQ_RA_block_2_layers(float * block_ori_data, float * mean, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * P0, float * P1, float * P_, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;

	// data_pos = block_ori_data;
	// for(size_t i=0; i<block_dim_0; i++){
	// 	for(size_t j=0; j<block_dim_1; j++){
	// 		for(size_t k=0; k<block_dim_2; k++){
	// 			sum += *data_pos;
	// 			data_pos ++;
	// 		}
	// 		data_pos += dim1_offset - block_dim_2;
	// 	}
	// 	data_pos += dim0_offset - block_dim_1 * dim1_offset;
	// }
	// size_t num_elements = block_dim_0 * block_dim_1 * block_dim_2;
	// if(num_elements > 0) mean[0] = sum / num_elements;
	// else mean[0] = 0.0;
	mean[0] = block_ori_data[0];

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = block_ori_data;
	float curData;
	float pred1D, pred2D, pred3D;
	double itvNum;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	// Process Row-0 data 0
	pred1D = mean[0];
	curData = cur_data_pos[0];
	diff = curData - pred1D;
	itvNum = fabs(diff)/realPrecision + 1;
	if (itvNum < intvCapacity){
		if (diff < 0) itvNum = -itvNum;
		type[0] = (int) (itvNum/2) + intvRadius;
		P1[0] = pred1D + 2 * (type[0] - intvRadius) * realPrecision;
		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(curData-P1[1])>realPrecision){	
			type[0] = 0;
			P1[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}		
	}
	else{
		type[0] = 0;
		P1[0] = curData;
		unpredictable_data[unpredictable_count ++] = curData;
	}

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	curData = cur_data_pos[1];
	diff = curData - pred1D;
	itvNum = fabs(diff)/realPrecision + 1;
	if (itvNum < intvCapacity){
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(curData-P1[1])>realPrecision){	
			type[1] = 0;
			P1[1] = curData;	
			unpredictable_data[unpredictable_count ++] = curData;
		}		
	}
	else{
		type[1] = 0;
		P1[1] = curData;
		unpredictable_data[unpredictable_count ++] = curData;
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*P1[j-1] - P1[j-2];
		curData = cur_data_pos[j];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[j])>realPrecision){	
				type[j] = 0;
				P1[j] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else{
			type[j] = 0;
			P1[j] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}
	cur_data_pos += dim1_offset;

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = P1[index-r3];
		curData = cur_data_pos[0];
		diff = curData - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P1[index] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[index])>realPrecision)
			{	
				type[index] = 0;
				P1[index] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else
		{
			type[index] = 0;
			P1[index] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}

		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];

			curData = cur_data_pos[j];
			diff = curData - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[index])>realPrecision)
				{	
					type[index] = 0;
					P1[index] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}				
			}
			else
			{
				type[index] = 0;
				P1[index] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = P1[0];
		curData = cur_data_pos[0];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P0[0])>realPrecision)
			{	
				type[index] = 0;
				P0[0] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else
		{
			type[index] = 0;
			P0[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = P0[j-1] + P1[j] - P1[j-1];
			curData = cur_data_pos[j];
			diff = curData - pred2D;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[j])>realPrecision)
				{	
					type[index] = 0;
					P0[j] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
			else
			{
				type[index] = 0;
				P0[j] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}

		cur_data_pos += dim1_offset;
	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
			curData = cur_data_pos[0];
			diff = curData - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[index2D])>realPrecision)
				{	
					type[index] = 0;
					P0[index2D] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}				
			}
			else
			{
				type[index] = 0;
				P0[index2D] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				if(i >= 2 && j>=2 && k >=2){
					// pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];

					pred3D = - P_[index2D - 2*r3 - 2]
							 + 2 * (P_[index2D - 2*r3 - 1] + P_[index2D - r3 - 2] + P1[index2D - 2*r3 - 2])
							 - 4 * (P_[index2D - r3 - 1] + P1[index2D - 2*r3 - 1] + P1[index2D - r3 - 2])
							 - (P_[index2D] + P1[index2D - 2*r3] + P0[index2D - 2])
							 - (P_[index2D - 2*r3] + P_[index2D - 2] + P0[index2D - 2*r3 - 2])
							 - 4 * (P1[index2D - r3] + P1[index2D - 1] + P0[index2D - r3 - 1])
							 + 2 * (P1[index2D] + P0[index2D - r3] + P0[index2D - 1])
							 + 2 * (P_[index2D - r3] + P_[index2D - 1] + P1[index2D - 2*r3] + P1[index2D - 2] + P0[index2D - 2*r3 - 1] + P0[index2D - r3 -2])
							 + 8 * (P1[index2D - r3 - 1]);
					// if(i==2 && j==2 && k==2){
					// 	printf("pred3D: %.4f\n", pred3D);
					// 	float * tmp_P;
					// 	for(int i=0; i<3; i++){
					// 		if(i==0) tmp_P = P0;
					// 		else if(i==1) tmp_P = P1;
					// 		else tmp_P = P_;
					// 		for(int j=0; j<3; j++){
					// 			for(int k=0; k<3; k++){
					// 				printf("%.4f ", tmp_P[index2D - j*r3 - k]);
					// 			}
					// 		}
					// 	}
					// 	printf("\n");

					// }

				}
				else{
					pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
				}
				curData = cur_data_pos[j];
				diff = curData - pred3D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
					
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[index2D])>realPrecision)
					{	
						type[index] = 0;
						P0[index2D] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}					
				}
				else
				{
					type[index] = 0;
					P0[index2D] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
		float *Pt;
		// Pt = P1;
		// P1 = P0;
		// P0 = Pt;
		Pt = P_;
		P_ = P1;
		P1 = P0;
		P0 = Pt;
	}

	return unpredictable_count;
}

float SZ_compress_float_3D_MDQ_RA_block_no_mean(float * block_ori_data, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * P0, float * P1, int * type, unsigned short * unpred_count, float * unpredictable_data){

	float pred_err = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;

	// data_pos = block_ori_data;
	// for(size_t i=0; i<block_dim_0; i++){
	// 	for(size_t j=0; j<block_dim_1; j++){
	// 		for(size_t k=0; k<block_dim_2; k++){
	// 			sum += *data_pos;
	// 			data_pos ++;
	// 		}
	// 		data_pos += dim1_offset - block_dim_2;
	// 	}
	// 	data_pos += dim0_offset - block_dim_1 * dim1_offset;
	// }
	// size_t num_elements = block_dim_0 * block_dim_1 * block_dim_2;
	// if(num_elements > 0) mean[0] = sum / num_elements;
	// else mean[0] = 0.0;

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = block_ori_data;
	float curData;
	float pred1D, pred2D, pred3D;
	double itvNum;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	// Process Row-0 data 0
	// make it unpredictable
	// pred1D = mean[0];
	// curData = cur_data_pos[0];
	// diff = curData - pred1D;
	// pred_err += diff * diff;
	// itvNum = fabs(diff)/realPrecision + 1;
	// if (itvNum < intvCapacity){
	// 	if (diff < 0) itvNum = -itvNum;
	// 	type[0] = (int) (itvNum/2) + intvRadius;
	// 	P1[0] = pred1D + 2 * (type[0] - intvRadius) * realPrecision;
	// 	//ganrantee comporession error against the case of machine-epsilon
	// 	if(fabs(curData-P1[1])>realPrecision){	
	// 		type[0] = 0;
	// 		P1[0] = curData;
	// 		unpredictable_data[unpredictable_count ++] = curData;
	// 	}		
	// }
	// else
	curData = cur_data_pos[0];
	{
		type[0] = 0;
		P1[0] = curData;
		unpredictable_data[unpredictable_count ++] = curData;
	}

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	curData = cur_data_pos[1];
	diff = curData - pred1D;
	pred_err += diff * diff;
	itvNum = fabs(diff)/realPrecision + 1;
	if (itvNum < intvCapacity){
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(curData-P1[1])>realPrecision){	
			type[1] = 0;
			P1[1] = curData;	
			unpredictable_data[unpredictable_count ++] = curData;
		}		
	}
	else{
		type[1] = 0;
		P1[1] = curData;
		unpredictable_data[unpredictable_count ++] = curData;
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*P1[j-1] - P1[j-2];
		curData = cur_data_pos[j];
		diff = curData - pred1D;
		pred_err += diff * diff;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[j])>realPrecision){	
				type[j] = 0;
				P1[j] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else{
			type[j] = 0;
			P1[j] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}
	cur_data_pos += dim1_offset;

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = P1[index-r3];
		curData = cur_data_pos[0];
		diff = curData - pred1D;
		pred_err += diff * diff;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P1[index] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[index])>realPrecision)
			{	
				type[index] = 0;
				P1[index] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else
		{
			type[index] = 0;
			P1[index] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}

		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];

			curData = cur_data_pos[j];
			diff = curData - pred2D;
			pred_err += diff * diff;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[index])>realPrecision)
				{	
					type[index] = 0;
					P1[index] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}				
			}
			else
			{
				type[index] = 0;
				P1[index] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = P1[0];
		curData = cur_data_pos[0];
		diff = curData - pred1D;
		pred_err += diff * diff;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P0[0])>realPrecision)
			{	
				type[index] = 0;
				P0[0] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else
		{
			type[index] = 0;
			P0[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = P0[j-1] + P1[j] - P1[j-1];
			curData = cur_data_pos[j];
			diff = curData - pred2D;
			pred_err += diff * diff;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[j])>realPrecision)
				{	
					type[index] = 0;
					P0[j] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
			else
			{
				type[index] = 0;
				P0[j] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}

		cur_data_pos += dim1_offset;
	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
			curData = cur_data_pos[0];
			diff = curData - pred2D;
			pred_err += diff * diff;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[index2D])>realPrecision)
				{	
					type[index] = 0;
					P0[index2D] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}				
			}
			else
			{
				type[index] = 0;
				P0[index2D] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
				curData = cur_data_pos[j];
				diff = curData - pred3D;
				pred_err += diff * diff;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
					
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[index2D])>realPrecision)
					{	
						type[index] = 0;
						P0[index2D] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}					
				}
				else
				{
					type[index] = 0;
					P0[index2D] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}
	*unpred_count = unpredictable_count;
	return pred_err;
}

size_t SZ_compress_float_3D_MDQ_RA_block(float * block_ori_data, float * mean, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * P0, float * P1, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;

	// data_pos = block_ori_data;
	// for(size_t i=0; i<block_dim_0; i++){
	// 	for(size_t j=0; j<block_dim_1; j++){
	// 		for(size_t k=0; k<block_dim_2; k++){
	// 			sum += *data_pos;
	// 			data_pos ++;
	// 		}
	// 		data_pos += dim1_offset - block_dim_2;
	// 	}
	// 	data_pos += dim0_offset - block_dim_1 * dim1_offset;
	// }
	// size_t num_elements = block_dim_0 * block_dim_1 * block_dim_2;
	// if(num_elements > 0) mean[0] = sum / num_elements;
	// else mean[0] = 0.0;
	mean[0] = block_ori_data[0];

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = block_ori_data;
	float curData;
	float pred1D, pred2D, pred3D;
	double itvNum;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	// Process Row-0 data 0
	pred1D = mean[0];
	curData = cur_data_pos[0];
	diff = curData - pred1D;
	itvNum = fabs(diff)/realPrecision + 1;
	if (itvNum < intvCapacity){
		if (diff < 0) itvNum = -itvNum;
		type[0] = (int) (itvNum/2) + intvRadius;
		P1[0] = pred1D + 2 * (type[0] - intvRadius) * realPrecision;
		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(curData-P1[1])>realPrecision){	
			type[0] = 0;
			P1[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}		
	}
	else{
		type[0] = 0;
		P1[0] = curData;
		unpredictable_data[unpredictable_count ++] = curData;
	}

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	curData = cur_data_pos[1];
	diff = curData - pred1D;
	itvNum = fabs(diff)/realPrecision + 1;
	if (itvNum < intvCapacity){
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(curData-P1[1])>realPrecision){	
			type[1] = 0;
			P1[1] = curData;	
			unpredictable_data[unpredictable_count ++] = curData;
		}		
	}
	else{
		type[1] = 0;
		P1[1] = curData;
		unpredictable_data[unpredictable_count ++] = curData;
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*P1[j-1] - P1[j-2];
		curData = cur_data_pos[j];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[j])>realPrecision){	
				type[j] = 0;
				P1[j] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else{
			type[j] = 0;
			P1[j] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}
	cur_data_pos += dim1_offset;

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = P1[index-r3];
		curData = cur_data_pos[0];
		diff = curData - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P1[index] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[index])>realPrecision)
			{	
				type[index] = 0;
				P1[index] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else
		{
			type[index] = 0;
			P1[index] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}

		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];

			curData = cur_data_pos[j];
			diff = curData - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[index])>realPrecision)
				{	
					type[index] = 0;
					P1[index] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}				
			}
			else
			{
				type[index] = 0;
				P1[index] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = P1[0];
		curData = cur_data_pos[0];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P0[0])>realPrecision)
			{	
				type[index] = 0;
				P0[0] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else
		{
			type[index] = 0;
			P0[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = P0[j-1] + P1[j] - P1[j-1];
			curData = cur_data_pos[j];
			diff = curData - pred2D;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[j])>realPrecision)
				{	
					type[index] = 0;
					P0[j] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
			else
			{
				type[index] = 0;
				P0[j] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}

		cur_data_pos += dim1_offset;
	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
			curData = cur_data_pos[0];
			diff = curData - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[index2D])>realPrecision)
				{	
					type[index] = 0;
					P0[index2D] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}				
			}
			else
			{
				type[index] = 0;
				P0[index2D] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
				curData = cur_data_pos[j];
				diff = curData - pred3D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
					
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[index2D])>realPrecision)
					{	
						type[index] = 0;
						P0[index2D] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}					
				}
				else
				{
					type[index] = 0;
					P0[index2D] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	return unpredictable_count;
}

size_t SZ_compress_float_3D_MDQ_RA_block_adaptive(float * block_ori_data, float * mean, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * P0, float * P1, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;

	// data_pos = block_ori_data;
	// for(size_t i=0; i<block_dim_0; i++){
	// 	for(size_t j=0; j<block_dim_1; j++){
	// 		for(size_t k=0; k<block_dim_2; k++){
	// 			sum += *data_pos;
	// 			data_pos ++;
	// 		}
	// 		data_pos += dim1_offset - block_dim_2;
	// 	}
	// 	data_pos += dim0_offset - block_dim_1 * dim1_offset;
	// }
	// size_t num_elements = block_dim_0 * block_dim_1 * block_dim_2;
	// if(num_elements > 0) mean[0] = sum / num_elements;
	// else mean[0] = 0.0;
	mean[0] = block_ori_data[0];

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = block_ori_data;
	float curData;
	float pred1D, pred2D, pred3D;
	double itvNum;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	// Process Row-0 data 0
	pred1D = mean[0];
	curData = cur_data_pos[0];
	diff = curData - pred1D;
	itvNum = fabs(diff)/realPrecision + 1;
	if (itvNum < intvCapacity){
		if (diff < 0) itvNum = -itvNum;
		type[0] = (int) (itvNum/2) + intvRadius;
		P1[0] = pred1D + 2 * (type[0] - intvRadius) * realPrecision;
		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(curData-P1[1])>realPrecision){	
			type[0] = 0;
			P1[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}		
	}
	else{
		type[0] = 0;
		P1[0] = curData;
		unpredictable_data[unpredictable_count ++] = curData;
	}

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	curData = cur_data_pos[1];
	diff = curData - pred1D;
	itvNum = fabs(diff)/realPrecision + 1;
	if (itvNum < intvCapacity){
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(curData-P1[1])>realPrecision){	
			type[1] = 0;
			P1[1] = curData;	
			unpredictable_data[unpredictable_count ++] = curData;
		}		
	}
	else{
		type[1] = 0;
		P1[1] = curData;
		unpredictable_data[unpredictable_count ++] = curData;
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*P1[j-1] - P1[j-2];
		curData = cur_data_pos[j];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[j])>realPrecision){	
				type[j] = 0;
				P1[j] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else{
			type[j] = 0;
			P1[j] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}
	cur_data_pos += dim1_offset;

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = P1[index-r3];
		curData = cur_data_pos[0];
		diff = curData - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P1[index] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[index])>realPrecision)
			{	
				type[index] = 0;
				P1[index] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else
		{
			type[index] = 0;
			P1[index] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}

		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];

			curData = cur_data_pos[j];
			diff = curData - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[index])>realPrecision)
				{	
					type[index] = 0;
					P1[index] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}				
			}
			else
			{
				type[index] = 0;
				P1[index] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = P1[0];
		curData = cur_data_pos[0];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P0[0])>realPrecision)
			{	
				type[index] = 0;
				P0[0] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else
		{
			type[index] = 0;
			P0[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = P0[j-1] + P1[j] - P1[j-1];
			curData = cur_data_pos[j];
			diff = curData - pred2D;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[j])>realPrecision)
				{	
					type[index] = 0;
					P0[j] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
			else
			{
				type[index] = 0;
				P0[j] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}

		cur_data_pos += dim1_offset;
	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
			curData = cur_data_pos[0];
			diff = curData - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[index2D])>realPrecision)
				{	
					type[index] = 0;
					P0[index2D] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}				
			}
			else
			{
				type[index] = 0;
				P0[index2D] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
				double tmp_realPrecision = realPrecision * (1 - ((j) % 10) * 0.1);
				if((i + j + k ) % 10 == 0){
					tmp_realPrecision = realPrecision * 0.5;
				}
				else{
					tmp_realPrecision = realPrecision;
				}

				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
				curData = cur_data_pos[j];
				diff = curData - pred3D;

				itvNum = fabs(diff)/tmp_realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * tmp_realPrecision;
					
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[index2D])>tmp_realPrecision)
					{	
						type[index] = 0;
						P0[index2D] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}					
				}
				else
				{
					type[index] = 0;
					P0[index2D] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}

			}
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	return unpredictable_count;
}

size_t SZ_compress_float_3D_MDQ_RA_block_3D_pred(float * block_ori_data, float * mean, float dense_pos, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * P0, float * P1, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;

	data_pos = block_ori_data;
	size_t mean_count = 0;
	for(size_t i=0; i<block_dim_0; i++){
		for(size_t j=0; j<block_dim_1; j++){
			for(size_t k=0; k<block_dim_2; k++){
				if(fabs(*data_pos - dense_pos) <= realPrecision){
					sum += *data_pos;
					mean_count ++;
				}
				data_pos ++;
			}
			data_pos += dim1_offset - block_dim_2;
		}
		data_pos += dim0_offset - block_dim_1 * dim1_offset;
	}
	size_t num_elements = block_dim_0 * block_dim_1 * block_dim_2;
	if(mean_count > 0) mean[0] = sum / mean_count;
	else mean[0] = 0;

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = block_ori_data;
	float curData;
	float pred1D, pred2D, pred3D;
	double itvNum;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	// Process Row-0 data 0
	curData = cur_data_pos[0];
	if(fabs(curData - mean[0]) <= realPrecision){
		type[0] = 1;
		P1[0] = mean[0];
	}
	else{
		pred1D = mean[0];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type[0] = (int) (itvNum/2) + intvRadius;
			P1[0] = pred1D + 2 * (type[0] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[0])>realPrecision){	
				type[0] = 0;
				P1[0] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}		
		}
		else{
			type[0] = 0;
			P1[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}

	/* Process Row-0 data 1*/
	curData = cur_data_pos[1];
	if(fabs(curData - mean[0]) <= realPrecision){
		type[1] = 1;
		P1[1] = mean[0];
	}
	else{
		pred1D = P1[0];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type[1] = (int) (itvNum/2) + intvRadius;
			P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[1])>realPrecision){	
				type[1] = 0;
				P1[1] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}		
		}
		else{
			type[1] = 0;
			P1[1] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		curData = cur_data_pos[j];
		if(fabs(curData - mean[0]) <= realPrecision){
			type[j] = 1;
			P1[j] = mean[0];
		}
		else{
			pred1D = 2*P1[j-1] - P1[j-2];
			diff = curData - pred1D;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity){
				if (diff < 0) itvNum = -itvNum;
				type[j] = (int) (itvNum/2) + intvRadius;
				P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[j])>realPrecision){	
					type[j] = 0;
					P1[j] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}			
			}
			else{
				type[j] = 0;
				P1[j] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
	}
	cur_data_pos += dim1_offset;

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;
		curData = cur_data_pos[0];
		if(fabs(curData - mean[0]) <= realPrecision){
			type[index] = 1;
			P1[index] = mean[0];
		}
		else{
			pred1D = P1[index-r3];
			diff = curData - pred1D;
			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
				
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[index])>realPrecision)
				{	
					type[index] = 0;
					P1[index] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}			
			}
			else
			{
				type[index] = 0;
				P1[index] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			curData = cur_data_pos[j];
			if(fabs(curData - mean[0]) <= realPrecision){
				type[index] = 1;
				P1[index] = mean[0];
			}
			else{
				pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];
				diff = curData - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P1[index] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
					
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P1[index])>realPrecision)
					{	
						type[index] = 0;
						P1[index] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}				
				}
				else
				{
					type[index] = 0;
					P1[index] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
		}
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r23;
		curData = cur_data_pos[0];
		if(fabs(curData - mean[0]) <= realPrecision){
			type[index] = 1;
			P0[0] = mean[0];
		}
		else{
			pred1D = P1[0];
			diff = curData - pred1D;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[0])>realPrecision)
				{	
					type[index] = 0;
					P0[0] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}			
			}
			else
			{
				type[index] = 0;
				P0[0] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			curData = cur_data_pos[j];
			if(fabs(curData - mean[0]) <= realPrecision){
				type[index] = 1;
				P0[j] = mean[0];
			}
			else{
				pred2D = P0[j-1] + P1[j] - P1[j-1];
				diff = curData - pred2D;
				itvNum = fabs(diff)/realPrecision + 1;
				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[j])>realPrecision)
					{	
						type[index] = 0;
						P0[j] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}
				}
				else
				{
					type[index] = 0;
					P0[j] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
		}

		cur_data_pos += dim1_offset;
	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			curData = cur_data_pos[0];
			index2D = i*r3;
			if(fabs(curData - mean[0]) <= realPrecision){
				type[index] = 1;
				P0[index2D] = mean[0];
			}
			else{			
				pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
				diff = curData - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[index2D])>realPrecision)
					{	
						type[index] = 0;
						P0[index2D] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}				
				}
				else
				{
					type[index] = 0;
					P0[index2D] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				curData = cur_data_pos[j];
				if(fabs(curData - mean[0]) <= realPrecision){
					type[index] = 1;
					P0[index2D] = mean[0];
				}
				else{
					pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
					diff = curData - pred3D;

					itvNum = fabs(diff)/realPrecision + 1;

					if (itvNum < intvCapacity)
					{
						if (diff < 0) itvNum = -itvNum;
						type[index] = (int) (itvNum/2) + intvRadius;
						P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
						
						//ganrantee comporession error against the case of machine-epsilon
						if(fabs(curData-P0[index2D])>realPrecision)
						{	
							type[index] = 0;
							P0[index2D] = curData;	
							unpredictable_data[unpredictable_count ++] = curData;
						}					
					}
					else
					{
						type[index] = 0;
						P0[index2D] = curData;
						unpredictable_data[unpredictable_count ++] = curData;
					}
				}
			}
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	return unpredictable_count;
}

size_t SZ_compress_float_3D_MDQ_RA_block_3D_pred_multi_means(float * block_ori_data, unsigned int mean_count, float * mean, float dense_pos, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * P0, float * P1, int * type, float * unpredictable_data){

	float sum = 0.0;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = block_ori_data;
	float curData;
	float pred1D, pred2D, pred3D;
	double itvNum;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	unsigned int mean_index;
	float mean_start = mean[0] - realPrecision;
	float mean_end = mean[mean_count - 1] + realPrecision;
	float mean_res;
	// Process Row-0 data 0
	curData = cur_data_pos[0];
	if(curData >= mean_start && curData <= mean_end){
		mean_index = ((curData - mean[0])/realPrecision + 1) / 2;
		if(mean_index >= mean_count) mean_index = mean_count - 1;
		while(mean_index < mean_count - 1 && mean[mean_index] < curData) mean_index ++;
		while(mean_index > 0 && mean[mean_index] > curData) mean_index --;
		if(mean_index < mean_count -1 && fabs(curData - mean[mean_index]) > fabs(curData - mean[mean_index+1])){
			mean_index ++;
		}
		type[0] = mean_index + 1;
		P1[0] = mean[mean_index];
	}
	else{
		pred1D = dense_pos;
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type[0] = (int) (itvNum/2) + intvRadius;
			P1[0] = pred1D + 2 * (type[0] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[0])>realPrecision){	
				type[0] = 0;
				P1[0] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}		
		}
		else{
			type[0] = 0;
			P1[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}

	/* Process Row-0 data 1*/
	curData = cur_data_pos[1];
	if(curData >= mean_start && curData <= mean_end){
		mean_index = ((curData - mean[0])/realPrecision + 1) / 2;
		if(mean_index >= mean_count) mean_index = mean_count - 1;
		while(mean_index < mean_count - 1 && mean[mean_index] < curData) mean_index ++;
		while(mean_index > 0 && mean[mean_index] > curData) mean_index --;
		if(mean_index < mean_count -1 && fabs(curData - mean[mean_index]) > fabs(curData - mean[mean_index+1])){
			mean_index ++;
		}
		type[1] = mean_index + 1;
		P1[1] = mean[mean_index];
	}
	else{
		pred1D = P1[0];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type[1] = (int) (itvNum/2) + intvRadius;
			P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[1])>realPrecision){	
				type[1] = 0;
				P1[1] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}		
		}
		else{
			type[1] = 0;
			P1[1] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		curData = cur_data_pos[j];
		if(curData >= mean_start && curData <= mean_end){
			mean_index = ((curData - mean[0])/realPrecision + 1) / 2;
			if(mean_index >= mean_count) mean_index = mean_count - 1;
			while(mean_index < mean_count - 1 && mean[mean_index] < curData) mean_index ++;
			while(mean_index > 0 && mean[mean_index] > curData) mean_index --;
			if(mean_index < mean_count -1 && fabs(curData - mean[mean_index]) > fabs(curData - mean[mean_index+1])){
				mean_index ++;
			}
			type[j] = mean_index + 1;
			P1[j] = mean[mean_index];
		}
		else{
			pred1D = 2*P1[j-1] - P1[j-2];
			diff = curData - pred1D;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity){
				if (diff < 0) itvNum = -itvNum;
				type[j] = (int) (itvNum/2) + intvRadius;
				P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[j])>realPrecision){	
					type[j] = 0;
					P1[j] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}			
			}
			else{
				type[j] = 0;
				P1[j] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
	}
	cur_data_pos += dim1_offset;

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;
		curData = cur_data_pos[0];
		if(curData >= mean_start && curData <= mean_end){
			mean_index = ((curData - mean[0])/realPrecision + 1) / 2;
			if(mean_index >= mean_count) mean_index = mean_count - 1;
			while(mean_index < mean_count - 1 && mean[mean_index] < curData) mean_index ++;
			while(mean_index > 0 && mean[mean_index] > curData) mean_index --;
			if(mean_index < mean_count -1 && fabs(curData - mean[mean_index]) > fabs(curData - mean[mean_index+1])){
				mean_index ++;
			}
			type[index] = mean_index + 1;
			P1[index] = mean[mean_index];
		}
		else{
			pred1D = P1[index-r3];
			diff = curData - pred1D;
			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
				
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[index])>realPrecision)
				{	
					type[index] = 0;
					P1[index] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}			
			}
			else
			{
				type[index] = 0;
				P1[index] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			curData = cur_data_pos[j];
			if(curData >= mean_start && curData <= mean_end){
				mean_index = ((curData - mean[0])/realPrecision + 1) / 2;
				if(mean_index >= mean_count) mean_index = mean_count - 1;
				while(mean_index < mean_count - 1 && mean[mean_index] < curData) mean_index ++;
				while(mean_index > 0 && mean[mean_index] > curData) mean_index --;
				if(mean_index < mean_count -1 && fabs(curData - mean[mean_index]) > fabs(curData - mean[mean_index+1])){
					mean_index ++;
				}
				type[index] = mean_index + 1;
				P1[index] = mean[mean_index];
			}
			else{
				pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];
				diff = curData - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P1[index] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
					
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P1[index])>realPrecision)
					{	
						type[index] = 0;
						P1[index] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}				
				}
				else
				{
					type[index] = 0;
					P1[index] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
		}
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r23;
		curData = cur_data_pos[0];
		if(curData >= mean_start && curData <= mean_end){
			mean_index = ((curData - mean[0])/realPrecision + 1) / 2;
			if(mean_index >= mean_count) mean_index = mean_count - 1;
			while(mean_index < mean_count - 1 && mean[mean_index] < curData) mean_index ++;
			while(mean_index > 0 && mean[mean_index] > curData) mean_index --;
			if(mean_index < mean_count -1 && fabs(curData - mean[mean_index]) > fabs(curData - mean[mean_index+1])){
				mean_index ++;
			}
			type[index] = mean_index + 1;
			P0[0] = mean[mean_index];
		}
		else{
			pred1D = P1[0];
			diff = curData - pred1D;
			itvNum = fabs(diff)/realPrecision + 1;
			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[0])>realPrecision)
				{	
					type[index] = 0;
					P0[0] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}			
			}
			else
			{
				type[index] = 0;
				P0[0] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			curData = cur_data_pos[j];
			if(curData >= mean_start && curData <= mean_end){
				mean_index = ((curData - mean[0])/realPrecision + 1) / 2;
				if(mean_index >= mean_count) mean_index = mean_count - 1;
				while(mean_index < mean_count - 1 && mean[mean_index] < curData) mean_index ++;
				while(mean_index > 0 && mean[mean_index] > curData) mean_index --;
				if(mean_index < mean_count -1 && fabs(curData - mean[mean_index]) > fabs(curData - mean[mean_index+1])){
					mean_index ++;
				}
				type[index] = mean_index + 1;
				P0[j] = mean[mean_index];
			}
			else{
				pred2D = P0[j-1] + P1[j] - P1[j-1];
				diff = curData - pred2D;
				itvNum = fabs(diff)/realPrecision + 1;
				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[j])>realPrecision)
					{	
						type[index] = 0;
						P0[j] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}
				}
				else
				{
					type[index] = 0;
					P0[j] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
		}

		cur_data_pos += dim1_offset;
	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			curData = cur_data_pos[0];
			index2D = i*r3;
			if(curData >= mean_start && curData <= mean_end){
				mean_index = ((curData - mean[0])/realPrecision + 1) / 2;
				if(mean_index >= mean_count) mean_index = mean_count - 1;
				while(mean_index < mean_count - 1 && mean[mean_index] < curData) mean_index ++;
				while(mean_index > 0 && mean[mean_index] > curData) mean_index --;
				if(mean_index < mean_count -1 && fabs(curData - mean[mean_index]) > fabs(curData - mean[mean_index+1])){
					mean_index ++;
				}
				type[index] = mean_index + 1;
				P0[index2D] = mean[mean_index];
			}
			else{			
				pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
				diff = curData - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[index2D])>realPrecision)
					{	
						type[index] = 0;
						P0[index2D] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}				
				}
				else
				{
					type[index] = 0;
					P0[index2D] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				curData = cur_data_pos[j];
				if(curData >= mean_start && curData <= mean_end){
					mean_index = ((curData - mean[0])/realPrecision + 1) / 2;
					if(mean_index >= mean_count) mean_index = mean_count - 1;
					while(mean_index < mean_count - 1 && mean[mean_index] < curData) mean_index ++;
					while(mean_index > 0 && mean[mean_index] > curData) mean_index --;
					if(mean_index < mean_count -1 && fabs(curData - mean[mean_index]) > fabs(curData - mean[mean_index+1])){
						mean_index ++;
					}
					type[index] = mean_index + 1;
					P0[index2D] = mean[mean_index];
				}
				else{
					pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
					diff = curData - pred3D;

					itvNum = fabs(diff)/realPrecision + 1;

					if (itvNum < intvCapacity)
					{
						if (diff < 0) itvNum = -itvNum;
						type[index] = (int) (itvNum/2) + intvRadius;
						P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
						
						//ganrantee comporession error against the case of machine-epsilon
						if(fabs(curData-P0[index2D])>realPrecision)
						{	
							type[index] = 0;
							P0[index2D] = curData;	
							unpredictable_data[unpredictable_count ++] = curData;
						}					
					}
					else
					{
						type[index] = 0;
						P0[index2D] = curData;
						unpredictable_data[unpredictable_count ++] = curData;
					}
				}
			}
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	return unpredictable_count;
}

size_t SZ_compress_float_3D_MDQ_RA_block_3D_pred_flush_after_compare(float * block_ori_data, float * mean, float dense_pos, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * P0, float * P1, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;

	data_pos = block_ori_data;
	size_t mean_count = 0;
	for(size_t i=0; i<block_dim_0; i++){
		for(size_t j=0; j<block_dim_1; j++){
			for(size_t k=0; k<block_dim_2; k++){
				if(fabs(*data_pos - dense_pos) <= realPrecision){
					sum += *data_pos;
					mean_count ++;
				}
				data_pos ++;
			}
			data_pos += dim1_offset - block_dim_2;
		}
		data_pos += dim0_offset - block_dim_1 * dim1_offset;
	}
	size_t num_elements = block_dim_0 * block_dim_1 * block_dim_2;
	if(mean_count > 0) mean[0] = sum / mean_count;
	else mean[0] = 0;

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = block_ori_data;
	float curData;
	float pred1D, pred2D, pred3D;
	double itvNum;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	// Process Row-0 data 0
	float pred_res;
	float mean_res;
	int type_;
	curData = cur_data_pos[0];
	pred1D = mean[0];
	diff = curData - pred1D;
	itvNum = fabs(diff)/realPrecision + 1;
	mean_res = mean[0];
	if (itvNum < intvCapacity){
		if (diff < 0) itvNum = -itvNum;
		type_ = (int) (itvNum/2) + intvRadius;
		pred_res = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		if(fabs(curData - mean_res) <= fabs(curData - pred_res)){
			type[0] = 1;
			P1[0] = mean_res;
		}	
		else{
			type[0] = type_;
			P1[0] = pred_res;
		}
		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(curData-P1[0])>realPrecision){	
			type[0] = 0;
			P1[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}		
	}
	else{
		if(fabs(curData - mean_res) <= realPrecision){
			type[0] = 1;
			P1[0] = mean_res;
		}
		else{
			type[0] = 0;
			P1[0] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}

	/* Process Row-0 data 1*/
	curData = cur_data_pos[1];
	pred1D = P1[0];
	diff = curData - pred1D;
	itvNum = fabs(diff)/realPrecision + 1;
	mean_res = mean[0];
	if (itvNum < intvCapacity){
		if (diff < 0) itvNum = -itvNum;
		type_ = (int) (itvNum/2) + intvRadius;
		pred_res = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		if(fabs(curData - mean_res) <= fabs(curData - pred_res)){
			type[1] = 1;
			P1[1] = mean_res;
		}	
		else{
			type[1] = type_;
			P1[1] = pred_res;
		}
		//ganrantee comporession error against the case of machine-epsilon
		if(fabs(curData-P1[1])>realPrecision){	
			type[1] = 0;
			P1[1] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}		
	}
	else{
		if(fabs(curData - mean_res) <= realPrecision){
			type[1] = 1;
			P1[1] = mean_res;
		}
		else{
			type[1] = 0;
			P1[1] = curData;
			unpredictable_data[unpredictable_count ++] = curData;
		}
	}

    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		curData = cur_data_pos[j];
		pred1D = 2*P1[j-1] - P1[j-2];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		mean_res = mean[0];
		if (itvNum < intvCapacity){
			if (diff < 0) itvNum = -itvNum;
			type_ = (int) (itvNum/2) + intvRadius;
			pred_res = pred1D + 2 * (type_ - intvRadius) * realPrecision;
			if(fabs(curData - mean_res) <= fabs(curData - pred_res)){
				type[j] = 1;
				P1[j] = mean_res;
			}
			else {
				type[j] = type_;
				P1[j] = pred_res;
			}
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[j])>realPrecision){	
				type[j] = 0;
				P1[j] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else{
			if(fabs(curData - mean_res) <= realPrecision){
				type[j] = 1;
				P1[j] = mean_res;
			}
			else{
				type[j] = 0;
				P1[j] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
	}
	cur_data_pos += dim1_offset;

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;
		curData = cur_data_pos[0];
		pred1D = P1[index-r3];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		mean_res = mean[0];
		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type_ = (int) (itvNum/2) + intvRadius;
			pred_res = pred1D + 2 * (type_ - intvRadius) * realPrecision;
			if(fabs(curData - mean_res) <= fabs(curData - pred_res)){
				type[index] = 1;
				P1[index] = mean_res;
			}
			else {
				type[index] = type_;
				P1[index] = pred_res;
			}
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P1[index])>realPrecision)
			{	
				type[index] = 0;
				P1[index] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else
		{
			if(fabs(curData - mean_res) <= realPrecision){
				type[index] = 1;
				P1[index] = mean_res;
			}
			else{
				type[index] = 0;
				P1[index] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}	
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			curData = cur_data_pos[j];
			pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];
			diff = curData - pred2D;
			itvNum = fabs(diff)/realPrecision + 1;
			mean_res = mean[0];
			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type_ = (int) (itvNum/2) + intvRadius;
				pred_res = pred2D + 2 * (type_ - intvRadius) * realPrecision;
				if(fabs(curData - mean_res) <= fabs(curData - pred_res)){
					type[index] = 1;
					P1[index] = mean_res;
				}
				else {
					type[index] = type_;
					P1[index] = pred_res;
				}
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P1[index])>realPrecision)
				{	
					type[index] = 0;
					P1[index] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}				
			}
			else
			{
				if(fabs(curData - mean_res) <= realPrecision){
					type[index] = 1;
					P1[index] = mean_res;
				}
				else{				
					type[index] = 0;
					P1[index] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
		}
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r23;
		curData = cur_data_pos[0];
		pred1D = P1[0];
		diff = curData - pred1D;
		itvNum = fabs(diff)/realPrecision + 1;
		mean_res = mean[0];
		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type_ = (int) (itvNum/2) + intvRadius;
			pred_res = pred1D + 2 * (type_ - intvRadius) * realPrecision;
			if(fabs(curData - mean_res) <= fabs(curData - pred_res)){
				type[index] = 1;
				P0[0] = mean_res;
			}
			else {
				type[index] = type_;
				P0[0] = pred_res;
			}
			//ganrantee comporession error against the case of machine-epsilon
			if(fabs(curData-P0[0])>realPrecision)
			{	
				type[index] = 0;
				P0[0] = curData;	
				unpredictable_data[unpredictable_count ++] = curData;
			}			
		}
		else
		{
			if(fabs(curData - mean_res) <= realPrecision){
				type[index] = 1;
				P0[0] = mean_res;
			}
			else{
				type[index] = 0;
				P0[0] = curData;
				unpredictable_data[unpredictable_count ++] = curData;
			}
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			curData = cur_data_pos[j];
			pred2D = P0[j-1] + P1[j] - P1[j-1];
			diff = curData - pred2D;
			itvNum = fabs(diff)/realPrecision + 1;
			mean_res = mean[0];
			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type_ = (int) (itvNum/2) + intvRadius;
				pred_res = pred2D + 2 * (type_ - intvRadius) * realPrecision;
				if(fabs(curData - mean_res) <= fabs(curData - pred_res)){
					type[index] = 1;
					P0[j] = mean_res;
				}
				else {
					type[index] = type_;
					P0[j] = pred_res;
				}
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[j])>realPrecision)
				{	
					type[index] = 0;
					P0[j] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
			else
			{
				if(fabs(curData - mean_res) <= realPrecision){
					type[index] = 1;
					P0[j] = mean_res;
				}
				else{
					type[index] = 0;
					P0[j] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}
		
		}

		cur_data_pos += dim1_offset;
	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			curData = cur_data_pos[0];
			index2D = i*r3;
			pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
			diff = curData - pred2D;
			itvNum = fabs(diff)/realPrecision + 1;
			mean_res = mean[0];
			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type_ = (int) (itvNum/2) + intvRadius;
				pred_res = pred2D + 2 * (type_ - intvRadius) * realPrecision;
				if(fabs(curData - mean_res) <= fabs(curData - pred_res)){
					type[index] = 1;
					P0[index2D] = mean_res;
				}
				else {
					type[index] = type_;
					P0[index2D] = pred_res;
				}
				//ganrantee comporession error against the case of machine-epsilon
				if(fabs(curData-P0[index2D])>realPrecision)
				{	
					type[index] = 0;
					P0[index2D] = curData;	
					unpredictable_data[unpredictable_count ++] = curData;
				}				
			}
			else
			{
				if(fabs(curData - mean_res) <= realPrecision){
					type[index] = 1;
					P0[index2D] = mean_res;
				}
				else{
					type[index] = 0;
					P0[index2D] = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				curData = cur_data_pos[j];
				pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
				diff = curData - pred3D;
				itvNum = fabs(diff)/realPrecision + 1;
				mean_res = mean[0];
				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type_ = (int) (itvNum/2) + intvRadius;
					pred_res = pred3D + 2 * (type_ - intvRadius) * realPrecision;
					if(fabs(curData - mean_res) <= fabs(curData - pred_res)){
						type[index] = 1;
						P0[index2D] = mean_res;
					}
					else {
						type[index] = type_;
						P0[index2D] = pred_res;
					}
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-P0[index2D])>realPrecision)
					{	
						type[index] = 0;
						P0[index2D] = curData;	
						unpredictable_data[unpredictable_count ++] = curData;
					}					
				}
				else
				{
					if(fabs(curData - mean_res) <= realPrecision){
						type[index] = 1;
						P0[index2D] = mean_res;
					}
					else{
						type[index] = 0;
						P0[index2D] = curData;
						unpredictable_data[unpredictable_count ++] = curData;
					}
				}
			}
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	return unpredictable_count;
}

unsigned char * SZ_compress_float_1D_MDQ_RA(float *oriData, size_t r1, double realPrecision, size_t * comp_size){
	SZ_Reset(allNodes, stateNum);	
	unsigned int quantization_intervals;
	float dense_pos;
	if(optQuantMode==1)
	{
		// quantization_intervals = optimize_intervals_float_1D(oriData, r1, realPrecision);
		quantization_intervals = optimize_intervals_and_compute_dense_position_float_1D(oriData, r1, realPrecision, &dense_pos);
		printf("number of bins: %d\nerror bound %.4f dense position %.4f\n", quantization_intervals, realPrecision, dense_pos);
		// quantization_intervals = 65536;
		updateQuantizationInfo(quantization_intervals);
		intvCapacity = quantization_intervals - 2;
	}	
	else{
		quantization_intervals = intvCapacity;
		intvCapacity = quantization_intervals - 2;
	}
	size_t num_x;
	size_t early_blockcount_x, late_blockcount_x;
	size_t split_index_x;

	size_t block_size = 512;
	COMPUTE_1D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);


	size_t num_elements = r1;
	size_t max_num_block_elements = early_blockcount_x;
	size_t num_blocks = num_x;

	int * result_type = (int *) malloc(max_num_block_elements * num_blocks * sizeof(int));
	size_t unpred_data_max_size = max_num_block_elements;
	// float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);
	DynamicFloatArray * result_unpredictable_data;
	new_DFA(&result_unpredictable_data, num_elements * 0.01);

	float * data_pos = oriData;
	unsigned short * unpredictable_count = (unsigned short *) malloc(num_blocks * sizeof(int));
	float * mean = malloc(num_blocks * sizeof(float));
	int * type;
	float * unpredictable_data;
	size_t total_unpred = 0;
	size_t index = 0;
	size_t max_unpred_count = 0;
	// printf("Block wise compression start: num_blocks %d num_block_elements %d\n", num_blocks, early_blockcount_x);
	fflush(stdout);
	size_t offset_x = 0;
	size_t type_offset = 0;
	size_t current_blockcount_x;
	for(size_t i=0; i<num_blocks; i++){

		offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
		data_pos = oriData + offset_x;
		type_offset = offset_x;

		type = result_type + type_offset;
		// unpredictable_data = result_unpredictable_data + i * unpred_data_max_size;

		current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
		unpredictable_count[i] = SZ_compress_float_1D_MDQ_RA_block_1D_pred(data_pos, mean + i, dense_pos, num_elements, current_blockcount_x, realPrecision, type, result_unpredictable_data);
		if(unpredictable_count[i] > max_unpred_count){
			max_unpred_count = unpredictable_count[i];
		}
		total_unpred += unpredictable_count[i];
		// printf("done, mean %.2f, unpredictable_count %d, 1st unpredictable_data: %.2f\n", mean, unpredictable_count, unpredictable_data[0]);
		// fflush(stdout);
	}
	// printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	// fflush(stdout);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;
	size_t enCodeSize = 0;

	intToBytes_bigEndian(result_pos, block_size);
	result_pos += 4;
	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;

	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	index = 0;
	size_t unpredictableEncodeSize;
	size_t totalEncodeSize = 0;
	unsigned short * block_pos = (unsigned short *) result_pos;
	unsigned char * block_start_pos = NULL;
	result_pos += num_blocks * sizeof(unsigned short);	// skip block size
	memcpy(result_pos, unpredictable_count, num_blocks * sizeof(unsigned short));
	result_pos += num_blocks * sizeof(unsigned short);	// skip unpredictable count
	memcpy(result_pos, mean, num_blocks * sizeof(float));
	result_pos += num_blocks * sizeof(float);			// skip mean
	// printf("Block wise encode start\n");
	// printf("compress offset to start: %ld\n", result_pos - result);
	// fflush(stdout);
	float * unpredictable_data_pos = result_unpredictable_data->array;
	for(int i=0; i<num_blocks; i++){
		// printf("i j k: %d %d %d\n", i, j, k);
		block_start_pos = result_pos;

		offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
		type_offset = offset_x;
		current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;

		if(unpredictable_count[index] > 0){
			unpredictableEncodeSize = unpredictable_count[index] * sizeof(float);
			// for(int j=0; j<unpredictable_count[index]; j++){
			// 	floatToBytes(result_pos, unpredictable_data_pos[i]);
			// 	result_pos += 4;
			// 	unpredictable_data_pos ++;
			// }
			memcpy(result_pos, unpredictable_data_pos, unpredictableEncodeSize);
			unpredictable_data_pos += unpredictable_count[index];
			result_pos += unpredictableEncodeSize;
		}
		type = result_type + type_offset;
		// caculate real block elements
		enCodeSize = 0;
		encode(type, current_blockcount_x, result_pos, &enCodeSize);

		result_pos += enCodeSize;
		*block_pos = result_pos - block_start_pos;
		block_pos ++;
	}
	totalEncodeSize = result_pos - result;
	printf("Total size %ld\n", totalEncodeSize);
	fflush(stdout);
	free(mean);
	free(unpredictable_count);
	free(result_type);
	SZ_ReleaseHuffman();
	free_DFA(result_unpredictable_data);
	*comp_size = totalEncodeSize;
	return result;
}

unsigned char * SZ_compress_float_3D_MDQ_nonblocked_multi_means(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size){
	unsigned int quantization_intervals;
	float dense_pos;
	float * means;
	unsigned int mean_count;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_and_compute_mean_intervals_float_3D(oriData, r1, r2, r3, realPrecision, &dense_pos, &mean_count, &means);
		printf("number of bins: %d\nerror bound %.20f dense position %.20f\n", quantization_intervals, realPrecision, dense_pos);
		// quantization_intervals = optimize_intervals_float_3D(oriData, r1, realPrecision);
		// printf("new number of bins: %d\n", quantization_intervals);
		// //dense_pos = realPrecision;
		// //dense_pos = 2.5867;
		// if(quantization_intervals < 64) quantization_intervals = 64;
		// printf("adjusted number of bins: %d\n", quantization_intervals);
		updateQuantizationInfo(quantization_intervals);
		intvCapacity = quantization_intervals - 2*((mean_count + 1)/2);
		intvRadius = intvCapacity/2 + 2*((mean_count + 1)/2);
	}	
	else{
		quantization_intervals = intvCapacity;
		intvCapacity = quantization_intervals - 2*((mean_count + 1)/2);
		intvRadius = intvCapacity/2 + 2*((mean_count + 1)/2);
	}
	printf("capacity %d radius %d\n", intvCapacity, intvRadius);

	// calculate block dims
	size_t num_blocks = 1;
	size_t num_elements = r1 * r2 * r3;

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	
	float *P0, *P1; // buffer
	size_t buffer_size = r2 * r3 * sizeof(float);
	P0 = (float *) malloc(buffer_size);
	P1 = (float *) malloc(buffer_size);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	// int unpred_data_max_size = ((int)(num_block_elements * 0.2) + 1) ;
	unsigned int unpred_data_max_size = num_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);

	// int unpredictable_count = 0;
	size_t total_unpred = 0;
	size_t index = 0;
	// NOTE: Currently max unpred count cannot exceed unsigned int limit
	unsigned int max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	float * unpredictable_data = result_unpredictable_data;
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	// float mean;
	// unsigned int unpredictable_count = SZ_compress_float_3D_MDQ_RA_block_3D_pred(data_pos, &mean, dense_pos, r1, r2, r3, r1, r2, r3, realPrecision, P0, P1, type, unpredictable_data);
	unsigned int unpredictable_count = SZ_compress_float_3D_MDQ_RA_block_3D_pred_multi_means(data_pos, mean_count, means, dense_pos, r1, r2, r3, r1, r2, r3, realPrecision, P0, P1, type, unpredictable_data);

	if(unpredictable_count > max_unpred_count){
		max_unpred_count = unpredictable_count;
	}
	total_unpred += unpredictable_count;
	
	//debug
	size_t flushed_count = 0;
	for(size_t i = 0;i<num_elements;i++){
		if(result_type[i] >= 1 && result_type[i] <= mean_count){
			flushed_count ++;
		}
	}
	printf("flushed count: %d\n", flushed_count);
	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	// fflush(stdout);
	free(P0);
	free(P1);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	memcpy(result_pos, &dense_pos, sizeof(float));
	result_pos += 4;
	intToBytes_bigEndian(result_pos, mean_count);
	result_pos += 4;
	memcpy(result_pos, means, mean_count * sizeof(float));
	result_pos += mean_count * sizeof(float);

	size_t totalEncodeSize = 0;
	memcpy(result_pos, &unpredictable_count, sizeof(unsigned int));
	result_pos += num_blocks * sizeof(unsigned int);
	size_t unpredictableEncodeSize = unpredictable_count * sizeof(float);
	memcpy(result_pos, unpredictable_data, unpredictableEncodeSize);
	result_pos += unpredictableEncodeSize;
	// memcpy(result_pos, &mean, sizeof(float));
	// result_pos += num_blocks * sizeof(float);
	
	encode(type, num_elements, result_pos, &enCodeSize);
	result_pos += enCodeSize;
	{
		int status;
		writeIntData_inBytes(result_type, num_elements, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/dataimage-generator/type_array.dat", &status);	
	}

	printf("type array size: %ld\n", enCodeSize);
	printf("--------- True bitrate %.4f ------------\n", enCodeSize * 32.0 / (num_elements * sizeof(float)));

	printf("type array size: %ld\n", enCodeSize);
	totalEncodeSize = result_pos - result;
	printf("Total size %ld\n", totalEncodeSize);
	free(result_unpredictable_data);
	free(result_type);
	SZ_ReleaseHuffman();

	*comp_size = totalEncodeSize;
	return result;
}

unsigned char * SZ_compress_float_3D_MDQ_nonblocked_ori(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size){
	unsigned int quantization_intervals;
	float dense_pos;
	if(optQuantMode==1)
	{
		
		quantization_intervals = optimize_intervals_and_compute_dense_position_float_3D(oriData, r1, r2, r3, realPrecision, &dense_pos);
		printf("number of bins: %d\nerror bound %.20f dense position %.20f\n", quantization_intervals, realPrecision, dense_pos);
		// quantization_intervals = optimize_intervals_float_3D(oriData, r1, r2, r3, realPrecision);
		// printf("new number of bins: %d\n", quantization_intervals);
		//dense_pos = realPrecision;
		//dense_pos = 2.5867;
		//quantization_intervals = 512;
		// if(quantization_intervals < 64) quantization_intervals = 64;
		updateQuantizationInfo(quantization_intervals);
	}	
	else{
		quantization_intervals = intvCapacity;
	}

	// calculate block dims
	size_t num_blocks = 1;
	size_t num_elements = r1 * r2 * r3;

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	
	float *P0, *P1; // buffer
	size_t buffer_size = r2 * r3 * sizeof(float);
	P0 = (float *) malloc(buffer_size);
	P1 = (float *) malloc(buffer_size);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	// int unpred_data_max_size = ((int)(num_block_elements * 0.2) + 1) ;
	unsigned int unpred_data_max_size = num_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);

	// int unpredictable_count = 0;
	size_t total_unpred = 0;
	size_t index = 0;
	// NOTE: Currently max unpred count cannot exceed unsigned int limit
	unsigned int max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	float * unpredictable_data = result_unpredictable_data;
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	float mean;
	unsigned int unpredictable_count = SZ_compress_float_3D_MDQ_RA_block(data_pos, &mean, r1, r2, r3, r1, r2, r3, realPrecision, P0, P1, type, unpredictable_data);
	// float *P_ = (float *) malloc(buffer_size);
	// unsigned int unpredictable_count = SZ_compress_float_3D_MDQ_RA_block_2_layers(data_pos, &mean, r1, r2, r3, r1, r2, r3, realPrecision, P0, P1, P_, type, unpredictable_data);
	// free(P_);

	if(unpredictable_count > max_unpred_count){
		max_unpred_count = unpredictable_count;
	}
	total_unpred += unpredictable_count;
	
	//debug
	size_t flushed_count = 0;
	for(size_t i = 0;i<num_elements;i++){
		if(result_type[i] == 1){
			flushed_count ++;
		}
	}
	printf("flushed count: %d\n", flushed_count);
	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	// fflush(stdout);
	free(P0);
	free(P1);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	size_t totalEncodeSize = 0;
	memcpy(result_pos, &unpredictable_count, sizeof(unsigned int));
	result_pos += num_blocks * sizeof(unsigned int);
	size_t unpredictableEncodeSize = unpredictable_count * sizeof(float);
	memcpy(result_pos, unpredictable_data, unpredictableEncodeSize);
	result_pos += unpredictableEncodeSize;
	memcpy(result_pos, &mean, sizeof(float));
	result_pos += num_blocks * sizeof(float);
	
	encode(type, num_elements, result_pos, &enCodeSize);
	result_pos += enCodeSize;
	printf("type array size: %ld\n", enCodeSize);
	printf("--------- True bitrate %.4f ------------\n", enCodeSize * 32.0 / (num_elements * sizeof(float)));
	{
		int status;
		writeIntData_inBytes(result_type, num_elements, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/dataimage-generator/type_array.dat", &status);	
	}

	printf("type array size: %ld\n", enCodeSize);
	totalEncodeSize = result_pos - result;
	printf("Total size %ld\n", totalEncodeSize);
	free(result_unpredictable_data);
	free(result_type);
	SZ_ReleaseHuffman();

	*comp_size = totalEncodeSize;
	return result;
}

unsigned char * SZ_compress_float_3D_MDQ_nonblocked_adaptive(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size){
	unsigned int quantization_intervals;
	float dense_pos;
	if(optQuantMode==1)
	{
		
		quantization_intervals = optimize_intervals_and_compute_dense_position_float_3D(oriData, r1, r2, r3, realPrecision, &dense_pos);
		printf("number of bins: %d\nerror bound %.20f dense position %.20f\n", quantization_intervals, realPrecision, dense_pos);
		// quantization_intervals = optimize_intervals_float_3D(oriData, r1, r2, r3, realPrecision);
		// printf("new number of bins: %d\n", quantization_intervals);
		//dense_pos = realPrecision;
		//dense_pos = 2.5867;
		//quantization_intervals = 512;
		// if(quantization_intervals < 64) quantization_intervals = 64;
		updateQuantizationInfo(quantization_intervals);
	}	
	else{
		quantization_intervals = intvCapacity;
	}

	// calculate block dims
	size_t num_blocks = 1;
	size_t num_elements = r1 * r2 * r3;

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	
	float *P0, *P1; // buffer
	size_t buffer_size = r2 * r3 * sizeof(float);
	P0 = (float *) malloc(buffer_size);
	P1 = (float *) malloc(buffer_size);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	// int unpred_data_max_size = ((int)(num_block_elements * 0.2) + 1) ;
	unsigned int unpred_data_max_size = num_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);

	// int unpredictable_count = 0;
	size_t total_unpred = 0;
	size_t index = 0;
	// NOTE: Currently max unpred count cannot exceed unsigned int limit
	unsigned int max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	float * unpredictable_data = result_unpredictable_data;
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	float mean;
	// unsigned int unpredictable_count = SZ_compress_float_3D_MDQ_RA_block_3D_pred(data_pos, &mean, dense_pos, r1, r2, r3, r1, r2, r3, realPrecision, P0, P1, type, unpredictable_data);
	unsigned int unpredictable_count = SZ_compress_float_3D_MDQ_RA_block_adaptive(data_pos, &mean, r1, r2, r3, r1, r2, r3, realPrecision, P0, P1, type, unpredictable_data);

	if(unpredictable_count > max_unpred_count){
		max_unpred_count = unpredictable_count;
	}
	total_unpred += unpredictable_count;
	
	//debug
	size_t flushed_count = 0;
	for(size_t i = 0;i<num_elements;i++){
		if(result_type[i] == 1){
			flushed_count ++;
		}
	}
	printf("flushed count: %d\n", flushed_count);
	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	// fflush(stdout);
	free(P0);
	free(P1);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	size_t totalEncodeSize = 0;
	memcpy(result_pos, &unpredictable_count, sizeof(unsigned int));
	result_pos += num_blocks * sizeof(unsigned int);
	size_t unpredictableEncodeSize = unpredictable_count * sizeof(float);
	memcpy(result_pos, unpredictable_data, unpredictableEncodeSize);
	result_pos += unpredictableEncodeSize;
	memcpy(result_pos, &mean, sizeof(float));
	result_pos += num_blocks * sizeof(float);
	
	encode(type, num_elements, result_pos, &enCodeSize);
	result_pos += enCodeSize;
	printf("type array size: %ld\n", enCodeSize);
	printf("--------- True bitrate %.4f ------------\n", enCodeSize * 32.0 / (num_elements * sizeof(float)));
	{
		int status;
		writeIntData_inBytes(result_type, num_elements, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/dataimage-generator/type_array.dat", &status);	
	}

	printf("type array size: %ld\n", enCodeSize);
	totalEncodeSize = result_pos - result;
	printf("Total size %ld\n", totalEncodeSize);
	free(result_unpredictable_data);
	free(result_type);
	SZ_ReleaseHuffman();

	*comp_size = totalEncodeSize;
	return result;
}

unsigned char * SZ_compress_float_3D_MDQ_nonblocked(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size){
	unsigned int quantization_intervals;
	float dense_pos;
	if(optQuantMode==1)
	{
		
		quantization_intervals = optimize_intervals_and_compute_dense_position_float_3D(oriData, r1, r2, r3, realPrecision, &dense_pos);
		printf("number of bins: %d\nerror bound %.20f dense position %.20f\n", quantization_intervals, realPrecision, dense_pos);
		// quantization_intervals = optimize_intervals_float_3D(oriData, r1, r2, r3, realPrecision);
		// printf("new number of bins: %d\n", quantization_intervals);
		//dense_pos = realPrecision;
		//dense_pos = 2.5867;
		//quantization_intervals = 512;
		// if(quantization_intervals < 64) quantization_intervals = 64;
		updateQuantizationInfo(quantization_intervals);
		intvCapacity = quantization_intervals - 2;
	}	
	else{
		quantization_intervals = intvCapacity;
		intvCapacity = quantization_intervals - 2;
	}

	// calculate block dims
	size_t num_blocks = 1;
	size_t num_elements = r1 * r2 * r3;

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	
	float *P0, *P1; // buffer
	size_t buffer_size = r2 * r3 * sizeof(float);
	P0 = (float *) malloc(buffer_size);
	P1 = (float *) malloc(buffer_size);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	// int unpred_data_max_size = ((int)(num_block_elements * 0.2) + 1) ;
	unsigned int unpred_data_max_size = num_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);

	// int unpredictable_count = 0;
	size_t total_unpred = 0;
	size_t index = 0;
	// NOTE: Currently max unpred count cannot exceed unsigned int limit
	unsigned int max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	float * unpredictable_data = result_unpredictable_data;
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	float mean;
	unsigned int unpredictable_count = SZ_compress_float_3D_MDQ_RA_block_3D_pred(data_pos, &mean, dense_pos, r1, r2, r3, r1, r2, r3, realPrecision, P0, P1, type, unpredictable_data);
	if(unpredictable_count > max_unpred_count){
		max_unpred_count = unpredictable_count;
	}
	total_unpred += unpredictable_count;
	
	//debug
	size_t flushed_count = 0;
	for(size_t i = 0;i<num_elements;i++){
		if(result_type[i] == 1){
			flushed_count ++;
		}
	}
	printf("flushed count: %d\n", flushed_count);
	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	// fflush(stdout);
	free(P0);
	free(P1);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	size_t totalEncodeSize = 0;
	memcpy(result_pos, &unpredictable_count, sizeof(unsigned int));
	result_pos += num_blocks * sizeof(unsigned int);
	size_t unpredictableEncodeSize = unpredictable_count * sizeof(float);
	memcpy(result_pos, unpredictable_data, unpredictableEncodeSize);
	result_pos += unpredictableEncodeSize;
	memcpy(result_pos, &mean, sizeof(float));
	result_pos += num_blocks * sizeof(float);
	
	encode(type, num_elements, result_pos, &enCodeSize);
	result_pos += enCodeSize;
	printf("type array size: %ld\n", enCodeSize);
	printf("--------- True bitrate %.4f ------------\n", enCodeSize * 32.0 / (num_elements * sizeof(float)));
	{
		int status;
		writeIntData_inBytes(result_type, num_elements, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/dataimage-generator/type_array.dat", &status);	
	}

	printf("type array size: %ld\n", enCodeSize);
	totalEncodeSize = result_pos - result;
	printf("Total size %ld\n", totalEncodeSize);
	free(result_unpredictable_data);
	free(result_type);
	SZ_ReleaseHuffman();

	*comp_size = totalEncodeSize;
	return result;
}

void SZ_blocked_regression(float * block_ori_data, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, float *params){

	float * data_pos = block_ori_data;
	float fx = 0.0;
	float fy = 0.0;
	float fz = 0.0;
	float f = 0;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	// TODO: use sum_x, sum_y to reduce multiplication operations
	float sum_x, sum_y; 
	float curData;
	for(size_t i=0; i<block_dim_0; i++){
		sum_x = 0;
		for(size_t j=0; j<block_dim_1; j++){
			sum_y = 0;
			for(size_t k=0; k<block_dim_2; k++){
				curData = *data_pos;
				// f += curData;
				// fx += curData * i;
				// fy += curData * j;
				// fz += curData * k;
				sum_y += curData;
				fz += curData * k;
				data_pos ++;
			}
			fy += sum_y * j;
			sum_x += sum_y;
			data_pos += dim1_offset - block_dim_2;
		}
		fx += sum_x * i;
		f += sum_x;
		data_pos += dim0_offset - block_dim_1 * dim1_offset;
	}
	float coeff = 1.0 / (block_dim_0 * block_dim_1 * block_dim_2);
	params[0] = (2 * fx / (block_dim_0 - 1) - f) * 6 * coeff / (block_dim_0 + 1);
	params[1] = (2 * fy / (block_dim_1 - 1) - f) * 6 * coeff / (block_dim_1 + 1);
	params[2] = (2 * fz / (block_dim_2 - 1) - f) * 6 * coeff / (block_dim_2 + 1);
	params[3] = f * coeff - ((block_dim_0 - 1) * params[0] / 2 + (block_dim_1 - 1) * params[1] / 2 + (block_dim_2 - 1) * params[2] / 2);
	return;
}

size_t SZ_compress_float_3D_MDQ_pred_by_regression(float * block_ori_data, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * reg_params, int * type, float * unpredictable_data){

	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	float curData;
	float pred;
	double itvNum;
	double diff;
	size_t index = 0;
	size_t unpredictable_count = 0;
	data_pos = block_ori_data;
	for(size_t i=0; i<block_dim_0; i++){
		for(size_t j=0; j<block_dim_1; j++){
			for(size_t k=0; k<block_dim_2; k++){
				curData = *data_pos;
				pred = reg_params[0] * i + reg_params[1] * j + reg_params[2] * k + reg_params[3];
				diff = curData - pred;
				itvNum = fabs(diff)/realPrecision + 1;
				if (itvNum < intvCapacity){
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					pred = pred + 2 * (type[index] - intvRadius) * realPrecision;
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-pred)>realPrecision){	
						type[index] = 0;
						unpredictable_data[unpredictable_count ++] = curData;
					}		
				}
				else{
					type[index] = 0;
					unpredictable_data[unpredictable_count ++] = curData;
				}
				index ++;	
				data_pos ++;
			}
			data_pos += dim1_offset - block_dim_2;
		}
		data_pos += dim0_offset - block_dim_1 * dim1_offset;
	}
	return unpredictable_count;
}

float SZ_compress_float_3D_MDQ_pred_by_regression_with_err(float * block_ori_data, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * reg_params, int * type, unsigned short * unpred_count, float * unpredictable_data){

	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	float curData;
	float pred;
	double itvNum;
	double diff;
	size_t index = 0;
	size_t unpredictable_count = 0;
	float pred_err = 0.0;
	data_pos = block_ori_data;
	for(size_t i=0; i<block_dim_0; i++){
		for(size_t j=0; j<block_dim_1; j++){
			for(size_t k=0; k<block_dim_2; k++){
				curData = *data_pos;
				pred = reg_params[0] * i + reg_params[1] * j + reg_params[2] * k + reg_params[3];
				diff = curData - pred;
				pred_err += diff * diff;
				itvNum = fabs(diff)/realPrecision + 1;
				if (itvNum < intvCapacity){
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					pred = pred + 2 * (type[index] - intvRadius) * realPrecision;
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData-pred)>realPrecision){	
						type[index] = 0;
						unpredictable_data[unpredictable_count ++] = curData;
					}		
				}
				else{
					type[index] = 0;
					unpredictable_data[unpredictable_count ++] = curData;
				}
				index ++;	
				data_pos ++;
			}
			data_pos += dim1_offset - block_dim_2;
		}
		data_pos += dim0_offset - block_dim_1 * dim1_offset;
	}
	*unpred_count = unpredictable_count;
	return pred_err;
}

float SZ_compress_float_3D_MDQ_strip_pred_by_regression_with_freq(float * block_ori_data, size_t dim_0, size_t dim_1, size_t dim_2, size_t strip_dim_0, size_t strip_dim_1, size_t strip_dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * reg_params, float * pb_pos, int * type, size_t * unpred_count, float * unpredictable_data){

	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	size_t strip_dim0_offset = strip_dim_1 * strip_dim_2;
	size_t strip_dim1_offset = strip_dim_2;
	float curData;
	float pred;
	double itvNum;
	double diff;
	size_t index = 0;
	size_t unpredictable_count = 0;
	size_t correct_pred_count = 0;
	data_pos = block_ori_data;
	for(size_t i=0; i<block_dim_0; i++){
		for(size_t j=0; j<block_dim_1; j++){
			for(size_t k=0; k<block_dim_2; k++){
				curData = *data_pos;
				pred = reg_params[0] * i + reg_params[1] * j + reg_params[2] * k + reg_params[3];
				diff = curData - pred;
				itvNum = fabs(diff)/realPrecision + 1;
				if (itvNum < intvCapacity){
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					pred = pred + 2 * (type[index] - intvRadius) * realPrecision;
					if(type[index] == intvRadius) correct_pred_count ++;
					//ganrantee comporession error against the case of machine-epsilon
					if(fabs(curData - pred)>realPrecision){	
						type[index] = 0;
						pred = curData;
						unpredictable_data[unpredictable_count ++] = curData;
					}		
				}
				else{
					type[index] = 0;
					pred = curData;
					unpredictable_data[unpredictable_count ++] = curData;
				}

				if((i == block_dim_0 - 1) || (j == block_dim_1 - 1) || (k == block_dim_2 - 1)){
					// assign value to strip surfaces
					pb_pos[i * strip_dim0_offset + j * strip_dim1_offset + k] = pred;
				}
				index ++;	
				data_pos ++;
			}
			data_pos += dim1_offset - block_dim_2;
		}
		data_pos += dim0_offset - block_dim_1 * dim1_offset;
	}
	*unpred_count = unpredictable_count;
	return (correct_pred_count * 1.0 / (block_dim_0 * block_dim_1 * block_dim_2));
}

// size_t SZ_compress_float_3D_MDQ_nonblocked_pred_with_blocked_regression_strip(float * block_ori_data, size_t dim_0, size_t dim_1, size_t dim_2, size_t strip_dim_0, size_t strip_dim_1, size_t strip_dim_2,
// 	size_t current_blockcount_x, size_t current_blockcount_y, size_t num_z, size_t split_index_z, size_t early_blockcount_z, size_t late_blockcount_z, double realPrecision, float sz_sample_correct_freq, 
// 	float * prediction_buffer, float * reg_params, int * type, float * unpredictable_data, size_t * block_reg_count, unsigned char * indicator){

// 	// prediction buffer is (current_block_count_x + 1) * (current_block_count_y + 1) * (current_block_count_z + 1)
// 	size_t strip_dim0_offset = strip_dim_1 * strip_dim_2;
// 	size_t strip_dim1_offset = strip_dim_2;
// 	size_t dim0_offset = dim_1 * dim_2;
// 	size_t dim1_offset = dim_2;
// 	float * data_pos = block_ori_data;
// 	size_t current_blockcount_z;
// 	float * reg_params_pos = reg_params;
// 	float * pb_pos = prediction_buffer + strip_dim0_offset + strip_dim1_offset + 1;
// 	size_t strip_unpredictable_count = 0;
// 	size_t unpredictable_count;
// 	size_t reg_count = 0;
// 	float reg_correct_freq;
// 	for(size_t k=0; k<num_z; k++){
// 		current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
// 		SZ_blocked_regression(data_pos, dim_0, dim_1, dim_2, current_blockcount_x, current_blockcount_y, current_blockcount_z, reg_params_pos);
// 		reg_correct_freq = SZ_compress_float_3D_MDQ_strip_pred_by_regression_with_freq(data_pos, dim_0, dim_1, dim_2, strip_dim_0, strip_dim_1, strip_dim_2, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, reg_params_pos, pb_pos, type, &unpredictable_count, unpredictable_data);
// 		// reg_correct_freq = 0;
// 		if(reg_correct_freq > sz_sample_correct_freq){
// 			// use regression, which is default in indicator
// 			strip_unpredictable_count += unpredictable_count;
// 			unpredictable_data += unpredictable_count;
// 			reg_params_pos += 4;
// 			reg_count ++;
// 		}
// 		else{
// 			// use SZ
// 			{
// 				// SZ predication
// 				unpredictable_count = 0;
// 				float * cur_pb_pos = pb_pos;
// 				float * cur_data_pos = data_pos;
// 				float curData;
// 				float pred3D;
// 				double itvNum, diff;
// 				size_t index = 0;
// 				for(size_t ii=0; ii<current_blockcount_x; ii++){
// 					for(size_t jj=0; jj<current_blockcount_y; jj++){
// 						for(size_t kk=0; kk<current_blockcount_z; kk++){
// 							if(ii==0 && jj==0 && kk==0){
// 								pred3D = 2;
// 							}
// 							pred3D = cur_pb_pos[-1] + cur_pb_pos[-strip_dim1_offset]+ cur_pb_pos[-strip_dim0_offset] - cur_pb_pos[-strip_dim1_offset - 1]
// 									 - cur_pb_pos[-strip_dim0_offset - 1] - cur_pb_pos[-strip_dim0_offset - strip_dim1_offset] + cur_pb_pos[-strip_dim0_offset - strip_dim1_offset - 1];
// 							curData = cur_data_pos[0];
// 							diff = curData - pred3D;
// 							itvNum = fabs(diff)/realPrecision + 1;
// 							if (itvNum < intvCapacity){
// 								if (diff < 0) itvNum = -itvNum;
// 								type[index] = (int) (itvNum/2) + intvRadius;
// 								cur_pb_pos[0] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
// 								//ganrantee comporession error against the case of machine-epsilon
// 								if(fabs(curData - cur_pb_pos[0])>realPrecision){	
// 									type[index] = 0;
// 									cur_pb_pos[0] = curData;	
// 									unpredictable_data[unpredictable_count ++] = curData;
// 								}					
// 							}
// 							else{
// 								type[index] = 0;
// 								cur_pb_pos[0] = curData;
// 								unpredictable_data[unpredictable_count ++] = curData;
// 							}

// 							index ++;
// 							cur_pb_pos ++;
// 							cur_data_pos ++;
// 						}
// 						cur_pb_pos += strip_dim1_offset - current_blockcount_z;
// 						cur_data_pos += dim1_offset - current_blockcount_z;
// 					}
// 					cur_pb_pos += strip_dim0_offset - current_blockcount_y * strip_dim1_offset;
// 					cur_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
// 				}
// 			}
// 			strip_unpredictable_count += unpredictable_count;
// 			unpredictable_data += unpredictable_count;
// 			// change indicator
// 			indicator[k] = 1;
// 		}// end SZ
// 		data_pos += current_blockcount_z;
// 		pb_pos += current_blockcount_z;
// 		type += current_blockcount_x * current_blockcount_y * current_blockcount_z;
// 	}
// 	*block_reg_count = reg_count;
// 	return strip_unpredictable_count;
// }

// void SZ_compress_regression_coefficients(float * coeff, size_t length, double realPrecision, unsigned char * bytes, size_t * bytes_size, DynamicByteArray * lead_num, DynamicByteArray * mid_byte, DynamicIntArray * residue_bit){
// 	float range;
// 	float median;
// 	computeRangeSize_float(coeff, length, &range, &median);
// 	short radius = getExponent_float(range/2);
// 	int reqLength;
// 	computeReqLength_float(realPrecision, radius, &reqLength, &median);
// 	int reqBytesLength = reqLength / 8;
// 	int resiBitsLength = reqLength % 8;
// 	FloatValueCompressElement * vce = (FloatValueCompressElement *) malloc(sizeof(FloatValueCompressElement));
// 	LossyCompressionElement * lce = (LossyCompressionElement *) malloc(sizeof(LossyCompressionElement));
// 	unsigned char preDataBytes[4] = {0, 0, 0, 0};
// 	// CHANGE REALPRECISION
// 	for(size_t i=0; i<length; i++){
// 		compressSingleFloatValue(vce, coeff[i], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
// 		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
// 		memcpy(preDataBytes, vce->curBytes, 4);
// 		addExactData(lead_num, mid_byte, residue_bit, lce);
// 	}

// 	free(FloatValueCompressElement);
// 	free(LossyCompressionElement);
// }

// void SZ_concat_coefficient_array(){

// }

unsigned char * SZ_compress_float_3D_MDQ_nonblocked_with_blocked_regression(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size){

	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_3D(oriData, r1, r2, r3, realPrecision);
		// quantization_intervals = optimize_intervals_and_compute_dense_position_float_3D(oriData, r1, r2, r3, realPrecision, &dense_pos);
		printf("number of bins: %d\n", quantization_intervals);
		updateQuantizationInfo(quantization_intervals);
	}	
	else{
		quantization_intervals = intvCapacity;
	}
	// calculate block dims
	size_t num_x, num_y, num_z;
	size_t block_size = 8;
	printf("block_size: %ld\n", block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	COLL_BASE_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;
	size_t num_elements = r1 * r2 * r3;

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;	

	int * result_type = (int *) malloc(num_elements * sizeof(int));
	size_t unpred_data_max_size = max_num_block_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);
	size_t total_unpred = 0;
	size_t unpredictable_count;
	size_t index = 0;
	size_t max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	size_t type_offset;
	// float * unpredictable_data = result_unpredictable_data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	float * reg_params = (float *) malloc(num_blocks * 4 * sizeof(float));
	float * reg_params_pos = reg_params;
	float * unpredictable_data = result_unpredictable_data;
	unsigned char * indicator = (unsigned char *) malloc(num_blocks * sizeof(unsigned char));
	memset(indicator, 0, num_blocks * sizeof(unsigned char));
	size_t current_block_elements;
	size_t reg_count = 0;
	size_t strip_dim_0 = (early_blockcount_x + 1);
	size_t strip_dim_1 = (early_blockcount_y + 1);
	size_t strip_dim_2 = r3 + 1;
	size_t strip_dim0_offset = strip_dim_1 * strip_dim_2;
	size_t strip_dim1_offset = strip_dim_2;
	size_t prediction_buffer_size = strip_dim_0 * strip_dim_1 * strip_dim_2 * sizeof(float);
	float * prediction_buffer = (float *) malloc(prediction_buffer_size);
	memset(prediction_buffer, 0, prediction_buffer_size);
	float * buffer1 = (float *) malloc((r2 + 1) * (r3 + 1) * sizeof(float));
	memset(buffer1, 0, (r2+1) * (r3+1) * sizeof(float));
	float * buffer2 = (float *) malloc((r2 + 1) * (r3 + 1) * sizeof(float));
	float * bottom_buffer = buffer1;
	float * bottom_buffer2 = buffer2;
	unsigned char * indicator_pos = indicator;
	float sz_sample_correct_freq = 0.5;
	float reg_correct_freq;
	// move regression part out
	size_t params_offset_b = num_blocks;
	size_t params_offset_c = 2*num_blocks;
	size_t params_offset_d = 3*num_blocks;
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
	
				data_pos = oriData + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;
				// SZ_blocked_regression(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, reg_params_pos);
				{
					float * cur_data_pos = data_pos;
					float fx = 0.0;
					float fy = 0.0;
					float fz = 0.0;
					float f = 0;
					float sum_x, sum_y; 
					float curData;
					for(size_t i=0; i<current_blockcount_x; i++){
						sum_x = 0;
						for(size_t j=0; j<current_blockcount_y; j++){
							sum_y = 0;
							for(size_t k=0; k<current_blockcount_z; k++){
								curData = *cur_data_pos;
								// f += curData;
								// fx += curData * i;
								// fy += curData * j;
								// fz += curData * k;
								sum_y += curData;
								fz += curData * k;
								cur_data_pos ++;
							}
							fy += sum_y * j;
							sum_x += sum_y;
							cur_data_pos += dim1_offset - current_blockcount_z;
						}
						fx += sum_x * i;
						f += sum_x;
						cur_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
					}
					float coeff = 1.0 / (current_blockcount_x * current_blockcount_y * current_blockcount_z);
					reg_params_pos[0] = (2 * fx / (current_blockcount_x - 1) - f) * 6 * coeff / (current_blockcount_x + 1);
					reg_params_pos[params_offset_b] = (2 * fy / (current_blockcount_y - 1) - f) * 6 * coeff / (current_blockcount_y + 1);
					reg_params_pos[params_offset_c] = (2 * fz / (current_blockcount_z - 1) - f) * 6 * coeff / (current_blockcount_z + 1);
					reg_params_pos[params_offset_d] = f * coeff - ((current_blockcount_x - 1) * reg_params_pos[0] / 2 + (current_blockcount_y - 1) * reg_params_pos[params_offset_b] / 2 + (current_blockcount_z - 1) * reg_params_pos[params_offset_c] / 2);
				}
				reg_params_pos ++;
			}
		}
	}

	//Compress coefficient arrays
	float medianValue_a, medianValue_b, medianValue_c, medianValue_d;
	double precision_a, precision_b, precision_c, precision_d;
	float rel_param_err = 0.0025;
	precision_a = rel_param_err * realPrecision / late_blockcount_x;
	precision_b = rel_param_err * realPrecision / late_blockcount_y;
	precision_c = rel_param_err * realPrecision / late_blockcount_z;
	precision_d = rel_param_err * realPrecision;

	unsigned char* leadArray_a = NULL, *midArray_a = NULL, *resiArray_a = NULL;
	int reqBytesLength_a = 0, resiBitsLength_a = 0; 
	compressExactDataArray(reg_params, precision_a, num_blocks, &leadArray_a, &midArray_a, &resiArray_a, &reqBytesLength_a, &resiBitsLength_a, &medianValue_a);	
	unsigned char* leadArray_b = NULL, *midArray_b = NULL, *resiArray_b = NULL;
	int reqBytesLength_b = 0, resiBitsLength_b = 0; 
	compressExactDataArray(reg_params+params_offset_b, precision_b, num_blocks, &leadArray_b, &midArray_b, &resiArray_b, &reqBytesLength_b, &resiBitsLength_b, &medianValue_b);	
	unsigned char* leadArray_c = NULL, *midArray_c = NULL, *resiArray_c = NULL;
	int reqBytesLength_c = 0, resiBitsLength_c = 0; 
	compressExactDataArray(reg_params+params_offset_c, precision_c, num_blocks, &leadArray_c, &midArray_c, &resiArray_c, &reqBytesLength_c, &resiBitsLength_c, &medianValue_c);			
	unsigned char* leadArray_d = NULL, *midArray_d = NULL, *resiArray_d = NULL;
	int reqBytesLength_d = 0, resiBitsLength_d = 0; 
	compressExactDataArray(reg_params+params_offset_d, precision_d, num_blocks, &leadArray_d, &midArray_d, &resiArray_d, &reqBytesLength_d, &resiBitsLength_d, &medianValue_d);	

	unsigned char *int_lead_a = (unsigned char*)malloc(num_blocks);
	unsigned char* pre_mid_a_p = midArray_a, *byte_mid_a = (unsigned char*)malloc(num_blocks*4);
	unsigned char* post_mid_a_p = byte_mid_a;
	unsigned char* int_resi_a = (unsigned char*)malloc(num_blocks); 
	
	unsigned char* int_lead_b = (unsigned char*)malloc(num_blocks);
	unsigned char* pre_mid_b_p = midArray_b, *byte_mid_b = (unsigned char*)malloc(num_blocks*4);
	unsigned char* post_mid_b_p = byte_mid_b;
	unsigned char* int_resi_b = (unsigned char*)malloc(num_blocks); 
	
	unsigned char* int_lead_c = (unsigned char*)malloc(num_blocks);
	unsigned char* pre_mid_c_p = midArray_c, *byte_mid_c = (unsigned char*)malloc(num_blocks*4);
	unsigned char* post_mid_c_p = byte_mid_c;
	unsigned char* int_resi_c = (unsigned char*)malloc(num_blocks); 
	
	unsigned char* int_lead_d = (unsigned char*)malloc(num_blocks);
	unsigned char* pre_mid_d_p = midArray_d, *byte_mid_d = (unsigned char*)malloc(num_blocks*4);
	unsigned char* post_mid_d_p = byte_mid_d;
	unsigned char* int_resi_d = (unsigned char*)malloc(num_blocks); 		 			
	
	size_t sel_block_index = 0;
	reg_params_pos = reg_params;	
		
	double tmp_realPrecision = realPrecision;
	for(size_t i=0; i<num_x; i++){
		current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
		offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
		// memset front surface
		for(size_t ii = 1; ii < current_blockcount_x + 1; ii ++){
			memset(prediction_buffer + ii * strip_dim0_offset, 0, strip_dim1_offset * sizeof(float));
		}
		for(size_t j=0; j<num_y; j++){
			offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
			current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
			data_pos = oriData + offset_x * dim0_offset + offset_y * dim1_offset;
			// copy bottom plane from plane buffer
			memcpy(prediction_buffer, bottom_buffer + offset_y * strip_dim1_offset, (current_blockcount_y + 1) * strip_dim1_offset * sizeof(float));
			type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset;
			type = result_type + type_offset;

			// prediction buffer is (current_block_count_x + 1) * (current_block_count_y + 1) * (current_block_count_z + 1)
			size_t current_blockcount_z;
			float * pb_pos = prediction_buffer + strip_dim0_offset + strip_dim1_offset + 1;
			size_t strip_unpredictable_count = 0;
			for(size_t k=0; k<num_z; k++){
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
				// SZ_blocked_regression(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, reg_params_pos);
				// reg_correct_freq = SZ_compress_float_3D_MDQ_strip_pred_by_regression_with_freq(data_pos, r1, r2, r3, strip_dim_0, strip_dim_1, strip_dim_2, current_blockcount_x, current_blockcount_y, current_blockcount_z, tmp_realPrecision, reg_params_pos, pb_pos, type, &unpredictable_count, unpredictable_data);
				{
					float curData;
					float pred;
					double itvNum;
					double diff;
					size_t index = 0;
					size_t block_unpredictable_count = 0;
					size_t correct_pred_count = 0;
					float * cur_data_pos = data_pos;
					for(size_t ii=0; ii<current_blockcount_x; ii++){
						for(size_t jj=0; jj<current_blockcount_y; jj++){
							for(size_t kk=0; kk<current_blockcount_z; kk++){
								curData = *cur_data_pos;
								pred = reg_params_pos[0] * ii + reg_params_pos[params_offset_b] * jj + reg_params_pos[params_offset_c] * kk + reg_params_pos[params_offset_d];
								diff = curData - pred;
								itvNum = fabs(diff)/tmp_realPrecision + 1;
								if (itvNum < intvCapacity){
									if (diff < 0) itvNum = -itvNum;
									type[index] = (int) (itvNum/2) + intvRadius;
									pred = pred + 2 * (type[index] - intvRadius) * tmp_realPrecision;
									if(type[index] == intvRadius) correct_pred_count ++;
									//ganrantee comporession error against the case of machine-epsilon
									if(fabs(curData - pred)>tmp_realPrecision){	
										type[index] = 0;
										pred = curData;
										unpredictable_data[block_unpredictable_count ++] = curData;
									}		
								}
								else{
									type[index] = 0;
									pred = curData;
									unpredictable_data[block_unpredictable_count ++] = curData;
								}

								if((ii == current_blockcount_x - 1) || (jj == current_blockcount_y - 1) || (kk == current_blockcount_z - 1)){
									// assign value to strip surfaces
									pb_pos[ii * strip_dim0_offset + jj * strip_dim1_offset + kk] = pred;
								}
								index ++;	
								cur_data_pos ++;
							}
							cur_data_pos += dim1_offset - current_blockcount_z;
						}
						cur_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
					}
					unpredictable_count = block_unpredictable_count;
					reg_correct_freq = (correct_pred_count * 1.0 / (current_blockcount_x * current_blockcount_y * current_blockcount_z));				
				}
				reg_params_pos ++;
				//compute the offset for mid_array
				int midByteSize_a = reqBytesLength_a - leadArray_a[sel_block_index];
				int midByteSize_b = reqBytesLength_b - leadArray_b[sel_block_index];
				int midByteSize_c = reqBytesLength_c - leadArray_c[sel_block_index];
				int midByteSize_d = reqBytesLength_d - leadArray_d[sel_block_index];
				
				if(reg_correct_freq > sz_sample_correct_freq){
					// use regression, which is default in indicator
					strip_unpredictable_count += unpredictable_count;
					unpredictable_data += unpredictable_count;
					
					int_lead_a[reg_count] = leadArray_a[sel_block_index];
					int_lead_b[reg_count] = leadArray_b[sel_block_index];					
					int_lead_c[reg_count] = leadArray_c[sel_block_index];
					int_lead_d[reg_count] = leadArray_d[sel_block_index];										
					
					int_resi_a[reg_count] = resiArray_a[sel_block_index];
					int_resi_b[reg_count] = resiArray_b[sel_block_index];
					int_resi_c[reg_count] = resiArray_c[sel_block_index];
					int_resi_d[reg_count] = resiArray_d[sel_block_index];															
					
					memcpy(post_mid_a_p, pre_mid_a_p, midByteSize_a);
					memcpy(post_mid_b_p, pre_mid_b_p, midByteSize_b);
					memcpy(post_mid_c_p, pre_mid_c_p, midByteSize_c);
					memcpy(post_mid_d_p, pre_mid_d_p, midByteSize_d);
					
					post_mid_a_p += midByteSize_a;
					post_mid_b_p += midByteSize_b;
					post_mid_c_p += midByteSize_c;
					post_mid_d_p += midByteSize_d;
					
					// reg_params_pos += 4;
					reg_count ++;
				}
				else{
					// use SZ
					// SZ predication
					unpredictable_count = 0;
					float * cur_pb_pos = pb_pos;
					float * cur_data_pos = data_pos;
					float curData;
					float pred3D;
					double itvNum, diff;
					size_t index = 0;
					for(size_t ii=0; ii<current_blockcount_x; ii++){
						for(size_t jj=0; jj<current_blockcount_y; jj++){
							for(size_t kk=0; kk<current_blockcount_z; kk++){

								pred3D = cur_pb_pos[-1] + cur_pb_pos[-strip_dim1_offset]+ cur_pb_pos[-strip_dim0_offset] - cur_pb_pos[-strip_dim1_offset - 1]
										 - cur_pb_pos[-strip_dim0_offset - 1] - cur_pb_pos[-strip_dim0_offset - strip_dim1_offset] + cur_pb_pos[-strip_dim0_offset - strip_dim1_offset - 1];
								curData = cur_data_pos[0];
								diff = curData - pred3D;
								itvNum = fabs(diff)/realPrecision + 1;
								if (itvNum < intvCapacity){
									if (diff < 0) itvNum = -itvNum;
									type[index] = (int) (itvNum/2) + intvRadius;
									cur_pb_pos[0] = pred3D + 2 * (type[index] - intvRadius) * tmp_realPrecision;
									//ganrantee comporession error against the case of machine-epsilon
									if(fabs(curData - cur_pb_pos[0])>tmp_realPrecision){	
										type[index] = 0;
										cur_pb_pos[0] = curData;	
										unpredictable_data[unpredictable_count ++] = curData;
									}					
								}
								else{
									type[index] = 0;
									cur_pb_pos[0] = curData;
									unpredictable_data[unpredictable_count ++] = curData;
								}

								index ++;
								cur_pb_pos ++;
								cur_data_pos ++;
							}
							cur_pb_pos += strip_dim1_offset - current_blockcount_z;
							cur_data_pos += dim1_offset - current_blockcount_z;
						}
						cur_pb_pos += strip_dim0_offset - current_blockcount_y * strip_dim1_offset;
						cur_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
					}
					strip_unpredictable_count += unpredictable_count;
					unpredictable_data += unpredictable_count;
					// change indicator
					indicator_pos[k] = 1;
				}// end SZ
				//TODO
				pre_mid_a_p += midByteSize_a;
				pre_mid_b_p += midByteSize_b;
				pre_mid_c_p += midByteSize_c;
				pre_mid_d_p += midByteSize_d;
				
				sel_block_index++;
				data_pos += current_blockcount_z;
				pb_pos += current_blockcount_z;
				type += current_blockcount_x * current_blockcount_y * current_blockcount_z;

			}
			// copy upper plane to bottom plane buffer
			memcpy(bottom_buffer2 + (offset_y + 1) * strip_dim1_offset, prediction_buffer + current_blockcount_x * strip_dim0_offset + strip_dim1_offset, current_blockcount_y * strip_dim1_offset * sizeof(float));
			// copy back surface to front surface
			for(size_t ii=1; ii<current_blockcount_x + 1; ii++){
				memcpy(prediction_buffer + 1 + ii * strip_dim0_offset, prediction_buffer + 1 + ii * strip_dim0_offset + current_blockcount_y * strip_dim1_offset, r3 * sizeof(float));
			}
			if(strip_unpredictable_count > max_unpred_count){
				max_unpred_count = strip_unpredictable_count;
			}
			total_unpred += strip_unpredictable_count;
			indicator_pos += num_z;
			// if(i == 16 && j == 25 && k == 34) {
			// 	printf("done, mean %.2f, unpredictable_count %d, 1st unpredictable_data: %.2f\n", mean[index], unpredictable_count[index], unpredictable_data[0]);
			// 	exit(0);
			// }
		}
		float * tmp;
		tmp = bottom_buffer;
		bottom_buffer = bottom_buffer2;
		bottom_buffer2 = tmp;
	}
	free(prediction_buffer);
	free(buffer1);
	free(buffer2);
	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	printf("sz_count: %ld reg_count: %ld\n", num_blocks - reg_count, reg_count);
	//reg_count = num_x * num_y * num_z;
	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) calloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int), 1);
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	intToBytes_bigEndian(result_pos, block_size);
	result_pos += 4;
	doubleToBytes(result_pos, realPrecision);
	// printf("%.17f\n", realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	memcpy(result_pos, indicator, num_blocks * sizeof(unsigned char));
	result_pos += num_blocks * sizeof(unsigned char);
	// {
	// 	size_t tmp_count = 0;
	// 	for(size_t i=0; i<num_blocks; i++){
	// 		tmp_count += indicator[i];
	// 	}
	// 	printf("total reg_count: %ld\n", tmp_count);
	// 	exit(0);
	// }
	//memcpy(result_pos, reg_params, reg_count * 4 * sizeof(float));
	//result_pos += reg_count * 4 * sizeof(float);
	// reorder reg_count
	// float * params_pos = (float *) result_pos;
	// result_pos += reg_count * 4 * sizeof(float);	
	// for(size_t i=0; i<4; i++){
	// 	for(size_t j=0; j<reg_count; j++){
	// 		params_pos[i*reg_count + j] = reg_params[j*4 + i];
	// 	}
	// }
	
	//convert the lead/mid/resi to byte stream 
	
	
	/*int sum_mid = 0, sum_lead = 0;
	for(int i = 0;i<num_blocks;i++)
	{
		if(int_lead_a[i] > 3 || int_lead_b[i] > 3 || int_lead_c[i] > 3 || int_lead_d[i] > 3)
		{
			printf("%d: %d\n", i, int_lead_a[i]);
			exit(0);
		}
		sum_mid += (reqBytesLength_d - int_lead_a[i]);
	}*/		
	
	
	if(reg_count>0)
	{
		// write coefficient info
		memcpy(result_pos, &medianValue_a, sizeof(float));
		result_pos += sizeof(float);
		memcpy(result_pos, &medianValue_b, sizeof(float));
		result_pos += sizeof(float);
		memcpy(result_pos, &medianValue_c, sizeof(float));
		result_pos += sizeof(float);
		memcpy(result_pos, &medianValue_d, sizeof(float));
		result_pos += sizeof(float);

		int reqLength_a, reqLength_b, reqLength_c, reqLength_d;
		reqLength_a = reqBytesLength_a * 8 + resiBitsLength_a;
		memcpy(result_pos, &reqLength_a, sizeof(int));
		result_pos += sizeof(int);
		reqLength_b = reqBytesLength_b * 8 + resiBitsLength_b;
		memcpy(result_pos, &reqLength_b, sizeof(int));
		result_pos += sizeof(int);
		reqLength_c = reqBytesLength_c * 8 + resiBitsLength_c;
		memcpy(result_pos, &reqLength_c, sizeof(int));
		result_pos += sizeof(int);
		reqLength_d = reqBytesLength_d * 8 + resiBitsLength_d;
		memcpy(result_pos, &reqLength_d, sizeof(int));
		result_pos += sizeof(int);

		size_t leadNumArray_size_a = convertIntArray2ByteArray_fast_2b_inplace(int_lead_a, reg_count, result_pos);
		result_pos += leadNumArray_size_a;
		size_t leadNumArray_size_b = convertIntArray2ByteArray_fast_2b_inplace(int_lead_b, reg_count, result_pos);
		result_pos += leadNumArray_size_b;
		size_t leadNumArray_size_c = convertIntArray2ByteArray_fast_2b_inplace(int_lead_c, reg_count, result_pos);
		result_pos += leadNumArray_size_c;
		size_t leadNumArray_size_d = convertIntArray2ByteArray_fast_2b_inplace(int_lead_d, reg_count, result_pos);
		result_pos += leadNumArray_size_d;
		
		size_t mid_byte_size_a = post_mid_a_p - byte_mid_a;
		memcpy(result_pos, &mid_byte_size_a, sizeof(size_t));
		result_pos += sizeof(size_t);
		memcpy(result_pos, byte_mid_a, mid_byte_size_a);
		result_pos += mid_byte_size_a;

		size_t mid_byte_size_b = post_mid_b_p - byte_mid_b;
		memcpy(result_pos, &mid_byte_size_b, sizeof(size_t));
		result_pos += sizeof(size_t);
		memcpy(result_pos, byte_mid_b, mid_byte_size_b);
		result_pos += mid_byte_size_b;
		
		size_t mid_byte_size_c = post_mid_c_p - byte_mid_c;
		memcpy(result_pos, &mid_byte_size_c, sizeof(size_t));
		result_pos += sizeof(size_t);
		memcpy(result_pos, byte_mid_c, mid_byte_size_c);
		result_pos += mid_byte_size_c;
		
		size_t mid_byte_size_d = post_mid_d_p - byte_mid_d;
		memcpy(result_pos, &mid_byte_size_d, sizeof(size_t));
		result_pos += sizeof(size_t);
		memcpy(result_pos, byte_mid_d, mid_byte_size_d);
		result_pos += mid_byte_size_d;		
		
		unsigned char* byteBuffer = NULL;
		size_t residualResiBits_size_a = convertIntArray2ByteArray_fast_dynamic(int_resi_a, resiBitsLength_a, reg_count, &byteBuffer);
		memcpy(result_pos, &residualResiBits_size_a, sizeof(size_t));
		result_pos += sizeof(size_t);
		memcpy(result_pos, byteBuffer, residualResiBits_size_a);
		result_pos += residualResiBits_size_a;
		free(byteBuffer);

		size_t residualResiBits_size_b = convertIntArray2ByteArray_fast_dynamic(int_resi_b, resiBitsLength_b, reg_count, &byteBuffer);
		memcpy(result_pos, &residualResiBits_size_b, sizeof(size_t));
		result_pos += sizeof(size_t);
		memcpy(result_pos, byteBuffer, residualResiBits_size_b);
		result_pos += residualResiBits_size_b;
		free(byteBuffer);
		
		size_t residualResiBits_size_c = convertIntArray2ByteArray_fast_dynamic(int_resi_c, resiBitsLength_c, reg_count, &byteBuffer);
		memcpy(result_pos, &residualResiBits_size_c, sizeof(size_t));
		result_pos += sizeof(size_t);
		memcpy(result_pos, byteBuffer, residualResiBits_size_c);
		result_pos += residualResiBits_size_c;
		free(byteBuffer);
		
		size_t residualResiBits_size_d = convertIntArray2ByteArray_fast_dynamic(int_resi_c, resiBitsLength_c, reg_count, &byteBuffer);
		memcpy(result_pos, &residualResiBits_size_d, sizeof(size_t));
		result_pos += sizeof(size_t);
		memcpy(result_pos, byteBuffer, residualResiBits_size_d);
		result_pos += residualResiBits_size_d;
		free(byteBuffer);		
	}
	
	//record the number of unpredictable data and also store them
	memcpy(result_pos, &total_unpred, sizeof(size_t));
	result_pos += sizeof(size_t);
	memcpy(result_pos, result_unpredictable_data, total_unpred * sizeof(float));
	result_pos += total_unpred * sizeof(float);
	// int status;
	// writeFloatData_inBytes(mean, num_blocks, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/mean12.dat", &status);
	size_t typeArray_size = 0;
	encode(result_type, num_elements, result_pos, &typeArray_size);
	result_pos += typeArray_size;
	printf("type array size: %ld\n", typeArray_size);
	printf("--------- True bitrate %.4f ------------\n", typeArray_size * 32.0 / (num_elements * sizeof(float)));
	size_t totalEncodeSize = result_pos - result;
	// printf("Total size %ld\n", totalEncodeSize);
	free(indicator);
	free(result_unpredictable_data);
	free(result_type);
	free(reg_params);

	free(leadArray_a);
	free(midArray_a);
	free(resiArray_a);
	free(leadArray_b);
	free(midArray_b);
	free(resiArray_b);
	free(leadArray_c);
	free(midArray_c);
	free(resiArray_c);
	free(leadArray_d);
	free(midArray_d);
	free(resiArray_d);
	free(int_lead_a);
	free(byte_mid_a);
	free(int_resi_a);
	free(int_lead_b);
	free(byte_mid_b);
	free(int_resi_b);
	free(int_lead_c);
	free(byte_mid_c);
	free(int_resi_c);
	free(int_lead_d);
	free(byte_mid_d);
	free(int_resi_d);

	SZ_ReleaseHuffman();
	*comp_size = totalEncodeSize;
	return result;
}

unsigned char * SZ_compress_float_3D_MDQ_RA_all_by_regression(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size){

	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_3D(oriData, r1, r2, r3, realPrecision);
		// quantization_intervals = optimize_intervals_and_compute_dense_position_float_3D(oriData, r1, r2, r3, realPrecision, &dense_pos);
		// printf("3D number of bins: %d\nerror bound %.20f dense position %.20f\n", quantization_intervals, realPrecision, dense_pos);
		// exit(0);		
		// if(quantization_intervals < 256) quantization_intervals = 256;
		printf("number of bins: %d\n", quantization_intervals);
		// //dense_pos = realPrecision;
		// //dense_pos = 2.5867;
		// printf("adjusted number of bins: %d\n", quantization_intervals);
		updateQuantizationInfo(quantization_intervals);
	}	
	else{
		quantization_intervals = intvCapacity;
	}

	// calculate block dims
	size_t num_x, num_y, num_z;
	size_t block_size = 8;
	printf("block_size: %ld\n", block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	COLL_BASE_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;
	size_t num_elements = r1 * r2 * r3;

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;

	float *P0, *P1; // buffer
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	size_t unpred_data_max_size = max_num_block_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);
	unsigned short * unpredictable_count = (unsigned short *) malloc(num_blocks * sizeof(int));
	size_t total_unpred = 0;
	size_t index = 0;
	size_t max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	// float * unpredictable_data = result_unpredictable_data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t type_offset = 0;
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	float * reg_params = (float *) malloc(num_blocks * 4 * sizeof(float));
	float * unpredictable_data_pos = result_unpredictable_data;
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
				data_pos = oriData + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
				// printf("x y z: %d %d %d\n", current_blockcount_x, current_blockcount_y, current_blockcount_z);
				// fflush(stdout);
				// if(current_blockcount_x != 8 || current_blockcount_y != 8 || current_blockcount_z != 8)
				// 	exit(0);
				type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
				type = result_type + type_offset;
				// printf("i j k: %d %d %d, offset %ld %ld %ld type offset %ld\n", i, j, k, offset_x, offset_y, offset_z, type_offset);

				index = i * num_y * num_z + j * num_z + k;
				SZ_blocked_regression(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, reg_params + index * 4);

				// unpredictable_data = result_unpredictable_data + index * unpred_data_max_size;
				unpredictable_count[index] = SZ_compress_float_3D_MDQ_pred_by_regression(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, reg_params + index * 4, type, unpredictable_data_pos);
				if(unpredictable_count[index] > max_unpred_count){
					max_unpred_count = unpredictable_count[index];
				}
				unpredictable_data_pos += unpredictable_count[index];
				total_unpred += unpredictable_count[index];
				// if(i == 16 && j == 25 && k == 34) {
				// 	printf("done, mean %.2f, unpredictable_count %d, 1st unpredictable_data: %.2f\n", mean[index], unpredictable_count[index], unpredictable_data[0]);
				// 	exit(0);
				// }
			}
		}
	}
	
	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	// fflush(stdout);
	// size_t typeArray_size;
	// unsigned char * typeArray;// = (unsigned char *) malloc(num_elements * sizeof(int));
	// encode_withTree(result_type, num_elements, &typeArray, &typeArray_size);
	// free(typeArray);
	// printf("typeArray_size: %ld\n", typeArray_size);
	// int status;
	// writeIntData_inBytes(result_type, num_elements, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/uf01_bs12.dat", &status);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	intToBytes_bigEndian(result_pos, block_size);
	result_pos += 4;
	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	size_t unpredictableEncodeSize;
	size_t totalEncodeSize = 0;
	unsigned char * block_start_pos = NULL;
	memcpy(result_pos, unpredictable_count, num_blocks * sizeof(unsigned short));
	result_pos += num_blocks * sizeof(unsigned short);
	memcpy(result_pos, result_unpredictable_data, total_unpred * sizeof(float));
	result_pos += total_unpred * sizeof(float);
	memcpy(result_pos, reg_params, num_blocks * 4 * sizeof(float));
	result_pos += num_blocks * 4 * sizeof(float);
	// int status;
	// writeFloatData_inBytes(mean, num_blocks, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/mean12.dat", &status);
	size_t typeArray_size = 0;
	encode(result_type, num_elements, result_pos, &typeArray_size);
	result_pos += typeArray_size;
	printf("type array size: %ld\n", typeArray_size);
	printf("--------- True bitrate %.4f ------------\n", typeArray_size * 32.0 / (num_elements * sizeof(float)));
	totalEncodeSize = result_pos - result;
	// printf("Total size %ld\n", totalEncodeSize);
	free(result_unpredictable_data);
	free(unpredictable_count);
	free(result_type);
	free(reg_params);
	SZ_ReleaseHuffman();
	*comp_size = totalEncodeSize;
	return result;
}

unsigned char * SZ_compress_float_3D_MDQ_RA_blocked_with_regression(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size){

	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_3D(oriData, r1, r2, r3, realPrecision);
		// quantization_intervals = optimize_intervals_and_compute_dense_position_float_3D(oriData, r1, r2, r3, realPrecision, &dense_pos);
		// printf("3D number of bins: %d\nerror bound %.20f dense position %.20f\n", quantization_intervals, realPrecision, dense_pos);
		// exit(0);		
		// if(quantization_intervals < 256) quantization_intervals = 256;
		printf("number of bins: %d\n", quantization_intervals);
		// //dense_pos = realPrecision;
		// //dense_pos = 2.5867;
		// printf("adjusted number of bins: %d\n", quantization_intervals);
		updateQuantizationInfo(quantization_intervals);
	}	
	else{
		quantization_intervals = intvCapacity;
	}

	// calculate block dims
	size_t num_x, num_y, num_z;
	size_t block_size = 8;
	printf("block_size: %ld\n", block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	COLL_BASE_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;
	size_t num_elements = r1 * r2 * r3;

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;

	float *P0, *P1; // buffer
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	size_t unpred_data_max_size = max_num_block_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);
	unsigned short * unpredictable_count = (unsigned short *) malloc(num_blocks * sizeof(int));
	size_t total_unpred = 0;
	size_t index = 0;
	size_t max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	// float * unpredictable_data = result_unpredictable_data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t type_offset = 0;
	size_t buffer_size = r2 * r3 * sizeof(float);
	P0 = (float *) malloc(buffer_size);
	P1 = (float *) malloc(buffer_size);
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	float * reg_params = (float *) malloc(num_blocks * 4 * sizeof(float));
	float * unpredictable_data_pos = result_unpredictable_data;
	unsigned char * indicator = (unsigned char *) malloc(num_blocks * sizeof(unsigned char));
	memset(indicator, 0, num_blocks * sizeof(unsigned char));
	unsigned short unpredictable_count_reg = 0;
	unsigned short unpredictable_count_sz = 0;
	float * unpred_buffer_sz = (float *) malloc(max_num_block_elements * sizeof(float));
	int * type_buffer_sz = (int *) malloc(max_num_block_elements * sizeof(int));
	float pred_err_sz, pred_err_reg;
	size_t current_block_elements;
	size_t sz_count = 0;
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
				data_pos = oriData + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
				// printf("x y z: %d %d %d\n", current_blockcount_x, current_blockcount_y, current_blockcount_z);
				// fflush(stdout);
				// if(current_blockcount_x != 8 || current_blockcount_y != 8 || current_blockcount_z != 8)
				// 	exit(0);
				type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
				type = result_type + type_offset;
				// printf("i j k: %d %d %d, offset %ld %ld %ld type offset %ld\n", i, j, k, offset_x, offset_y, offset_z, type_offset);

				index = i * num_y * num_z + j * num_z + k;
				SZ_blocked_regression(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, reg_params + index * 4);

				if(i == 0 && j== 16 && k == 47){
					printf("sz_count: %ld\n", sz_count);
				}
				pred_err_reg = SZ_compress_float_3D_MDQ_pred_by_regression_with_err(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, reg_params + index * 4, type, unpredictable_count + index, unpredictable_data_pos);
				pred_err_sz = SZ_compress_float_3D_MDQ_RA_block_no_mean(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, P0, P1, type_buffer_sz, &unpredictable_count_sz, unpred_buffer_sz);
				if(pred_err_reg > pred_err_sz){
					unpredictable_count[index] = unpredictable_count_sz;
					current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					memcpy(unpredictable_data_pos, unpred_buffer_sz, unpredictable_count_sz * sizeof(float));
					memcpy(type, type_buffer_sz, current_block_elements * sizeof(int));
					indicator[index] = 1;
					sz_count ++;
				}
				if(unpredictable_count[index] > max_unpred_count){
					max_unpred_count = unpredictable_count[index];
				}
				unpredictable_data_pos += unpredictable_count[index];
				total_unpred += unpredictable_count[index];
				// if(i == 16 && j == 25 && k == 34) {
				// 	printf("done, mean %.2f, unpredictable_count %d, 1st unpredictable_data: %.2f\n", mean[index], unpredictable_count[index], unpredictable_data[0]);
				// 	exit(0);
				// }
			}
		}
	}
	
	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	printf("sz_count: %ld reg_count: %ld\n", sz_count, num_blocks - sz_count);
	// fflush(stdout);
	free(unpred_buffer_sz);
	free(type_buffer_sz);
	free(P0);
	free(P1);
	// size_t typeArray_size;
	// unsigned char * typeArray;// = (unsigned char *) malloc(num_elements * sizeof(int));
	// encode_withTree(result_type, num_elements, &typeArray, &typeArray_size);
	// free(typeArray);
	// printf("typeArray_size: %ld\n", typeArray_size);
	// int status;
	// writeIntData_inBytes(result_type, num_elements, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/uf01_bs12.dat", &status);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	intToBytes_bigEndian(result_pos, block_size);
	result_pos += 4;
	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	size_t unpredictableEncodeSize;
	size_t totalEncodeSize = 0;
	unsigned char * block_start_pos = NULL;
	memcpy(result_pos, indicator, num_blocks * sizeof(unsigned char));
	result_pos += num_blocks * sizeof(unsigned char);
	memcpy(result_pos, unpredictable_count, num_blocks * sizeof(unsigned short));
	result_pos += num_blocks * sizeof(unsigned short);
	memcpy(result_pos, result_unpredictable_data, total_unpred * sizeof(float));
	result_pos += total_unpred * sizeof(float);
	memcpy(result_pos, reg_params, num_blocks * 4 * sizeof(float));
	result_pos += num_blocks * 4 * sizeof(float);
	// int status;
	// writeFloatData_inBytes(mean, num_blocks, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/mean12.dat", &status);
	size_t typeArray_size = 0;
	encode(result_type, num_elements, result_pos, &typeArray_size);
	result_pos += typeArray_size;
	printf("type array size: %ld\n", typeArray_size);
	printf("--------- True bitrate %.4f ------------\n", typeArray_size * 32.0 / (num_elements * sizeof(float)));
	totalEncodeSize = result_pos - result;
	// printf("Total size %ld\n", totalEncodeSize);
	free(indicator);
	free(result_unpredictable_data);
	free(unpredictable_count);
	free(result_type);
	free(reg_params);
	SZ_ReleaseHuffman();
	*comp_size = totalEncodeSize;
	return result;
}


unsigned char * SZ_compress_float_3D_MDQ_RA(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size){

	unsigned int quantization_intervals;
	float dense_pos;
	if(optQuantMode==1)
	{
		// quantization_intervals = optimize_intervals_float_3D(oriData, r1, realPrecision);
		quantization_intervals = optimize_intervals_and_compute_dense_position_float_3D(oriData, r1, r2, r3, realPrecision, &dense_pos);
		printf("3D number of bins: %d\nerror bound %.20f dense position %.20f\n", quantization_intervals, realPrecision, dense_pos);
		// exit(0);		
		// if(quantization_intervals < 256) quantization_intervals = 256;
		// printf("new number of bins: %d\n", quantization_intervals);
		// //dense_pos = realPrecision;
		// //dense_pos = 2.5867;
		// printf("adjusted number of bins: %d\n", quantization_intervals);
		updateQuantizationInfo(quantization_intervals);
		intvCapacity = quantization_intervals - 2;
	}	
	else{
		quantization_intervals = intvCapacity;
		intvCapacity = quantization_intervals - 2;
	}

	// calculate block dims
	size_t num_x, num_y, num_z;
	size_t block_size = 12;
	printf("block_size: %ld\n", block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	COLL_BASE_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	// printf("nums: %d %d %d\n", num_x, num_y, num_z);
	// printf("splits: %d %d %d\n", split_index_x, split_index_y, split_index_z);
	// printf("early counts: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// printf("late counts: %d %d %d\n", late_blockcount_x, late_blockcount_y, late_blockcount_z);
	// exit(0);
	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;
	size_t num_elements = r1 * r2 * r3;
	// printf("max_num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	
	// printf("malloc blockinfo array start\n");
	// fflush(stdout);

	float *P0, *P1; // buffer
	size_t buffer_size = r2 * r3 * sizeof(float);
	P0 = (float *) malloc(buffer_size);
	P1 = (float *) malloc(buffer_size);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	// int unpred_data_max_size = ((int)(num_block_elements * 0.2) + 1) ;
	size_t unpred_data_max_size = max_num_block_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);

	// int unpredictable_count = 0;
	unsigned short * unpredictable_count = (unsigned short *) malloc(num_blocks * sizeof(int));
	float * mean = malloc(num_blocks * sizeof(float));
	size_t total_unpred = 0;
	size_t index = 0;
	size_t max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	float * unpredictable_data = result_unpredictable_data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t type_offset = 0;
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
				data_pos = oriData + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
				// printf("x y z: %d %d %d\n", current_blockcount_x, current_blockcount_y, current_blockcount_z);
				// fflush(stdout);
				// if(current_blockcount_x != 8 || current_blockcount_y != 8 || current_blockcount_z != 8)
				// 	exit(0);
				type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
				type = result_type + type_offset;
				// printf("i j k: %d %d %d, offset %ld %ld %ld type offset %ld\n", i, j, k, offset_x, offset_y, offset_z, type_offset);

				index = i * num_y * num_z + j * num_z + k;
				unpredictable_data = result_unpredictable_data + index * unpred_data_max_size;
				unpredictable_count[index] = SZ_compress_float_3D_MDQ_RA_block(data_pos, mean + index, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, P0, P1, type, unpredictable_data);
				// unpredictable_count[index] = SZ_compress_float_3D_MDQ_RA_block_3D_pred(data_pos, mean + index, dense_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, P0, P1, type, unpredictable_data);
				// unpredictable_count[index] = SZ_compress_float_3D_MDQ_RA_block_1D_pred(data_pos, mean + index, dense_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, type, unpredictable_data);
				// unpredictable_count[index] = SZ_compress_float_3D_MDQ_RA_block_3D_pred_flush_after_compare(data_pos, mean + index, dense_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, P0, P1, type, unpredictable_data);
				if(unpredictable_count[index] > max_unpred_count){
					max_unpred_count = unpredictable_count[index];
				}
				total_unpred += unpredictable_count[index];
				// if(i == 16 && j == 25 && k == 34) {
				// 	printf("done, mean %.2f, unpredictable_count %d, 1st unpredictable_data: %.2f\n", mean[index], unpredictable_count[index], unpredictable_data[0]);
				// 	exit(0);
				// }
			}
		}
	}
	
	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	// fflush(stdout);
	free(P0);
	free(P1);
	// size_t typeArray_size;
	// unsigned char * typeArray;// = (unsigned char *) malloc(num_elements * sizeof(int));
	// encode_withTree(result_type, num_elements, &typeArray, &typeArray_size);
	// free(typeArray);
	// printf("typeArray_size: %ld\n", typeArray_size);
	// int status;
	// writeIntData_inBytes(result_type, num_elements, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/uf01_bs12.dat", &status);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	intToBytes_bigEndian(result_pos, block_size);
	result_pos += 4;
	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	size_t unpredictableEncodeSize;
	size_t totalEncodeSize = 0;
	unsigned short * block_pos = (unsigned short *) result_pos;
	unsigned short * print_bs = block_pos;
	unsigned char * block_start_pos = NULL;
	result_pos += num_blocks * sizeof(unsigned short); // skip block size
	memcpy(result_pos, unpredictable_count, num_blocks * sizeof(unsigned short));
	result_pos += num_blocks * sizeof(unsigned short);
	memcpy(result_pos, mean, num_blocks * sizeof(float));
	result_pos += num_blocks * sizeof(float);
	
	size_t typeArray_size = 0;
	
	//debug
	size_t flushed_count = 0;
	for(size_t i = 0;i<num_elements;i++){
		if(result_type[i] == 1){
			flushed_count ++;
		}
	}
	printf("flushed count: %d\n", flushed_count);
	// int status;
	// writeFloatData_inBytes(mean, num_blocks, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/mean12.dat", &status);

	size_t current_block_elements;
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				// printf("i j k: %d %d %d\n", i, j, k);
				index = i * num_y * num_z + j * num_z + k;
				block_start_pos = result_pos;
				// memcpy(result_pos, mean+index, 4);
				// result_pos += 4;

				if(unpredictable_count[index] > 0){
					unpredictable_data = result_unpredictable_data + index * unpred_data_max_size;
					unpredictableEncodeSize = unpredictable_count[index] * sizeof(float);
					memcpy(result_pos, unpredictable_data, unpredictableEncodeSize);
					result_pos += unpredictableEncodeSize;
				}
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
				type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;

				type = result_type + type_offset;
				current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
				enCodeSize = 0;
				// printf("%ld %ld %ld\n", i, j, k);
				// fflush(stdout);
				// if( i == 0 && j == 1 && k == 30){
				// 	for(int i=0; i<current_block_elements; i++){
				// 		printf("%d ", type[i]);
				// 	}
				// 	printf("\n");
				// 	exit(0);
				// }
				encode(type, current_block_elements, result_pos, &enCodeSize);
				typeArray_size += enCodeSize;
				result_pos += enCodeSize;
				*block_pos = result_pos - block_start_pos;
				block_pos ++;
			}
		}
	}
	printf("type array size: %ld\n", typeArray_size);
	printf("--------- True bitrate %.4f ------------\n", typeArray_size * 32.0 / (num_elements * sizeof(float)));
	totalEncodeSize = result_pos - result;
	// printf("Total size %ld\n", totalEncodeSize);
	free(mean);
	free(result_unpredictable_data);
	free(unpredictable_count);
	free(result_type);
	SZ_ReleaseHuffman();

	// writeUShortData_inBytes(print_bs, num_blocks, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/uf01_bs12.dat", &status);
	*comp_size = totalEncodeSize;
	return result;
}

unsigned char * SZ_compress_float_3D_MDQ_RA_multi_means(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size){

	unsigned int quantization_intervals;
	float dense_pos;
	float * means;
	unsigned int mean_count;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_and_compute_mean_intervals_float_3D(oriData, r1, r2, r3, realPrecision, &dense_pos, &mean_count, &means);
		printf("number of bins: %d\nerror bound %.20f dense position %.20f\n", quantization_intervals, realPrecision, dense_pos);
		// quantization_intervals = optimize_intervals_float_3D(oriData, r1, realPrecision);
		// printf("new number of bins: %d\n", quantization_intervals);
		// //dense_pos = realPrecision;
		// //dense_pos = 2.5867;
		// if(quantization_intervals < 64) quantization_intervals = 64;
		// printf("adjusted number of bins: %d\n", quantization_intervals);
		updateQuantizationInfo(quantization_intervals);
		if(quantization_intervals < 256) quantization_intervals = 256;
		intvCapacity = quantization_intervals - 2*((mean_count + 1)/2);
		intvRadius = intvCapacity/2 + 2*((mean_count + 1)/2);
	}	
	else{
		quantization_intervals = intvCapacity;
		intvCapacity = quantization_intervals - 2*((mean_count + 1)/2);
		intvRadius = intvCapacity/2 + 2*((mean_count + 1)/2);
	}
	printf("capacity %d radius %d\n", intvCapacity, intvRadius);

	// calculate block dims
	size_t num_x, num_y, num_z;
	size_t block_size = 12;
	COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	COLL_BASE_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;
	size_t num_elements = r1 * r2 * r3;

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	
	float *P0, *P1; // buffer
	size_t buffer_size = r2 * r3 * sizeof(float);
	P0 = (float *) malloc(buffer_size);
	P1 = (float *) malloc(buffer_size);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	size_t unpred_data_max_size = max_num_block_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);

	unsigned short * unpredictable_count = (unsigned short *) malloc(num_blocks * sizeof(int));
	size_t total_unpred = 0;
	size_t index = 0;
	size_t max_unpred_count = 0;
	float * data_pos = oriData;
	int * type = result_type;
	float * unpredictable_data = result_unpredictable_data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t type_offset = 0;
	// printf("Block wise compression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
				data_pos = oriData + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
				// printf("x y z: %d %d %d\n", current_blockcount_x, current_blockcount_y, current_blockcount_z);
				// fflush(stdout);
				// if(current_blockcount_x != 8 || current_blockcount_y != 8 || current_blockcount_z != 8)
				// 	exit(0);
				type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
				type = result_type + type_offset;
				// printf("i j k: %d %d %d, offset %ld %ld %ld type offset %ld\n", i, j, k, offset_x, offset_y, offset_z, type_offset);
				// if(i == 0 && j == 3 && k == 54){

				// }
				index = i * num_y * num_z + j * num_z + k;
				unpredictable_data = result_unpredictable_data + index * unpred_data_max_size;
				unpredictable_count[index] = SZ_compress_float_3D_MDQ_RA_block_3D_pred_multi_means(data_pos, mean_count, means, dense_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, P0, P1, type, unpredictable_data);
				if(unpredictable_count[index] > max_unpred_count){
					max_unpred_count = unpredictable_count[index];
				}
				total_unpred += unpredictable_count[index];
				// if(i == 16 && j == 25 && k == 34) {
				// 	printf("done, mean %.2f, unpredictable_count %d, 1st unpredictable_data: %.2f\n", mean[index], unpredictable_count[index], unpredictable_data[0]);
				// 	exit(0);
				// }
			}
		}
	}
	
	printf("Block wise compression end, unpredictable num %d, num_elements %ld, max unpred count %d\n", total_unpred, num_elements, max_unpred_count);
	// fflush(stdout);
	free(P0);
	free(P1);

	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	intToBytes_bigEndian(result_pos, block_size);
	result_pos += 4;
	doubleToBytes(result_pos, realPrecision);
	result_pos += 8;
	intToBytes_bigEndian(result_pos, quantization_intervals);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, treeByteSize);
	result_pos += 4;
	intToBytes_bigEndian(result_pos, nodeCount);
	result_pos += 4;
	memcpy(result_pos, treeBytes, treeByteSize);
	result_pos += treeByteSize;
	free(treeBytes);

	size_t unpredictableEncodeSize;
	size_t totalEncodeSize = 0;
	memcpy(result_pos, &dense_pos, sizeof(float));
	result_pos += 4;
	intToBytes_bigEndian(result_pos, mean_count);
	result_pos += 4;
	memcpy(result_pos, means, mean_count * sizeof(float));
	result_pos += mean_count * sizeof(float);
	unsigned short * block_pos = (unsigned short *) result_pos;
	unsigned char * block_start_pos = NULL;
	result_pos += num_blocks * sizeof(unsigned short); // skip block size
	memcpy(result_pos, unpredictable_count, num_blocks * sizeof(unsigned short));
	result_pos += num_blocks * sizeof(unsigned short);
	
	size_t typeArray_size = 0;
	
	//debug
	size_t flushed_count = 0;
	for(size_t i = 0;i<num_elements;i++){
		if(result_type[i] > 0 && result_type[i] <= mean_count){
			flushed_count ++;
		}
	}
	printf("flushed count: %d\n", flushed_count);
		
	size_t current_block_elements;
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				// printf("i j k: %d %d %d\n", i, j, k);
				index = i * num_y * num_z + j * num_z + k;
				block_start_pos = result_pos;
				// memcpy(result_pos, mean+index, 4);
				// result_pos += 4;

				if(unpredictable_count[index] > 0){
					unpredictable_data = result_unpredictable_data + index * unpred_data_max_size;
					unpredictableEncodeSize = unpredictable_count[index] * sizeof(float);
					memcpy(result_pos, unpredictable_data, unpredictableEncodeSize);
					result_pos += unpredictableEncodeSize;
				}
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
				type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;

				type = result_type + type_offset;
				current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
				enCodeSize = 0;
				// printf("%ld %ld %ld\n", i, j, k);
				// fflush(stdout);
				// if( i == 0 && j == 1 && k == 30){
				// 	for(int i=0; i<current_block_elements; i++){
				// 		printf("%d ", type[i]);
				// 	}
				// 	printf("\n");
				// 	exit(0);
				// }
				encode(type, current_block_elements, result_pos, &enCodeSize);
				typeArray_size += enCodeSize;
				result_pos += enCodeSize;
				*block_pos = result_pos - block_start_pos;
				block_pos ++;
			}
		}
	}
	printf("type array size: %ld\n", typeArray_size);
	totalEncodeSize = result_pos - result;
	// printf("Total size %ld\n", totalEncodeSize);
	free(means);
	free(result_unpredictable_data);
	free(unpredictable_count);
	free(result_type);
	SZ_ReleaseHuffman();

	*comp_size = totalEncodeSize;
	return result;
}



