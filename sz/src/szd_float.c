/**
 *  @file szd_float.c
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
 *  @brief 
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "szd_float.h"
#include "TightDataPointStorageF.h"
#include "sz.h"
#include "Huffman.h"
#include "szd_float_pwr.h"
//#include "rw.h"

/**
 * 
 * 
 * @return status SUCCESSFUL (SZ_SCES) or not (other error codes) f
 * */
int SZ_decompress_args_float(float** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize)
{
	int status = SZ_SCES;
	size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
	
	//unsigned char* tmpBytes;
	size_t targetUncompressSize = dataLength <<2; //i.e., *4
	//tmpSize must be "much" smaller than dataLength
	size_t i, tmpSize = 8+MetaDataByteLength+SZ_SIZE_TYPE;
	unsigned char* szTmpBytes;	
	
	if(cmpSize!=8+4+MetaDataByteLength && cmpSize!=8+8+MetaDataByteLength) //4,8 means two posibilities of SZ_SIZE_TYPE
	{
		int isZlib = isZlibFormat(cmpBytes[0], cmpBytes[1]);
		if(isZlib)
			szMode = SZ_BEST_COMPRESSION;
		else
			szMode = SZ_BEST_SPEED;		
		if(szMode==SZ_BEST_SPEED)
		{
			tmpSize = cmpSize;
			szTmpBytes = cmpBytes;	
		}
		else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
		{
			if(targetUncompressSize<MIN_ZLIB_DEC_ALLOMEM_BYTES) //Considering the minimum size
				targetUncompressSize = MIN_ZLIB_DEC_ALLOMEM_BYTES; 
			tmpSize = zlib_uncompress5(cmpBytes, (unsigned long)cmpSize, &szTmpBytes, (unsigned long)targetUncompressSize+4+MetaDataByteLength+SZ_SIZE_TYPE);//		(unsigned long)targetUncompressSize+8: consider the total length under lossless compression mode is actually 3+4+1+targetUncompressSize
			//szTmpBytes = (unsigned char*)malloc(sizeof(unsigned char)*tmpSize);
			//memcpy(szTmpBytes, tmpBytes, tmpSize);
			//free(tmpBytes); //release useless memory		
		}
		else
		{
			printf("Wrong value of szMode in the double compressed bytes.\n");
			status = SZ_MERR;
			return status;
		}	
	}
	else
		szTmpBytes = cmpBytes;
	//TODO: convert szTmpBytes to data array.
	TightDataPointStorageF* tdps = NULL;
	//int isConstant = 0, isRandomAccess = 0; 
	//isSameOrRandomAccessMode(szTmpBytes, &isConstant, &isRandomAccess);
	
	int errBoundMode = new_TightDataPointStorageF_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
	
	//writeByteData(tdps->typeArray, tdps->typeArray_size, "decompress-typebytes.tbt");
	int dim = computeDimension(r5,r4,r3,r2,r1);	
	int floatSize = sizeof(float);
	if(tdps->isLossless)
	{
		*newData = (float*)malloc(floatSize*dataLength);
		if(sysEndianType==BIG_ENDIAN_SYSTEM)
		{
			memcpy(*newData, szTmpBytes+4+MetaDataByteLength+SZ_SIZE_TYPE, dataLength*floatSize);
		}
		else
		{
			unsigned char* p = szTmpBytes+4+MetaDataByteLength+SZ_SIZE_TYPE;
			for(i=0;i<dataLength;i++,p+=floatSize)
				(*newData)[i] = bytesToFloat(p);
		}		
	}
	else if (tdps->allSameData) {
		float value = bytesToFloat(tdps->exactMidBytes);
		*newData = (float*)malloc(sizeof(float)*dataLength);
		for (i = 0; i < dataLength; i++)
			(*newData)[i] = value;
	}
	else
	{
		if(tdps->raBytes_size > 0) //random access mode
		{
			if (dim == 1)
				decompressDataSeries_float_1D_RA(newData, r1, tdps->raBytes);
			else if(dim == 2)
				decompressDataSeries_float_2D_nonblocked_with_blocked_regression(newData, r2, r1, tdps->raBytes+SZ_SIZE_TYPE);
			else if(dim == 3)
				//decompressDataSeries_float_3D_nonblocked(newData, r3, r2, r1, tdps->raBytes);
				decompressDataSeries_float_3D_nonblocked_with_blocked_regression(newData, r3, r2, r1, tdps->raBytes+SZ_SIZE_TYPE);
				// decompressDataSeries_float_3D_all_by_interpolation(newData, r3, r2, r1, tdps->raBytes+SZ_SIZE_TYPE);
			else if(dim == 4)
				;
			else
			{
				printf("Error: currently support only at most 4 dimensions!\n");
				status = SZ_DERR;
			}	
		}
		else
		{
			if (dim == 1)
				getSnapshotData_float_1D(newData,r1,tdps, errBoundMode);
			else if (dim == 2)
				getSnapshotData_float_2D(newData,r2,r1,tdps, errBoundMode);
			else if (dim == 3)
				getSnapshotData_float_3D(newData,r3,r2,r1,tdps, errBoundMode);
			else if (dim == 4)
				getSnapshotData_float_4D(newData,r4,r3,r2,r1,tdps, errBoundMode);
			else
			{
				printf("Error: currently support only at most 4 dimensions!\n");
				status = SZ_DERR;
			}			
		}
	}
	

	free_TightDataPointStorageF(tdps);
	if(szMode!=SZ_BEST_SPEED && cmpSize!=8+MetaDataByteLength+SZ_SIZE_TYPE)
		free(szTmpBytes);
	SZ_ReleaseHuffman();	
	return status;
}

void decompressDataSeries_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
								// in resiMidBits, p is to track the
								// byte_index of resiMidBits, l is for
								// leadNum
	unsigned char* leadNum;
	double interval = tdps->realPrecision*2;
	
	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);

	//sdi:Debug
	//writeUShortData(type, dataSeriesLength, "decompressStateBytes.sb");

	unsigned char preBytes[4];
	unsigned char curBytes[4];
	
	memset(preBytes, 0, 4);

	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits; 
	unsigned char leadingNum;	
	float medianValue, exactData, predValue;
	
	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;
	
	int type_;
	for (i = 0; i < dataSeriesLength; i++) {
		type_ = type[i];
		switch (type_) {
		case 0:
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data	
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}
			
			exactData = bytesToFloat(curBytes);
			(*data)[i] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
			break;
		default:
			//predValue = 2 * (*data)[i-1] - (*data)[i-2];
			predValue = (*data)[i-1];
			(*data)[i] = predValue + (type_-intvRadius)*interval;
			break;
		}
		//printf("%.30G\n",(*data)[i]);
	}
	free(leadNum);
	free(type);
	return;
}

void decompressDataSeries_float_2D(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	//printf("tdps->intervals=%d, intvRadius=%d\n", tdps->intervals, intvRadius);
	
	size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	size_t dataSeriesLength = r1*r2;
	//	printf ("%d %d\n", r1, r2);

	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);

	unsigned char preBytes[4];
	unsigned char curBytes[4];

	memset(preBytes, 0, 4);

	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits; 
	unsigned char leadingNum;	
	float medianValue, exactData, predValue;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;
	
	float pred1D, pred2D;
	size_t ii, jj;

	/* Process Row-0, data 0 */

	// compute resiBits
	resiBits = 0;
	if (resiBitsLength != 0) {
		int kMod8 = k % 8;
		int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
		if (rightMovSteps > 0) {
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
		} else if (rightMovSteps < 0) {
			int code1 = getLeftMovingCode(kMod8);
			int code2 = getRightMovingCode(kMod8, resiBitsLength);
			int leftMovSteps = -rightMovSteps;
			rightMovSteps = 8 - leftMovSteps;
			resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
			p++;
			resiBits = resiBits
					| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
		} else // rightMovSteps == 0
		{
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code);
			p++;
		}
		k += resiBitsLength;
	}

	// recover the exact data
	memset(curBytes, 0, 4);
	leadingNum = leadNum[l++];
	memcpy(curBytes, preBytes, leadingNum);
	for (j = leadingNum; j < reqBytesLength; j++)
		curBytes[j] = tdps->exactMidBytes[curByteIndex++];
	if (resiBitsLength != 0) {
		unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
		curBytes[reqBytesLength] = resiByte;
	}

	exactData = bytesToFloat(curBytes);
	(*data)[0] = exactData + medianValue;
	memcpy(preBytes,curBytes,4);

	/* Process Row-0, data 1 */
	type_ = type[1]; 
	if (type_ != 0)
	{
		pred1D = (*data)[0];
		(*data)[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else
	{
		// compute resiBits
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
			} else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data
		memset(curBytes, 0, 4);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}

		exactData = bytesToFloat(curBytes);
		(*data)[1] = exactData + medianValue;
		memcpy(preBytes,curBytes,4);
	}

	/* Process Row-0, data 2 --> data r2-1 */
	for (jj = 2; jj < r2; jj++)
	{
		type_ = type[jj];
		if (type_ != 0)
		{
			pred1D = 2*(*data)[jj-1] - (*data)[jj-2];				
			(*data)[jj] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[jj] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}
	}

	size_t index;
	/* Process Row-1 --> Row-r1-1 */
	for (ii = 1; ii < r1; ii++)
	{
		/* Process row-ii data 0 */
		index = ii*r2;

		type_ = type[index];
		if (type_ != 0)
		{
			pred1D = (*data)[index-r2];		
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process row-ii data 1 --> r2-1*/
		for (jj = 1; jj < r2; jj++)
		{
			index = ii*r2+jj;
			pred2D = (*data)[index-1] + (*data)[index-r2] - (*data)[index-r2-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}
	}

	free(leadNum);
	free(type);
	return;
}

void decompressDataSeries_float_3D(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	size_t dataSeriesLength = r1*r2*r3;
	size_t r23 = r2*r3;
//	printf ("%d %d %d\n", r1, r2, r3);
	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);
	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);

	unsigned char preBytes[4];
	unsigned char curBytes[4];

	memset(preBytes, 0, 4);
	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits;
	unsigned char leadingNum;
	float medianValue, exactData, predValue;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;
	
	float pred1D, pred2D, pred3D;
	size_t ii, jj, kk;

	///////////////////////////	Process layer-0 ///////////////////////////
	/* Process Row-0 data 0*/
	// compute resiBits
	resiBits = 0;
	if (resiBitsLength != 0) {
		int kMod8 = k % 8;
		int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
		if (rightMovSteps > 0) {
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
		} else if (rightMovSteps < 0) {
			int code1 = getLeftMovingCode(kMod8);
			int code2 = getRightMovingCode(kMod8, resiBitsLength);
			int leftMovSteps = -rightMovSteps;
			rightMovSteps = 8 - leftMovSteps;
			resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
			p++;
			resiBits = resiBits
					| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
		} else // rightMovSteps == 0
		{
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code);
			p++;
		}
		k += resiBitsLength;
	}

	// recover the exact data
	memset(curBytes, 0, 4);
	leadingNum = leadNum[l++];
	memcpy(curBytes, preBytes, leadingNum);
	for (j = leadingNum; j < reqBytesLength; j++)
		curBytes[j] = tdps->exactMidBytes[curByteIndex++];
	if (resiBitsLength != 0) {
		unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
		curBytes[reqBytesLength] = resiByte;
	}
	exactData = bytesToFloat(curBytes);
	(*data)[0] = exactData + medianValue;
	memcpy(preBytes,curBytes,4);

	/* Process Row-0, data 1 */
	pred1D = (*data)[0];

	type_ = type[1];
	if (type_ != 0)
	{
		(*data)[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else
	{
		// compute resiBits
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
			} else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data
		memset(curBytes, 0, 4);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}

		exactData = bytesToFloat(curBytes);
		(*data)[1] = exactData + medianValue;
		memcpy(preBytes,curBytes,4);
	}
	/* Process Row-0, data 2 --> data r3-1 */
	for (jj = 2; jj < r3; jj++)
	{
		pred1D = 2*(*data)[jj-1] - (*data)[jj-2];

		type_ = type[jj];
		if (type_ != 0)
		{
			(*data)[jj] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[jj] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}
	}

	size_t index;
	/* Process Row-1 --> Row-r2-1 */
	for (ii = 1; ii < r2; ii++)
	{
		/* Process row-ii data 0 */
		index = ii*r3;
		pred1D = (*data)[index-r3];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process row-ii data 1 --> r3-1*/
		for (jj = 1; jj < r3; jj++)
		{
			index = ii*r3+jj;
			pred2D = (*data)[index-1] + (*data)[index-r3] - (*data)[index-r3-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}
	}

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (kk = 1; kk < r1; kk++)
	{
		/* Process Row-0 data 0*/
		index = kk*r23;
		pred1D = (*data)[index-r23];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process Row-0 data 1 --> data r3-1 */
		for (jj = 1; jj < r3; jj++)
		{
			index = kk*r23+jj;
			pred2D = (*data)[index-1] + (*data)[index-r23] - (*data)[index-r23-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}

		/* Process Row-1 --> Row-r2-1 */
		for (ii = 1; ii < r2; ii++)
		{
			/* Process Row-i data 0 */
			index = kk*r23 + ii*r3;
			pred2D = (*data)[index-r3] + (*data)[index-r23] - (*data)[index-r23-r3];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (jj = 1; jj < r3; jj++)
			{
				index = kk*r23 + ii*r3 + jj;
				pred3D = (*data)[index-1] + (*data)[index-r3] + (*data)[index-r23]
					- (*data)[index-r3-1] - (*data)[index-r23-r3] - (*data)[index-r23-1] + (*data)[index-r23-r3-1];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
				}
				else
				{
					// compute resiBits
					resiBits = 0;
					if (resiBitsLength != 0) {
						int kMod8 = k % 8;
						int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
						if (rightMovSteps > 0) {
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
						} else if (rightMovSteps < 0) {
							int code1 = getLeftMovingCode(kMod8);
							int code2 = getRightMovingCode(kMod8, resiBitsLength);
							int leftMovSteps = -rightMovSteps;
							rightMovSteps = 8 - leftMovSteps;
							resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
							p++;
							resiBits = resiBits
									| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
						} else // rightMovSteps == 0
						{
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code);
							p++;
						}
						k += resiBitsLength;
					}

					// recover the exact data
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}
			}
		}
	}

	free(leadNum);
	free(type);
	return;
}


void decompressDataSeries_float_4D(float** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageF* tdps)
{
	updateQuantizationInfo(tdps->intervals);
	size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	size_t dataSeriesLength = r1*r2*r3*r4;
	size_t r234 = r2*r3*r4;
	size_t r34 = r3*r4;
//	printf ("%d %d %d %d\n", r1, r2, r3, r4);
	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);
	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);

	unsigned char preBytes[4];
	unsigned char curBytes[4];

	memset(preBytes, 0, 4);
	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits;
	unsigned char leadingNum;
	float medianValue, exactData, predValue;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;

	float pred1D, pred2D, pred3D;
	size_t ii, jj, kk, ll;
	size_t index;

	for (ll = 0; ll < r1; ll++)
	{

		///////////////////////////	Process layer-0 ///////////////////////////
		/* Process Row-0 data 0*/
		index = ll*r234;

		// compute resiBits
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
			} else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data
		memset(curBytes, 0, 4);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}
		exactData = bytesToFloat(curBytes);
		(*data)[index] = exactData + medianValue;
		memcpy(preBytes,curBytes,4);

		/* Process Row-0, data 1 */
		index = ll*r234+1;

		pred1D = (*data)[index-1];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process Row-0, data 2 --> data r4-1 */
		for (jj = 2; jj < r4; jj++)
		{
			index = ll*r234+jj;

			pred1D = 2*(*data)[index-1] - (*data)[index-2];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}

		/* Process Row-1 --> Row-r3-1 */
		for (ii = 1; ii < r3; ii++)
		{
			/* Process row-ii data 0 */
			index = ll*r234+ii*r4;

			pred1D = (*data)[index-r4];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}

			/* Process row-ii data 1 --> r4-1*/
			for (jj = 1; jj < r4; jj++)
			{
				index = ll*r234+ii*r4+jj;

				pred2D = (*data)[index-1] + (*data)[index-r4] - (*data)[index-r4-1];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
				}
				else
				{
					// compute resiBits
					resiBits = 0;
					if (resiBitsLength != 0) {
						int kMod8 = k % 8;
						int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
						if (rightMovSteps > 0) {
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
						} else if (rightMovSteps < 0) {
							int code1 = getLeftMovingCode(kMod8);
							int code2 = getRightMovingCode(kMod8, resiBitsLength);
							int leftMovSteps = -rightMovSteps;
							rightMovSteps = 8 - leftMovSteps;
							resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
							p++;
							resiBits = resiBits
									| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
						} else // rightMovSteps == 0
						{
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code);
							p++;
						}
						k += resiBitsLength;
					}

					// recover the exact data
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}
			}
		}

		///////////////////////////	Process layer-1 --> layer-r2-1 ///////////////////////////

		for (kk = 1; kk < r2; kk++)
		{
			/* Process Row-0 data 0*/
			index = ll*r234+kk*r34;

			pred1D = (*data)[index-r34];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}

			/* Process Row-0 data 1 --> data r4-1 */
			for (jj = 1; jj < r4; jj++)
			{
				index = ll*r234+kk*r34+jj;

				pred2D = (*data)[index-1] + (*data)[index-r34] - (*data)[index-r34-1];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
				}
				else
				{
					// compute resiBits
					resiBits = 0;
					if (resiBitsLength != 0) {
						int kMod8 = k % 8;
						int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
						if (rightMovSteps > 0) {
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
						} else if (rightMovSteps < 0) {
							int code1 = getLeftMovingCode(kMod8);
							int code2 = getRightMovingCode(kMod8, resiBitsLength);
							int leftMovSteps = -rightMovSteps;
							rightMovSteps = 8 - leftMovSteps;
							resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
							p++;
							resiBits = resiBits
									| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
						} else // rightMovSteps == 0
						{
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code);
							p++;
						}
						k += resiBitsLength;
					}

					// recover the exact data
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}
			}

			/* Process Row-1 --> Row-r3-1 */
			for (ii = 1; ii < r3; ii++)
			{
				/* Process Row-i data 0 */
				index = ll*r234+kk*r34+ii*r4;

				pred2D = (*data)[index-r4] + (*data)[index-r34] - (*data)[index-r34-r4];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
				}
				else
				{
					// compute resiBits
					resiBits = 0;
					if (resiBitsLength != 0) {
						int kMod8 = k % 8;
						int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
						if (rightMovSteps > 0) {
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
						} else if (rightMovSteps < 0) {
							int code1 = getLeftMovingCode(kMod8);
							int code2 = getRightMovingCode(kMod8, resiBitsLength);
							int leftMovSteps = -rightMovSteps;
							rightMovSteps = 8 - leftMovSteps;
							resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
							p++;
							resiBits = resiBits
									| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
						} else // rightMovSteps == 0
						{
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code);
							p++;
						}
						k += resiBitsLength;
					}

					// recover the exact data
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}

				/* Process Row-i data 1 --> data r4-1 */
				for (jj = 1; jj < r4; jj++)
				{
					index = ll*r234+kk*r34+ii*r4+jj;

					pred3D = (*data)[index-1] + (*data)[index-r4] + (*data)[index-r34]
							- (*data)[index-r4-1] - (*data)[index-r34-r4] - (*data)[index-r34-1] + (*data)[index-r34-r4-1];


					type_ = type[index];
					if (type_ != 0)
					{
						(*data)[index] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
					}
					else
					{
						// compute resiBits
						resiBits = 0;
						if (resiBitsLength != 0) {
							int kMod8 = k % 8;
							int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
							if (rightMovSteps > 0) {
								int code = getRightMovingCode(kMod8, resiBitsLength);
								resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
							} else if (rightMovSteps < 0) {
								int code1 = getLeftMovingCode(kMod8);
								int code2 = getRightMovingCode(kMod8, resiBitsLength);
								int leftMovSteps = -rightMovSteps;
								rightMovSteps = 8 - leftMovSteps;
								resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
								p++;
								resiBits = resiBits
										| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
							} else // rightMovSteps == 0
							{
								int code = getRightMovingCode(kMod8, resiBitsLength);
								resiBits = (tdps->residualMidBits[p] & code);
								p++;
							}
							k += resiBitsLength;
						}

						// recover the exact data
						memset(curBytes, 0, 4);
						leadingNum = leadNum[l++];
						memcpy(curBytes, preBytes, leadingNum);
						for (j = leadingNum; j < reqBytesLength; j++)
							curBytes[j] = tdps->exactMidBytes[curByteIndex++];
						if (resiBitsLength != 0) {
							unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
							curBytes[reqBytesLength] = resiByte;
						}

						exactData = bytesToFloat(curBytes);
						(*data)[index] = exactData + medianValue;
						memcpy(preBytes,curBytes,4);
					}
				}
			}

		}
	}

	free(leadNum);
	free(type);
	return;
}

void getSnapshotData_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps, int errBoundMode)
{	
	SZ_Reset(allNodes, stateNum);
	size_t i;

	if (tdps->rtypeArray == NULL) {
		if(errBoundMode < PW_REL)
			decompressDataSeries_float_1D(data, dataSeriesLength, tdps);
		else 
		{
			//decompressDataSeries_float_1D_pwr(data, dataSeriesLength, tdps);
			// decompressDataSeries_float_1D_pwrgroup(data, dataSeriesLength, tdps);
			decompressDataSeries_float_1D_pwr_pre_log(data, dataSeriesLength, tdps);
		}
		return;
	} else {
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		// insert the reserved values
		//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
		//		dataSeriesLength, rtypeArray);
		int* rtypes;
		int validLength = computeBitNumRequired(dataSeriesLength);
		decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
		size_t count = 0;
		for (i = 0; i < dataSeriesLength; i++) {
			if (rtypes[i] == 1)
				(*data)[i] = tdps->reservedValue;
			else
				count++;
		}
		// get the decompressed data
		float* decmpData;
		if(errBoundMode < PW_REL)
			decompressDataSeries_float_1D(&decmpData, dataSeriesLength, tdps);
		else 
			decompressDataSeries_float_1D_pwr(&decmpData, dataSeriesLength, tdps);
		// insert the decompressed data
		size_t k = 0;
		for (i = 0; i < dataSeriesLength; i++) {
			if (rtypes[i] == 0) {
				(*data)[i] = decmpData[k++];
			}
		}
		free(decmpData);
		free(rtypes);
	}
}

void getSnapshotData_float_2D(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps, int errBoundMode) 
{
	SZ_Reset(allNodes, stateNum);
	size_t i;
	size_t dataSeriesLength = r1*r2;
	if (tdps->rtypeArray == NULL) {
		if(errBoundMode < PW_REL)
			decompressDataSeries_float_2D(data, r1, r2, tdps);
		else 
			// decompressDataSeries_float_2D_pwr(data, r1, r2, tdps);
			decompressDataSeries_float_2D_pwr_pre_log(data, r1, r2, tdps);
		return;
	} else {
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		// insert the reserved values
		//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
		//		dataSeriesLength, rtypeArray);
		int* rtypes;
		int validLength = computeBitNumRequired(dataSeriesLength);
		decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
		size_t count = 0;
		for (i = 0; i < dataSeriesLength; i++) {
			if (rtypes[i] == 1)
				(*data)[i] = tdps->reservedValue;
			else
				count++;
		}
		// get the decompressed data
		float* decmpData;
		if(errBoundMode < PW_REL)
			decompressDataSeries_float_2D(&decmpData, r1, r2, tdps);
		else 
			decompressDataSeries_float_2D_pwr(&decmpData, r1, r2, tdps);
		// insert the decompressed data
		size_t k = 0;
		for (i = 0; i < dataSeriesLength; i++) {
			if (rtypes[i] == 0) {
				(*data)[i] = decmpData[k++];
			}
		}
		free(decmpData);
		free(rtypes);
	}
}

void getSnapshotData_float_3D(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps, int errBoundMode)
{
	SZ_Reset(allNodes, stateNum);
	size_t i;
	size_t dataSeriesLength = r1*r2*r3;
	if (tdps->rtypeArray == NULL) {
		if(errBoundMode < PW_REL)
			decompressDataSeries_float_3D(data, r1, r2, r3, tdps);
		else 
			// decompressDataSeries_float_3D_pwr(data, r1, r2, r3, tdps);
			decompressDataSeries_float_3D_pwr_pre_log(data, r1, r2, r3, tdps);
		return;
	} else {
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		// insert the reserved values
		//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
		//		dataSeriesLength, rtypeArray);
		int* rtypes;
		int validLength = computeBitNumRequired(dataSeriesLength);
		decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
		size_t count = 0;
		for (i = 0; i < dataSeriesLength; i++) {
			if (rtypes[i] == 1)
				(*data)[i] = tdps->reservedValue;
			else
				count++;
		}
		// get the decompressed data
		float* decmpData;
		if(errBoundMode < PW_REL)
			decompressDataSeries_float_3D(&decmpData, r1, r2, r3, tdps);
		else 
			decompressDataSeries_float_3D_pwr(&decmpData, r1, r2, r3, tdps);
		// insert the decompressed data
		size_t k = 0;
		for (i = 0; i < dataSeriesLength; i++) {
			if (rtypes[i] == 0) {
				(*data)[i] = decmpData[k++];
			}
		}
		free(decmpData);
		free(rtypes);
	}
}

void getSnapshotData_float_4D(float** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageF* tdps, int errBoundMode)
{
	SZ_Reset(allNodes, stateNum);
	size_t i;
	size_t dataSeriesLength = r1*r2*r3*r4;

	if (tdps->rtypeArray == NULL) {
		if(errBoundMode < PW_REL)
			decompressDataSeries_float_4D(data, r1, r2, r3, r4, tdps);
		else
			decompressDataSeries_float_3D_pwr_pre_log(data, r1*r2, r3, r4, tdps);
			//ToDO
			//decompressDataSeries_float_4D_pwr(data, r1, r2, r3, r4, tdps);
		return;
	} else {
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		int* rtypes;
		int validLength = computeBitNumRequired(dataSeriesLength);
		decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
		size_t count = 0;
		for (i = 0; i < dataSeriesLength; i++) {
			if (rtypes[i] == 1)
				(*data)[i] = tdps->reservedValue;
			else
				count++;
		}
		// get the decompressed data
		float* decmpData;
		if(errBoundMode < PW_REL)
			decompressDataSeries_float_4D(&decmpData, r1, r2, r3, r4, tdps);
		else
			decompressDataSeries_float_3D_pwr_pre_log(&decmpData, r1*r2, r3, r4, tdps);
			//ToDO
			//decompressDataSeries_float_4D_pwr(&decompData, r1, r2, r3, r4, tdps);
		// insert the decompressed data
		size_t k = 0;
		for (i = 0; i < dataSeriesLength; i++) {
			if (rtypes[i] == 0) {
				(*data)[i] = decmpData[k++];
			}
		}
		free(decmpData);
		free(rtypes);
	}
}

size_t decompressDataSeries_float_1D_RA_block(float * data, float mean, size_t dim_0, size_t block_dim_0, double realPrecision, int * type, float * unpredictable_data){

	size_t unpredictable_count = 0;
	
	float * cur_data_pos = data;
	size_t type_index = 0;
	int type_;
	float last_over_thres = mean;
	for(size_t i=0; i<block_dim_0; i++){
		type_ = type[type_index];
		if(type_ == 0){
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
			last_over_thres = cur_data_pos[0];
		}
		else{
			cur_data_pos[0] = last_over_thres + 2 * (type_ - intvRadius) * realPrecision;
			last_over_thres = cur_data_pos[0];
		}

		type_index ++;
		cur_data_pos ++;
	}

	return unpredictable_count;
}

size_t decompressDataSeries_float_1D_RA_block_1D_pred(float * data, float mean, size_t dim_0, size_t block_dim_0, double realPrecision, int * type, float * unpredictable_data){
	
	size_t unpredictable_count = 0;
	
	float * cur_data_pos = data;
	size_t type_index = 0;
	int type_;
	float last_over_thres = mean;
	for(size_t i=0; i<block_dim_0; i++){
		type_ = type[type_index];
		if(type_ == 0){
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
			last_over_thres = cur_data_pos[0];
		}
		else if(type_ == 1){
			cur_data_pos[0] = mean;
		}
		else{
			cur_data_pos[0] = last_over_thres + 2 * (type_ - intvRadius) * realPrecision;
			last_over_thres = cur_data_pos[0];
		}

		type_index ++;
		cur_data_pos ++;
	}

	return unpredictable_count;
}

size_t decompressDataSeries_float_3D_RA_block_2_layers(float * data, float mean, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	// printf("SZ_compress_float_3D_MDQ_RA_block real dim: %d %d %d\n", real_block_dims[0], real_block_dims[1], real_block_dims[2]);
	// fflush(stdout);

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = data;
	float * last_row_pos;
	float curData;
	float pred1D, pred2D, pred3D;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	int type_;
	// Process Row-0 data 0
	pred1D = mean;
	type_ = type[0];
	// printf("Type 0 %d, mean %.4f\n", type_, mean);
	if (type_ != 0){
		cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
	}

	/* Process Row-0 data 1*/
	pred1D = cur_data_pos[0];
	type_ = type[1];
	if (type_ != 0){
		cur_data_pos[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[1] = unpredictable_data[unpredictable_count ++];
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*cur_data_pos[j-1] - cur_data_pos[j-2];
		type_ = type[j];
		if (type_ != 0){
			cur_data_pos[j] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
		}
	}

	last_row_pos = cur_data_pos;
	cur_data_pos += dim1_offset;
	// printf("SZ_compress_float_3D_MDQ_RA_block row 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = last_row_pos[0];
		type_ = type[index];
		if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = cur_data_pos[j-1] + last_row_pos[j] - last_row_pos[j-1];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f last_row_data %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], last_row_pos[j], last_row_pos[j-1], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	// printf("SZ_compress_float_3D_MDQ_RA_block layer 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);
	// exit(0);

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer %d done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = cur_data_pos[- dim0_offset];
		type_ = type[index];
		if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = cur_data_pos[j-1] + cur_data_pos[j - dim0_offset] - cur_data_pos[j - 1 - dim0_offset];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], cur_data_pos[j - dim0_offset], cur_data_pos[j - 1 - dim0_offset], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;

		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row 0 done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }

	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			// if(idx == 63 && idy == 63 && idz == 63){
			// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row %d done, cur_data_pos: %ld\n", i-1, cur_data_pos - data);
			// 	fflush(stdout);
			// }
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = last_row_pos[0] + cur_data_pos[- dim0_offset] - last_row_pos[- dim0_offset];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[0] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
//				if(k==63&&i==43&&j==27)
//					printf("i=%d\n", i);
				//index = k*r2*r3 + i*r3 + j;
				index ++;
				index2D = i*r3 + j;
				if(i >= 2 && j>= 2 && k >= 2){
					// pred3D = cur_data_pos[j-1] + last_row_pos[j]+ cur_data_pos[j - dim0_offset] - last_row_pos[j-1] - last_row_pos[j - dim0_offset] - cur_data_pos[j-1 - dim0_offset] + last_row_pos[j-1 - dim0_offset];
					float * tmp_pos = cur_data_pos + j;
					pred3D = - tmp_pos[-2*dim0_offset - 2*dim1_offset - 2]
							 + 2 * (tmp_pos[-2*dim0_offset - 2*dim1_offset - 1] + tmp_pos[-2*dim0_offset - dim1_offset - 2] + tmp_pos[- dim0_offset - 2*dim1_offset - 2])
							 - 4 * (tmp_pos[-2*dim0_offset - dim1_offset - 1] + tmp_pos[- dim0_offset - 2*dim1_offset - 1] + tmp_pos[-dim0_offset - dim1_offset - 2])
							 - (tmp_pos[-2*dim0_offset] + tmp_pos[-dim0_offset - 2*dim1_offset] + tmp_pos[- 2])
							 - (tmp_pos[-2*dim0_offset - 2*dim1_offset] + tmp_pos[-2*dim0_offset - 2] + tmp_pos[- 2*dim1_offset - 2])
							 - 4 * (tmp_pos[-dim0_offset - dim1_offset] + tmp_pos[-dim0_offset - 1] + tmp_pos[- dim1_offset - 1])
							 + 2 * (tmp_pos[-dim0_offset] + tmp_pos[- dim1_offset] + tmp_pos[- 1])
							 + 2 * (tmp_pos[-2*dim0_offset - dim1_offset] + tmp_pos[-2*dim0_offset - 1] + tmp_pos[-dim0_offset - 2*dim1_offset] + tmp_pos[-dim0_offset - 2] + tmp_pos[- 2*dim1_offset - 1] + tmp_pos[- dim1_offset -2])
							 + 8 * (tmp_pos[-dim0_offset - dim1_offset - 1]);
					// if(i==2 && j==2 && k==2){
					// 	printf("pred3D: %.4f\n", pred3D);
					// 	for(int i=0; i<3; i++){
					// 		for(int j=0; j<3; j++){
					// 			for(int k=0; k<3; k++){
					// 				printf("%.4f ", tmp_pos[-i*dim0_offset - j*dim1_offset - k]);
					// 			}
					// 		}
					// 	}
					// 	printf("\n");
					// }
					// exit(0);
				}
				else
					pred3D = cur_data_pos[j-1] + last_row_pos[j]+ cur_data_pos[j - dim0_offset] - last_row_pos[j-1] - last_row_pos[j - dim0_offset] - cur_data_pos[j-1 - dim0_offset] + last_row_pos[j-1 - dim0_offset];
				type_ = type[index];
				if (type_ != 0){
					cur_data_pos[j] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
				}
				else{
					cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
				}
			}
			last_row_pos = cur_data_pos;
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
	}

	return unpredictable_count;
}

size_t decompressDataSeries_float_3D_RA_block_no_mean(float * data, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	// printf("SZ_compress_float_3D_MDQ_RA_block real dim: %d %d %d\n", real_block_dims[0], real_block_dims[1], real_block_dims[2]);
	// fflush(stdout);

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = data;
	float * last_row_pos;
	float curData;
	float pred1D, pred2D, pred3D;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	int type_;
	// Process Row-0 data 0
	type_ = type[0];
	// printf("Type 0 %d, mean %.4f\n", type_, mean);
	if (type_ != 0){
		printf("first data should be unpredictable!\n");
		exit(0);
	}
	else{
		cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
	}

	/* Process Row-0 data 1*/
	pred1D = cur_data_pos[0];
	type_ = type[1];
	if (type_ != 0){
		cur_data_pos[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[1] = unpredictable_data[unpredictable_count ++];
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*cur_data_pos[j-1] - cur_data_pos[j-2];
		type_ = type[j];
		if (type_ != 0){
			cur_data_pos[j] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
		}
	}

	last_row_pos = cur_data_pos;
	cur_data_pos += dim1_offset;
	// printf("SZ_compress_float_3D_MDQ_RA_block row 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = last_row_pos[0];
		type_ = type[index];
		if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = cur_data_pos[j-1] + last_row_pos[j] - last_row_pos[j-1];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f last_row_data %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], last_row_pos[j], last_row_pos[j-1], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	// printf("SZ_compress_float_3D_MDQ_RA_block layer 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);
	// exit(0);

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer %d done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = cur_data_pos[- dim0_offset];
		type_ = type[index];
		if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = cur_data_pos[j-1] + cur_data_pos[j - dim0_offset] - cur_data_pos[j - 1 - dim0_offset];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], cur_data_pos[j - dim0_offset], cur_data_pos[j - 1 - dim0_offset], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;

		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row 0 done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }

	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			// if(idx == 63 && idy == 63 && idz == 63){
			// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row %d done, cur_data_pos: %ld\n", i-1, cur_data_pos - data);
			// 	fflush(stdout);
			// }
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = last_row_pos[0] + cur_data_pos[- dim0_offset] - last_row_pos[- dim0_offset];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[0] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
//				if(k==63&&i==43&&j==27)
//					printf("i=%d\n", i);
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				pred3D = cur_data_pos[j-1] + last_row_pos[j]+ cur_data_pos[j - dim0_offset] - last_row_pos[j-1] - last_row_pos[j - dim0_offset] - cur_data_pos[j-1 - dim0_offset] + last_row_pos[j-1 - dim0_offset];
				type_ = type[index];
				if (type_ != 0){
					cur_data_pos[j] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
				}
				else{
					cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
				}
			}
			last_row_pos = cur_data_pos;
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
	}

	return unpredictable_count;
}

size_t decompressDataSeries_float_3D_RA_block(float * data, float mean, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	// printf("SZ_compress_float_3D_MDQ_RA_block real dim: %d %d %d\n", real_block_dims[0], real_block_dims[1], real_block_dims[2]);
	// fflush(stdout);

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = data;
	float * last_row_pos;
	float curData;
	float pred1D, pred2D, pred3D;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	int type_;
	// Process Row-0 data 0
	pred1D = mean;
	type_ = type[0];
	// printf("Type 0 %d, mean %.4f\n", type_, mean);
	if (type_ != 0){
		cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
	}

	/* Process Row-0 data 1*/
	pred1D = cur_data_pos[0];
	type_ = type[1];
	if (type_ != 0){
		cur_data_pos[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[1] = unpredictable_data[unpredictable_count ++];
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*cur_data_pos[j-1] - cur_data_pos[j-2];
		type_ = type[j];
		if (type_ != 0){
			cur_data_pos[j] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
		}
	}

	last_row_pos = cur_data_pos;
	cur_data_pos += dim1_offset;
	// printf("SZ_compress_float_3D_MDQ_RA_block row 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = last_row_pos[0];
		type_ = type[index];
		if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = cur_data_pos[j-1] + last_row_pos[j] - last_row_pos[j-1];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f last_row_data %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], last_row_pos[j], last_row_pos[j-1], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	// printf("SZ_compress_float_3D_MDQ_RA_block layer 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);
	// exit(0);

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer %d done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = cur_data_pos[- dim0_offset];
		type_ = type[index];
		if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = cur_data_pos[j-1] + cur_data_pos[j - dim0_offset] - cur_data_pos[j - 1 - dim0_offset];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], cur_data_pos[j - dim0_offset], cur_data_pos[j - 1 - dim0_offset], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;

		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row 0 done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }

	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			// if(idx == 63 && idy == 63 && idz == 63){
			// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row %d done, cur_data_pos: %ld\n", i-1, cur_data_pos - data);
			// 	fflush(stdout);
			// }
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = last_row_pos[0] + cur_data_pos[- dim0_offset] - last_row_pos[- dim0_offset];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[0] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
//				if(k==63&&i==43&&j==27)
//					printf("i=%d\n", i);
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				pred3D = cur_data_pos[j-1] + last_row_pos[j]+ cur_data_pos[j - dim0_offset] - last_row_pos[j-1] - last_row_pos[j - dim0_offset] - cur_data_pos[j-1 - dim0_offset] + last_row_pos[j-1 - dim0_offset];
				type_ = type[index];
				if (type_ != 0){
					cur_data_pos[j] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
				}
				else{
					cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
				}
			}
			last_row_pos = cur_data_pos;
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
	}

	return unpredictable_count;
}

size_t decompressDataSeries_float_3D_RA_block_adaptive(float * data, float mean, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	// printf("SZ_compress_float_3D_MDQ_RA_block real dim: %d %d %d\n", real_block_dims[0], real_block_dims[1], real_block_dims[2]);
	// fflush(stdout);

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = data;
	float * last_row_pos;
	float curData;
	float pred1D, pred2D, pred3D;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	int type_;
	// Process Row-0 data 0
	pred1D = mean;
	type_ = type[0];
	// printf("Type 0 %d, mean %.4f\n", type_, mean);
	if (type_ != 0){
		cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
	}

	/* Process Row-0 data 1*/
	pred1D = cur_data_pos[0];
	type_ = type[1];
	if (type_ != 0){
		cur_data_pos[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[1] = unpredictable_data[unpredictable_count ++];
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*cur_data_pos[j-1] - cur_data_pos[j-2];
		type_ = type[j];
		if (type_ != 0){
			cur_data_pos[j] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
		}
	}

	last_row_pos = cur_data_pos;
	cur_data_pos += dim1_offset;
	// printf("SZ_compress_float_3D_MDQ_RA_block row 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = last_row_pos[0];
		type_ = type[index];
		if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = cur_data_pos[j-1] + last_row_pos[j] - last_row_pos[j-1];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f last_row_data %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], last_row_pos[j], last_row_pos[j-1], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	// printf("SZ_compress_float_3D_MDQ_RA_block layer 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);
	// exit(0);

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer %d done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = cur_data_pos[- dim0_offset];
		type_ = type[index];
		if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = cur_data_pos[j-1] + cur_data_pos[j - dim0_offset] - cur_data_pos[j - 1 - dim0_offset];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], cur_data_pos[j - dim0_offset], cur_data_pos[j - 1 - dim0_offset], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;

		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row 0 done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }

	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			// if(idx == 63 && idy == 63 && idz == 63){
			// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row %d done, cur_data_pos: %ld\n", i-1, cur_data_pos - data);
			// 	fflush(stdout);
			// }
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = last_row_pos[0] + cur_data_pos[- dim0_offset] - last_row_pos[- dim0_offset];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[0] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
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

//				if(k==63&&i==43&&j==27)
//					printf("i=%d\n", i);
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				pred3D = cur_data_pos[j-1] + last_row_pos[j]+ cur_data_pos[j - dim0_offset] - last_row_pos[j-1] - last_row_pos[j - dim0_offset] - cur_data_pos[j-1 - dim0_offset] + last_row_pos[j-1 - dim0_offset];
				type_ = type[index];
				if (type_ != 0){
					cur_data_pos[j] = pred3D + 2 * (type_ - intvRadius) * tmp_realPrecision;
				}
				else{
					cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
				}
			}
			last_row_pos = cur_data_pos;
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
	}

	return unpredictable_count;
}

size_t decompressDataSeries_float_3D_RA_block_3D_pred_multi_means(float * data, unsigned int mean_count, float * means, float dense_pos, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, int * type, float * unpredictable_data){

	float sum = 0.0;
	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	// printf("SZ_compress_float_3D_MDQ_RA_block real dim: %d %d %d\n", real_block_dims[0], real_block_dims[1], real_block_dims[2]);
	// fflush(stdout);

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = data;
	float * last_row_pos;
	float curData;
	float pred1D, pred2D, pred3D;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	int type_;
	// Process Row-0 data 0
	pred1D = dense_pos;
	type_ = type[0];
	// printf("Type 0 %d, mean %.4f\n", type_, mean);
	if(type_ == 0){
		cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
	}
	else if(type_ <= mean_count){
		cur_data_pos[0] = means[type_ - 1];
	}
	else{
		cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}

	/* Process Row-0 data 1*/
	pred1D = cur_data_pos[0];
	type_ = type[1];
	if(type_ == 0){
		cur_data_pos[1] = unpredictable_data[unpredictable_count ++];
	}
	else if (type_ <= mean_count){
		cur_data_pos[1] = means[type_ - 1];
	}
	else{
		cur_data_pos[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*cur_data_pos[j-1] - cur_data_pos[j-2];
		type_ = type[j];
		if(type_ == 0){
			cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
		}
		else if (type_ <= mean_count){
			cur_data_pos[j] = means[type_ - 1];
		}
		else{
			cur_data_pos[j] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
	}

	last_row_pos = cur_data_pos;
	cur_data_pos += dim1_offset;
	// printf("SZ_compress_float_3D_MDQ_RA_block row 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = last_row_pos[0];
		type_ = type[index];
		if(type_ == 0){
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		else if (type_ <= mean_count){
			cur_data_pos[0] = means[type_ - 1];
		}
		else{
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = cur_data_pos[j-1] + last_row_pos[j] - last_row_pos[j-1];
			type_ = type[index];
			if(type_ == 0){
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			else if (type_ <= mean_count){
				cur_data_pos[j] = means[type_ - 1];
			}
			else{
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			// printf("pred2D %.2f cur_data %.2f last_row_data %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], last_row_pos[j], last_row_pos[j-1], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	// printf("SZ_compress_float_3D_MDQ_RA_block layer 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);
	// exit(0);

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer %d done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = cur_data_pos[- dim0_offset];
		type_ = type[index];
		if(type_ == 0){
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		else if (type_ <= mean_count){
			cur_data_pos[0] = means[type_ - 1];
		}
		else{
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = cur_data_pos[j-1] + cur_data_pos[j - dim0_offset] - cur_data_pos[j - 1 - dim0_offset];
			type_ = type[index];
			if(type_ == 0){
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			else if (type_ <= mean_count){
				cur_data_pos[j] = means[type_ - 1];
			}
			else{
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			// printf("pred2D %.2f cur_data %.2f %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], cur_data_pos[j - dim0_offset], cur_data_pos[j - 1 - dim0_offset], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;

		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row 0 done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }

	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			// if(idx == 63 && idy == 63 && idz == 63){
			// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row %d done, cur_data_pos: %ld\n", i-1, cur_data_pos - data);
			// 	fflush(stdout);
			// }
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = last_row_pos[0] + cur_data_pos[- dim0_offset] - last_row_pos[- dim0_offset];
			type_ = type[index];
			if(type_ == 0){
				cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
			}
			else if (type_ <= mean_count){
				cur_data_pos[0] = means[type_ - 1];
			}
			else{
				cur_data_pos[0] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
//				if(k==63&&i==43&&j==27)
//					printf("i=%d\n", i);
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				pred3D = cur_data_pos[j-1] + last_row_pos[j]+ cur_data_pos[j - dim0_offset] - last_row_pos[j-1] - last_row_pos[j - dim0_offset] - cur_data_pos[j-1 - dim0_offset] + last_row_pos[j-1 - dim0_offset];
				type_ = type[index];
				if(type_ == 0){
					cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
				}
				else if (type_ <= mean_count){
					cur_data_pos[j] = means[type_ - 1];
				}
				else{
					cur_data_pos[j] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
				}
			}
			last_row_pos = cur_data_pos;
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
	}

	return unpredictable_count;
}

size_t decompressDataSeries_float_2D_RA_block_2D_pred(float * data, float mean, size_t dim_0, size_t dim_1, size_t block_dim_0, size_t block_dim_1, double realPrecision, int * type, float * unpredictable_data){

	float * data_pos;
	size_t dim0_offset = dim_1;
	// printf("SZ_compress_float_3D_MDQ_RA_block real dim: %d %d %d\n", real_block_dims[0], real_block_dims[1], real_block_dims[2]);
	// fflush(stdout);

	size_t unpredictable_count = 0;
	size_t r1, r2;
	r1 = block_dim_0;
	r2 = block_dim_1;

	float * cur_data_pos = data;
	float * last_row_pos;
	float curData;
	float pred1D, pred2D;
	double diff;
	size_t i, j;
	int type_;
	// Process Row-0 data 0
	pred1D = mean;
	type_ = type[0];
	// printf("Type 0 %d, mean %.4f\n", type_, mean);
	if(type_ == 1){
		cur_data_pos[0] = mean;
	}
	else if (type_ != 0){
		cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
	}

	/* Process Row-0 data 1*/
	pred1D = cur_data_pos[0];
	type_ = type[1];
	if(type_ == 1){
		cur_data_pos[1] = mean;
	}
	else if (type_ != 0){
		cur_data_pos[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[1] = unpredictable_data[unpredictable_count ++];
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r2; j++){
		pred1D = 2*cur_data_pos[j-1] - cur_data_pos[j-2];
		type_ = type[j];
		if(type_ == 1){
			cur_data_pos[j] = mean;
		}
		else if (type_ != 0){
			cur_data_pos[j] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
		}
	}

	last_row_pos = cur_data_pos;
	cur_data_pos += dim0_offset;
	// printf("SZ_compress_float_3D_MDQ_RA_block row 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r1; i++)
	{
		/* Process row-i data 0 */
		index = i*r2;	
		type_ = type[index];
		if(type_ == 1){
			cur_data_pos[0] = mean;
		}
		else if (type_ != 0){
			pred1D = last_row_pos[0];
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r2; j++)
		{
			index = i*r2+j;
			pred2D = cur_data_pos[j-1] + last_row_pos[j] - last_row_pos[j-1];
			type_ = type[index];
			if(type_ == 1){
				cur_data_pos[j] = mean;
			}
			else if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f last_row_data %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], last_row_pos[j], last_row_pos[j-1], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim0_offset;
	}
	return unpredictable_count;
}

size_t decompressDataSeries_float_2D_RA_block(float * data, float mean, size_t dim_0, size_t dim_1, size_t block_dim_0, size_t block_dim_1, double realPrecision, int * type, float * unpredictable_data){

	float * data_pos;
	size_t dim0_offset = dim_1;
	// printf("SZ_compress_float_3D_MDQ_RA_block real dim: %d %d %d\n", real_block_dims[0], real_block_dims[1], real_block_dims[2]);
	// fflush(stdout);

	size_t unpredictable_count = 0;
	size_t r1, r2;
	r1 = block_dim_0;
	r2 = block_dim_1;

	float * cur_data_pos = data;
	float * last_row_pos;
	float curData;
	float pred1D, pred2D;
	double diff;
	size_t i, j;
	int type_;
	// Process Row-0 data 0
	pred1D = mean;
	type_ = type[0];
	// printf("Type 0 %d, mean %.4f\n", type_, mean);
	if (type_ != 0){
		cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
	}

	/* Process Row-0 data 1*/
	pred1D = cur_data_pos[0];
	type_ = type[1];
	if (type_ != 0){
		cur_data_pos[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[1] = unpredictable_data[unpredictable_count ++];
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r2; j++){
		pred1D = 2*cur_data_pos[j-1] - cur_data_pos[j-2];
		type_ = type[j];
		if (type_ != 0){
			cur_data_pos[j] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
		}
	}

	last_row_pos = cur_data_pos;
	cur_data_pos += dim0_offset;
	// printf("SZ_compress_float_3D_MDQ_RA_block row 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r1; i++)
	{
		/* Process row-i data 0 */
		index = i*r2;	
		type_ = type[index];
		if (type_ != 0){
			pred1D = last_row_pos[0];
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r2; j++)
		{
			index = i*r2+j;
			pred2D = cur_data_pos[j-1] + last_row_pos[j] - last_row_pos[j-1];
			type_ = type[index];
			if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f last_row_data %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], last_row_pos[j], last_row_pos[j-1], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim0_offset;
	}
	return unpredictable_count;
}

void decompressDataSeries_float_2D_nonblocked(float** data, size_t r1, size_t r2, unsigned char* comp_data){
	// calculate block dims
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t num_elements = r1 * r2;
	size_t dim0_offset = r2;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2;
	// intvRadius = (int)((tdps->intervals - 1)/ 2);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;

	unsigned int unpred_count = *((unsigned int *)comp_data_pos);
	comp_data_pos += 4;
	float * unpredictable_data = (float *) comp_data_pos;
	comp_data_pos += unpred_count * sizeof(float);
	
	float mean = *((float *) comp_data_pos);
	comp_data_pos += 4;

	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	int * type = (int *) malloc(num_elements * sizeof(int));
	float * data_pos = *data;
	decode(comp_data_pos, num_elements, root, type);
	size_t cur_unpred_count = decompressDataSeries_float_2D_RA_block_2D_pred(data_pos, mean, r1, r2, r1, r2, realPrecision, type, unpredictable_data);

	if(cur_unpred_count != unpred_count){
		printf("Check bugs, unpredictable_count is not the same: %d %d\n", unpred_count, cur_unpred_count);
		for(size_t i=0; i<512; i++){
			printf("%d ", type[i]);
		}
		printf("\n");
		exit(0);
	}
	free(type);

}


void decompressDataSeries_float_2D_RA(float** data, size_t r1, size_t r2, unsigned char* comp_data){
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t dim0_offset = r2;
	size_t num_elements = r1 * r2;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	// calculate block dims
	size_t num_x, num_y;
	SZ_COMPUTE_2D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	SZ_COMPUTE_2D_NUMBER_OF_BLOCKS(r2, num_y, block_size);

	size_t split_index_x, split_index_y;
	size_t early_blockcount_x, early_blockcount_y;
	size_t late_blockcount_x, late_blockcount_y;
	SZ_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	SZ_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y;
	size_t num_blocks = num_x * num_y;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2;
	// intvRadius = (int)((tdps->intervals - 1)/ 2);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;
	unsigned short * block_pos = (unsigned short *) comp_data_pos;
	// skip block index here
	comp_data_pos += num_blocks * sizeof(unsigned short);
	unsigned short * unpred_count = (unsigned short *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(unsigned short);
	float * mean_pos = (float *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(float);

	size_t unpredictable_count;

	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	int * type = (int *) malloc(max_num_block_elements * sizeof(int));
	size_t unpred_data_max_size = max_num_block_elements;
	float * unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float));
	float mean;
	unsigned char * tmp;
	size_t unpredictableEncodeSize;
	size_t index = 0;
	float * data_pos = *data;
	size_t offset_x, offset_y;
	size_t current_blockcount_x, current_blockcount_y;
	size_t cur_unpred_count;
	// printf("decompress offset to start: %ld\n", comp_data_pos - comp_data);
	// fflush(stdout);
	size_t  current_block_elements;
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
			offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
			data_pos = *data + offset_x * dim0_offset + offset_y;

			current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
			current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;

			tmp = comp_data_pos;
			mean = *mean_pos;
			// mean = 0.004558146;
			// mean = *((float *) tmp);
			// tmp += 4;
			unpredictable_count = *(unpred_count);
			if(unpredictable_count > 0){
				unpredictableEncodeSize = unpredictable_count * sizeof(float);
				memcpy(unpredictable_data, tmp, unpredictableEncodeSize);
				tmp += unpredictableEncodeSize;
			}
			current_block_elements = current_blockcount_x * current_blockcount_y;
			decode(tmp, current_block_elements, root, type);

			cur_unpred_count = decompressDataSeries_float_2D_RA_block_2D_pred(data_pos, mean, r1, r2, current_blockcount_x, current_blockcount_y, realPrecision, type, unpredictable_data);

			if(cur_unpred_count != unpredictable_count){
				printf("Check bugs, unpredictable_count is not the same: compress %d current %d\n", unpredictable_count, cur_unpred_count);
				printf("Current index: %d %d\n\n", i, j);
				size_t count = 0;
				for(int i=0; i<current_block_elements; i++){
					// printf("%d ", type[i]);
					if(type[i] == 0){
						count ++;
					}
				}
				printf("type count: %ld\n", count);
				printf("dist to decomp start: %ld\n", tmp - comp_data);
				printf("\n");
				exit(0);
			}

			comp_data_pos += *block_pos;
			block_pos ++;
			unpred_count ++;
			mean_pos ++;

		}
	}
	free(type);
	free(unpredictable_data);
}

size_t decompressDataSeries_float_3D_RA_block_3D_pred(float * data, float mean, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, int * type, float * unpredictable_data){

	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	// printf("SZ_compress_float_3D_MDQ_RA_block real dim: %d %d %d\n", real_block_dims[0], real_block_dims[1], real_block_dims[2]);
	// fflush(stdout);

	size_t unpredictable_count = 0;
	size_t r1, r2, r3;
	r1 = block_dim_0;
	r2 = block_dim_1;
	r3 = block_dim_2;

	float * cur_data_pos = data;
	float * last_row_pos;
	float curData;
	float pred1D, pred2D, pred3D;
	double diff;
	size_t i, j, k;
	size_t r23 = r2*r3;
	int type_;
	// Process Row-0 data 0
	pred1D = mean;
	type_ = type[0];
	// printf("Type 0 %d, mean %.4f\n", type_, mean);
	if(type_ == 1){
		cur_data_pos[0] = mean;
	}
	else if (type_ != 0){
		cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
	}

	/* Process Row-0 data 1*/
	pred1D = cur_data_pos[0];
	type_ = type[1];
	if(type_ == 1){
		cur_data_pos[1] = mean;
	}
	else if (type_ != 0){
		cur_data_pos[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else{
		cur_data_pos[1] = unpredictable_data[unpredictable_count ++];
	}
    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++){
		pred1D = 2*cur_data_pos[j-1] - cur_data_pos[j-2];
		type_ = type[j];
		if(type_ == 1){
			cur_data_pos[j] = mean;
		}
		else if (type_ != 0){
			cur_data_pos[j] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
		}
	}

	last_row_pos = cur_data_pos;
	cur_data_pos += dim1_offset;
	// printf("SZ_compress_float_3D_MDQ_RA_block row 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);

	/* Process Row-1 --> Row-r2-1 */
	size_t index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = last_row_pos[0];
		type_ = type[index];
		if(type_ == 1){
			cur_data_pos[0] = mean;
		}
		else if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;
			pred2D = cur_data_pos[j-1] + last_row_pos[j] - last_row_pos[j-1];
			type_ = type[index];
			if(type_ == 1){
				cur_data_pos[j] = mean;
			}
			else if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f last_row_data %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], last_row_pos[j], last_row_pos[j-1], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;
	}
	cur_data_pos += dim0_offset - r2 * dim1_offset;

	// printf("SZ_compress_float_3D_MDQ_RA_block layer 0 done, cur_data_pos: %ld\n", cur_data_pos - block_ori_data);
	// fflush(stdout);
	// exit(0);

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer %d done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }
		/* Process Row-0 data 0*/
		index = k*r23;
		pred1D = cur_data_pos[- dim0_offset];
		type_ = type[index];
		if(type_ == 1){
			cur_data_pos[0] = mean;
		}
		else if (type_ != 0){
			cur_data_pos[0] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else{
			cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
		}
	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;
			index ++;
			pred2D = cur_data_pos[j-1] + cur_data_pos[j - dim0_offset] - cur_data_pos[j - 1 - dim0_offset];
			type_ = type[index];
			if(type_ == 1){
				cur_data_pos[j] = mean;
			}
			else if (type_ != 0){
				cur_data_pos[j] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
			}
			// printf("pred2D %.2f cur_data %.2f %.2f %.2f, result %.2f\n", pred2D, cur_data_pos[j-1], cur_data_pos[j - dim0_offset], cur_data_pos[j - 1 - dim0_offset], cur_data_pos[j]);
			// getchar();
		}
		last_row_pos = cur_data_pos;
		cur_data_pos += dim1_offset;

		// if(idx == 63 && idy == 63 && idz == 63){
		// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row 0 done, cur_data_pos: %ld\n", k-1, cur_data_pos - data);
		// 	fflush(stdout);
		// }

	    /* Process Row-1 --> Row-r2-1 */
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			// if(idx == 63 && idy == 63 && idz == 63){
			// 	printf("SZ_compress_float_3D_MDQ_RA_block layer row %d done, cur_data_pos: %ld\n", i-1, cur_data_pos - data);
			// 	fflush(stdout);
			// }
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = last_row_pos[0] + cur_data_pos[- dim0_offset] - last_row_pos[- dim0_offset];
			type_ = type[index];
			if(type_ == 1){
				cur_data_pos[0] = mean;
			}
			else if (type_ != 0){
				cur_data_pos[0] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				cur_data_pos[0] = unpredictable_data[unpredictable_count ++];
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
//				if(k==63&&i==43&&j==27)
//					printf("i=%d\n", i);
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				pred3D = cur_data_pos[j-1] + last_row_pos[j]+ cur_data_pos[j - dim0_offset] - last_row_pos[j-1] - last_row_pos[j - dim0_offset] - cur_data_pos[j-1 - dim0_offset] + last_row_pos[j-1 - dim0_offset];
				type_ = type[index];
				if(type_ == 1){
					cur_data_pos[j] = mean;
				}
				else if (type_ != 0){
					cur_data_pos[j] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
				}
				else{
					cur_data_pos[j] = unpredictable_data[unpredictable_count ++];
				}
			}
			last_row_pos = cur_data_pos;
			cur_data_pos += dim1_offset;
		}
		cur_data_pos += dim0_offset - r2 * dim1_offset;
	}

	return unpredictable_count;
}

void decompressDataSeries_float_1D_RA(float** data, size_t r1, unsigned char * comp_data){

	size_t num_elements = r1;
	*data = (float*)malloc(sizeof(float)*num_elements);

	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;
	
	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	size_t num_x;
	size_t early_blockcount_x, late_blockcount_x;
	size_t split_index_x;

	SZ_COMPUTE_1D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	SZ_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);

	size_t max_num_block_elements = early_blockcount_x;
	size_t num_blocks = num_x;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2;
	// intvRadius = (int)((tdps->intervals - 1)/ 2);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;
	unsigned short * block_pos = (unsigned short *) comp_data_pos;
	// skip block index here
	comp_data_pos += num_blocks * sizeof(unsigned short);
	unsigned short * unpred_count_pos = (unsigned short *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(unsigned short);
	float * mean_pos = (float *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(float);

	unsigned int unpredictable_count;
	int * type = (int *) malloc(max_num_block_elements * sizeof(int));
	// int unpred_data_max_size = ((int)(num_block_elements * 0.2) + 1);
	size_t unpred_data_max_size = max_num_block_elements;
	float * unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float));
	float mean;
	unsigned char * tmp;
	unsigned int unpredictableEncodeSize;
	size_t index = 0;
	float * data_pos = *data;
	// printf("decompress offset to start: %ld\n", comp_data_pos - tdps->data);
	// fflush(stdout);
	// exit(0);

	size_t offset_x = 0;
	size_t type_offset = 0;
	size_t current_blockcount_x;
	size_t cur_unpred_count;
	for(size_t i=0; i<num_blocks; i++){
		offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
		data_pos = *data + offset_x;
		current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;

		tmp = comp_data_pos;
		mean = *mean_pos;
		unpredictable_count = *(unpred_count_pos);
		if(unpredictable_count > 0){
			unpredictableEncodeSize = unpredictable_count * sizeof(float);
			memcpy(unpredictable_data, tmp, unpredictableEncodeSize);
			tmp += unpredictableEncodeSize;
		}
		// caculate real block elements
		decode(tmp, current_blockcount_x, root, type);

		cur_unpred_count = decompressDataSeries_float_1D_RA_block_1D_pred(data_pos, mean, num_elements, current_blockcount_x, realPrecision, type, unpredictable_data);
		if(cur_unpred_count != unpredictable_count){
			printf("Check bugs, unpredictable_count is not the same: %d %d\n", unpredictable_count, cur_unpred_count);
			printf("Current index: %d\n\n", i);
			for(size_t i=0; i<current_blockcount_x; i++){
				printf("%d ", type[i]);
			}
			printf("\n");
			exit(0);
		}

		comp_data_pos += *block_pos;
		block_pos ++;
		unpred_count_pos ++;
		mean_pos ++;

		// printf("block comp done, data_offset from %ld to %ld: diff %ld\n", *data, data_pos, data_pos - *data);
		// fflush(stdout);

	}
	free(type);
	free(unpredictable_data);
}

void decompressDataSeries_float_3D_nonblocked_ori(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// calculate block dims
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t num_elements = r1 * r2 * r3;
	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2;
	// intvRadius = (int)((tdps->intervals - 1)/ 2);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;

	size_t unpred_count = *((unsigned int *)comp_data_pos);
	comp_data_pos += 4;
	float * unpredictable_data = (float *) comp_data_pos;
	comp_data_pos += unpred_count * sizeof(float);
	
	float mean = *((float *) comp_data_pos);
	comp_data_pos += 4;

	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	int * type = (int *) malloc(num_elements * sizeof(int));
	float * data_pos = *data;
	decode(comp_data_pos, num_elements, root, type);
	size_t cur_unpred_count = decompressDataSeries_float_3D_RA_block(data_pos, mean, r1, r2, r3, r1, r2, r3, realPrecision, type, unpredictable_data);
	// size_t cur_unpred_count = decompressDataSeries_float_3D_RA_block_2_layers(data_pos, mean, r1, r2, r3, r1, r2, r3, realPrecision, type, unpredictable_data);

	if(cur_unpred_count != unpred_count){
		printf("Check bugs, unpredictable_count is not the same: %d %d\n", unpred_count, cur_unpred_count);
		for(size_t i=0; i<512; i++){
			printf("%d ", type[i]);
		}
		printf("\n");
		exit(0);
	}
	free(type);

}

void decompressDataSeries_float_3D_nonblocked_adaptive(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// calculate block dims
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t num_elements = r1 * r2 * r3;
	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2;
	// intvRadius = (int)((tdps->intervals - 1)/ 2);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;

	size_t unpred_count = *((unsigned int *)comp_data_pos);
	comp_data_pos += 4;
	float * unpredictable_data = (float *) comp_data_pos;
	comp_data_pos += unpred_count * sizeof(float);
	
	float mean = *((float *) comp_data_pos);
	comp_data_pos += 4;

	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	int * type = (int *) malloc(num_elements * sizeof(int));
	float * data_pos = *data;
	decode(comp_data_pos, num_elements, root, type);
	// size_t cur_unpred_count = decompressDataSeries_float_3D_RA_block_3D_pred(data_pos, mean, r1, r2, r3, r1, r2, r3, realPrecision, type, unpredictable_data);
	size_t cur_unpred_count = decompressDataSeries_float_3D_RA_block_adaptive(data_pos, mean, r1, r2, r3, r1, r2, r3, realPrecision, type, unpredictable_data);

	if(cur_unpred_count != unpred_count){
		printf("Check bugs, unpredictable_count is not the same: %d %d\n", unpred_count, cur_unpred_count);
		for(size_t i=0; i<512; i++){
			printf("%d ", type[i]);
		}
		printf("\n");
		exit(0);
	}
	free(type);

}


void decompressDataSeries_float_3D_nonblocked(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// calculate block dims
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t num_elements = r1 * r2 * r3;
	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2;
	// intvRadius = (int)((tdps->intervals - 1)/ 2);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;

	size_t unpred_count = *((unsigned int *)comp_data_pos);
	comp_data_pos += 4;
	float * unpredictable_data = (float *) comp_data_pos;
	comp_data_pos += unpred_count * sizeof(float);
	
	float mean = *((float *) comp_data_pos);
	comp_data_pos += 4;

	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	int * type = (int *) malloc(num_elements * sizeof(int));
	float * data_pos = *data;
	decode(comp_data_pos, num_elements, root, type);
	size_t cur_unpred_count = decompressDataSeries_float_3D_RA_block_3D_pred(data_pos, mean, r1, r2, r3, r1, r2, r3, realPrecision, type, unpredictable_data);

	if(cur_unpred_count != unpred_count){
		printf("Check bugs, unpredictable_count is not the same: %d %d\n", unpred_count, cur_unpred_count);
		for(size_t i=0; i<512; i++){
			printf("%d ", type[i]);
		}
		printf("\n");
		exit(0);
	}
	free(type);

}

void decompressDataSeries_float_1D_RA_all_by_regression(float** data, size_t r1, unsigned char* comp_data){
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t num_elements = r1;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	// calculate block dims
	size_t num_x = r1 / block_size;
	size_t num_blocks = num_x;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2;
	// intvRadius = (int)((tdps->intervals - 1)/ 2);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;
	size_t unpredictable_count;
	memcpy(&unpredictable_count, comp_data_pos, sizeof(size_t));
	comp_data_pos += sizeof(size_t);
	float * unpred_data = (float *) comp_data_pos;
	comp_data_pos += unpredictable_count * sizeof(float);
	float * reg_params = (float *) comp_data_pos;
	comp_data_pos += num_blocks * 2 * sizeof(float);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	decode(comp_data_pos, num_elements, root, result_type);

	// printf("decompress offset to start: %ld\n", comp_data_pos - tdps->data);
	// fflush(stdout);
	int type_;
	float pred;
	int * type = result_type;
	float * unpredictable_data = unpred_data;
	size_t dec_unpred_count = 0;
	float * data_pos = *data;
	for(size_t i=0; i<num_x; i++){
		for(size_t ii=0; ii<block_size; ii++){
			type_ = type[0];
			if (type_ != 0){
				pred = reg_params[0] * ii + reg_params[num_blocks];
				*data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
			}
			else{
				*data_pos = unpredictable_data[dec_unpred_count ++];
			}
			data_pos ++;
			type ++;
		}
	}
	// residue
	for(size_t i=0; i<r1 - num_x*block_size; i++){
		*data_pos = unpredictable_data[dec_unpred_count ++];
		data_pos ++;
	}

	free(result_type);
}

size_t decompressDataSeries_float_3D_RA_block_all_by_regression(float * block_ori_data, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, double realPrecision, float * reg_params, int * type, float * unpredictable_data){

	float * data_pos;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;
	float curData;
	float pred;
	int type_;
	size_t index = 0;
	size_t unpredictable_count = 0;
	data_pos = block_ori_data;
	for(size_t i=0; i<block_dim_0; i++){
		for(size_t j=0; j<block_dim_1; j++){
			for(size_t k=0; k<block_dim_2; k++){
				type_ = type[index];
				if (type_ != 0){
					pred = reg_params[0] * i + reg_params[1] * j + reg_params[2] * k + reg_params[3];
					*data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
				}
				else{
					*data_pos = unpredictable_data[unpredictable_count ++];
				}
				// if(*data_pos - *(tmp_data + (data_pos - tmp_dec_data)) < -1.3046569824 || *data_pos - *(tmp_data + (data_pos - tmp_dec_data)) > 1.3046569824){
					// printf("DEC REG PPPPPPP\n");
				// }
				index ++;	
				data_pos ++;
			}
			data_pos += dim1_offset - block_dim_2;
		}
		data_pos += dim0_offset - block_dim_1 * dim1_offset;
	}
	return unpredictable_count;
}

void decompressDataSeries_float_3D_RA_all_by_regression(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	size_t num_elements = r1 * r2 * r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	// calculate block dims
	size_t num_x, num_y, num_z;
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	SZ_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	SZ_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	SZ_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2;
	// intvRadius = (int)((tdps->intervals - 1)/ 2);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;
	unsigned short * unpred_count = (unsigned short *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(unsigned short);
	size_t total_unpred = 0;
	for(size_t i=0; i<num_blocks; i++){
		total_unpred += unpred_count[i];
	}
	float * unpred_data = (float *) comp_data_pos;
	comp_data_pos += total_unpred * sizeof(float);
	float * reg_params = (float *) comp_data_pos;
	comp_data_pos += num_blocks * 4 * sizeof(float);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	decode(comp_data_pos, num_elements, root, result_type);
	int * type;
	size_t unpredictable_count;
	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	size_t index = 0;
	float * data_pos = *data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t cur_unpred_count;
	size_t type_offset = 0;
	// printf("decompress offset to start: %ld\n", comp_data_pos - tdps->data);
	// fflush(stdout);
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
				data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

				type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
				type = result_type + type_offset;
				size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
				unpredictable_count = * unpred_count;
				cur_unpred_count = decompressDataSeries_float_3D_RA_block_all_by_regression(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, reg_params, type, unpred_data);

				if(cur_unpred_count != unpredictable_count){
					printf("Check bugs, unpredictable_count is not the same: %d %d\n", unpredictable_count, cur_unpred_count);
					printf("Current index: %d %d %d\n\n", i, j, k);
					for(int i=0; i<current_block_elements; i++){
						printf("%d ", type[i]);
					}
					printf("\n");
					exit(0);
				}

				unpred_data += unpredictable_count;
				unpred_count ++;
				reg_params += 4;

				// printf("block comp done, data_offset from %ld to %ld: diff %ld\n", *data, data_pos, data_pos - *data);
				// fflush(stdout);

			}
		}
	}
	free(result_type);
}

size_t decompressDataSeries_float_3D_blocked_nonblock_pred(float * data, size_t dim_0, size_t dim_1, size_t dim_2, size_t block_dim_0, size_t block_dim_1, size_t block_dim_2, size_t idx, size_t idy, size_t idz, double realPrecision, int * type, float * unpredictable_data){

	// float * tmp_tmp_pos = tmp_pos;

	float * data_pos;
	float pred;
	size_t dim0_offset = dim_1 * dim_2;
	size_t dim1_offset = dim_2;

	data_pos = data;
	size_t index = 0;
	int type_;
	// d111 is current data
	size_t unpredictable_count = 0;
	float d000, d001, d010, d011, d100, d101, d110;
	for(size_t i=0; i<block_dim_0; i++){
		for(size_t j=0; j<block_dim_1; j++){
			for(size_t k=0; k<block_dim_2; k++){
				d000 = d001 = d010 = d011 = d100 = d101 = d110 = 1;
				if(idx == 0 && i == 0){
					d000 = d001 = d010 = d011 = 0;
				}
				if(idy == 0 && j == 0){
					d000 = d001 = d100 = d101 = 0;
				}
				if(idz == 0 && k == 0){
					d000 = d010 = d100 = d110 = 0;
				}
				if(d000){
					d000 = data_pos[- dim0_offset - dim1_offset - 1];
				}
				if(d001){
					d001 = data_pos[- dim0_offset - dim1_offset];
				}
				if(d010){
					d010 = data_pos[- dim0_offset - 1];
				}
				if(d011){
					d011 = data_pos[- dim0_offset];
				}
				if(d100){
					d100 = data_pos[- dim1_offset - 1];
				}
				if(d101){
					d101 = data_pos[- dim1_offset];
				}
				if(d110){
					d110 = data_pos[- 1];
				}
				type_ = type[index];
				if (type_ != 0){
					pred = d110 + d101 + d011 - d100 - d010 - d001 + d000;
					*data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
				}
				else{
					*data_pos = unpredictable_data[unpredictable_count ++];
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
void decompressDataSeries_float_2D_nonblocked_with_blocked_regression(float** data, size_t r1, size_t r2, unsigned char* comp_data){
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);
	size_t dim0_offset = r2;
	size_t num_elements = r1 * r2;

	*data = (float*)malloc(sizeof(float)*num_elements);
	// tmp_dec_data = *data;

	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	// calculate block dims
	size_t num_x, num_y;
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);

	size_t split_index_x, split_index_y;
	size_t early_blockcount_x, early_blockcount_y;
	size_t late_blockcount_x, late_blockcount_y;
	SZ_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	SZ_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y;
	size_t num_blocks = num_x * num_y;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);
	comp_data_pos += 4 + tree_size;

	float mean;
	unsigned char use_mean;
	memcpy(&use_mean, comp_data_pos, sizeof(unsigned char));
	comp_data_pos += 1;
	memcpy(&mean, comp_data_pos, sizeof(float));
	comp_data_pos += 4;
	size_t reg_count = 0;
	// unsigned char * indicator = (unsigned char *) comp_data_pos;
	// comp_data_pos += num_blocks * sizeof(unsigned char);
	unsigned char * indicator;
	size_t indicator_bitlength = (num_blocks - 1)/8 + 1;
	convertByteArray2IntArray_fast_1b(num_blocks, comp_data_pos, indicator_bitlength, &indicator);
	comp_data_pos += indicator_bitlength;
	for(size_t i=0; i<num_blocks; i++){
		if(!indicator[i]) reg_count ++;
	}
	printf("reg_count: %ld\n", reg_count);
	// float * reg_params_buf = (float *) malloc(reg_count * 4 * sizeof(float));
	float* dec_a = NULL, *dec_b = NULL, *dec_c = NULL;
	if(reg_count > 0){
		float * medians = (float *) comp_data_pos;
		float medianValue_a = *medians;
		float medianValue_b = *(medians + 1);
		float medianValue_c = *(medians + 2);
		comp_data_pos += 3 * sizeof(float);

		int * reqLength = (int *) comp_data_pos;
		int reqLength_a = *reqLength;
		int reqLength_b = *(reqLength + 1);
		int reqLength_c = *(reqLength + 2);
		comp_data_pos += 3 * sizeof(int);
		
		//reconstruct leading_number array in the form of meaningful integers....
		size_t leadNumArray_size = (reg_count-1)/4+1;
		
		unsigned char* leadNum_a = NULL, *leadNum_b = NULL, *leadNum_c = NULL;
		
		unsigned char* leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_a);
		comp_data_pos += leadNumArray_size;
		
		leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_b);
		comp_data_pos += leadNumArray_size;
		
		leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_c);	
		comp_data_pos += leadNumArray_size;
				
		//reconstruct mid bytes...
		size_t * mid_byte_size_a = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_a = comp_data_pos;
		comp_data_pos += *mid_byte_size_a;
		
		size_t * mid_byte_size_b = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_b = comp_data_pos;	
		comp_data_pos += *mid_byte_size_b;	
		
		size_t * mid_byte_size_c = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_c = comp_data_pos;	
		comp_data_pos += *mid_byte_size_c;			
		
		//reconstruct the residualMidBits		
		size_t * resiMidBites_a_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_a = comp_data_pos;
		comp_data_pos += *resiMidBites_a_size;

		size_t * resiMidBites_b_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_b = comp_data_pos;
		comp_data_pos += *resiMidBites_b_size;
		
		size_t * resiMidBites_c_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_c = comp_data_pos;
		comp_data_pos += *resiMidBites_c_size;
				
		//perform the decompression using the reconstructed leadNum, exactMidBytes and residualMidBits....
		decompressExactDataArray_float(leadNum_a, exactMidBytes_a, residualMidBits_a, reg_count, reqLength_a, medianValue_a, &dec_a);
		decompressExactDataArray_float(leadNum_b, exactMidBytes_b, residualMidBits_b, reg_count, reqLength_b, medianValue_b, &dec_b);
		decompressExactDataArray_float(leadNum_c, exactMidBytes_c, residualMidBits_c, reg_count, reqLength_c, medianValue_c, &dec_c);
		free(leadNum_a);
		free(leadNum_b);
		free(leadNum_c);
		
	}
	size_t total_unpred;
	memcpy(&total_unpred, comp_data_pos, sizeof(size_t));
	comp_data_pos += sizeof(size_t);
	float * unpred_data = (float *) comp_data_pos;
	comp_data_pos += total_unpred * sizeof(float);

	int * result_type = (int *) malloc(num_elements * sizeof(int));
	decode(comp_data_pos, num_elements, root, result_type);
	int * type;
	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	size_t index = 0;
	float * data_pos = *data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t cur_unpred_count;
	// size_t type_offset = 0;
	// size_t decomp_unpred = 0;
	// printf("decompress offset to start: %ld\n", comp_data_pos - tdps->data);
	// fflush(stdout);
	float * dec_a_pos = dec_a;
	float * dec_b_pos = dec_b;
	float * dec_c_pos = dec_c;
	unsigned char * indicator_pos = indicator;
	if(use_mean){
		type = result_type;
		for(size_t i=0; i<num_x; i++){
			for(size_t j=0; j<num_y; j++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				data_pos = *data + offset_x * dim0_offset + offset_y;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;

				size_t current_block_elements = current_blockcount_x * current_blockcount_y;
				if(*indicator_pos){
					// decompress by SZ
					// cur_unpred_count = decompressDataSeries_float_3D_blocked_nonblock_pred(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, i, j, k, realPrecision, type, unpred_data);
					float * block_data_pos = data_pos;
					float pred;
					size_t index = 0;
					int type_;
					// d11 is current data
					size_t unpredictable_count = 0;
					float d00, d01, d10;
					for(size_t ii=0; ii<current_blockcount_x; ii++){
						for(size_t jj=0; jj<current_blockcount_y; jj++){
							type_ = type[index];
							if(type_ == intvRadius){
								*block_data_pos = mean;
							}
							else if(type_ == 0){
								*block_data_pos = unpred_data[unpredictable_count ++];
							}
							else{
								d00 = d01 = d10 = 1;
								if(i == 0 && ii == 0){
									d00 = d01 = 0;
								}
								if(j == 0 && jj == 0){
									d00 = d10 = 0;
								}
								if(d00){
									d00 = block_data_pos[- dim0_offset - 1];
								}
								if(d01){
									d01 = block_data_pos[- dim0_offset];
								}
								if(d10){
									d10 = block_data_pos[- 1];
								}
								if(type_ < intvRadius) type_ += 1;
								pred = d10 + d01 - d00;
								*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
							}
							index ++;
							block_data_pos ++;
						}
						block_data_pos += dim0_offset - current_blockcount_y;
					}
					cur_unpred_count = unpredictable_count;
				}
				else{
					// decompress by regression
					// cur_unpred_count = decompressDataSeries_float_3D_RA_block_all_by_regression(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, reg_params, type, unpred_data);
					{
						float * block_data_pos = data_pos;
						float curData;
						float pred;
						int type_;
						size_t index = 0;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								type_ = type[index];
								if (type_ != 0){
									pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0];
									*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
								}
								else{
									*block_data_pos = unpred_data[unpredictable_count ++];
								}
								// if(*block_data_pos - *(tmp_data + (block_data_pos - tmp_dec_data)) < -1.3046569824 || *block_data_pos - *(tmp_data + (block_data_pos - tmp_dec_data)) > 1.3046569824){
									// printf("DEC REG PPPPPPP\n");
								// }
								index ++;	
								block_data_pos ++;
							}
							block_data_pos += dim0_offset - current_blockcount_y;
						}
						cur_unpred_count = unpredictable_count;
					}
					// reg_params += 4;
					dec_a_pos ++, dec_b_pos ++, dec_c_pos ++;
				}

				type += current_block_elements;
				indicator_pos ++;
				unpred_data += cur_unpred_count;
				// decomp_unpred += cur_unpred_count;
				// printf("block comp done, data_offset from %ld to %ld: diff %ld\n", *data, data_pos, data_pos - *data);
				// fflush(stdout);
			}
		}
	}
	else{
		type = result_type;
		for(size_t i=0; i<num_x; i++){
			for(size_t j=0; j<num_y; j++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				data_pos = *data + offset_x * dim0_offset + offset_y;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;

				size_t current_block_elements = current_blockcount_x * current_blockcount_y;
				if(*indicator_pos){
					// decompress by SZ
					// cur_unpred_count = decompressDataSeries_float_3D_blocked_nonblock_pred(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, i, j, k, realPrecision, type, unpred_data);
					float * block_data_pos = data_pos;
					float pred;
					size_t index = 0;
					int type_;
					// d11 is current data
					size_t unpredictable_count = 0;
					float d00, d01, d10;
					for(size_t ii=0; ii<current_blockcount_x; ii++){
						for(size_t jj=0; jj<current_blockcount_y; jj++){
							type_ = type[index];
							if(type_ == 0){
								*block_data_pos = unpred_data[unpredictable_count ++];
							}
							else{
								d00 = d01 = d10 = 1;
								if(i == 0 && ii == 0){
									d00 = d01 = 0;
								}
								if(j == 0 && jj == 0){
									d00 = d10 = 0;
								}
								if(d00){
									d00 = block_data_pos[- dim0_offset - 1];
								}
								if(d01){
									d01 = block_data_pos[- dim0_offset];
								}
								if(d10){
									d10 = block_data_pos[- 1];
								}
								pred = d10 + d01 - d00;
								*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
							}
							index ++;
							block_data_pos ++;
						}
						block_data_pos += dim0_offset - current_blockcount_y;
					}
					cur_unpred_count = unpredictable_count;
				}
				else{
					// decompress by regression
					// cur_unpred_count = decompressDataSeries_float_3D_RA_block_all_by_regression(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, reg_params, type, unpred_data);
					{
						float * block_data_pos = data_pos;
						float curData;
						float pred;
						int type_;
						size_t index = 0;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								type_ = type[index];
								if (type_ != 0){
									pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0];
									*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
								}
								else{
									*block_data_pos = unpred_data[unpredictable_count ++];
								}
								// if(*block_data_pos - *(tmp_data + (block_data_pos - tmp_dec_data)) < -1.3046569824 || *block_data_pos - *(tmp_data + (block_data_pos - tmp_dec_data)) > 1.3046569824){
									// printf("DEC REG PPPPPPP\n");
								// }
								index ++;	
								block_data_pos ++;
							}
							block_data_pos += dim0_offset - current_blockcount_y;
						}
						cur_unpred_count = unpredictable_count;
					}
					// reg_params += 4;
					dec_a_pos ++, dec_b_pos ++, dec_c_pos ++;
				}

				type += current_block_elements;
				indicator_pos ++;
				unpred_data += cur_unpred_count;
				// decomp_unpred += cur_unpred_count;
				// printf("block comp done, data_offset from %ld to %ld: diff %ld\n", *data, data_pos, data_pos - *data);
				// fflush(stdout);
			}
		}
	}
	if(NULL != dec_a) free(dec_a);
	if(NULL != dec_b) free(dec_b);
	if(NULL != dec_c) free(dec_c);

	free(indicator);
	free(result_type);
}
void decompressDataSeries_float_3D_nonblocked_with_blocked_regression(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);
	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	size_t num_elements = r1 * r2 * r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	// tmp_dec_data = *data;

	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	// calculate block dims
	size_t num_x, num_y, num_z;
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	SZ_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	SZ_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	SZ_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);
	comp_data_pos += 4 + tree_size;

	float mean;
	unsigned char use_mean;
	memcpy(&use_mean, comp_data_pos, sizeof(unsigned char));
	comp_data_pos += 1;
	memcpy(&mean, comp_data_pos, sizeof(float));
	comp_data_pos += sizeof(float);
	size_t reg_count = 0;
	// unsigned char * indicator = (unsigned char *) comp_data_pos;
	// comp_data_pos += num_blocks * sizeof(unsigned char);
	unsigned char * indicator;
	size_t indicator_bitlength = (num_blocks - 1)/8 + 1;
	convertByteArray2IntArray_fast_1b(num_blocks, comp_data_pos, indicator_bitlength, &indicator);
	comp_data_pos += indicator_bitlength;
	for(size_t i=0; i<num_blocks; i++){
		if(!indicator[i]) reg_count ++;
	}
	// printf("reg_count: %ld\n", reg_count);
	// float * reg_params_buf = (float *) malloc(reg_count * 4 * sizeof(float));
	float* dec_a = NULL, *dec_b = NULL, *dec_c = NULL, *dec_d = NULL;
	if(reg_count > 0){
		float * medians = (float *) comp_data_pos;
		float medianValue_a = *medians;
		float medianValue_b = *(medians + 1);
		float medianValue_c = *(medians + 2);
		float medianValue_d = *(medians + 3);
		comp_data_pos += 4 * sizeof(float);

		int * reqLength = (int *) comp_data_pos;
		int reqLength_a = *reqLength;
		int reqLength_b = *(reqLength + 1);
		int reqLength_c = *(reqLength + 2);
		int reqLength_d = *(reqLength + 3);
		comp_data_pos += 4 * sizeof(int);
		
		//reconstruct leading_number array in the form of meaningful integers....
		size_t leadNumArray_size = (reg_count-1)/4+1;
		
		unsigned char* leadNum_a = NULL, *leadNum_b = NULL, *leadNum_c = NULL, *leadNum_d = NULL;
		
		unsigned char* leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_a);
		comp_data_pos += leadNumArray_size;
		
		leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_b);
		comp_data_pos += leadNumArray_size;
		
		leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_c);	
		comp_data_pos += leadNumArray_size;
		
		leadNumArray = comp_data_pos;
		convertByteArray2IntArray_fast_2b(reg_count, leadNumArray, leadNumArray_size, &leadNum_d);			
		comp_data_pos += leadNumArray_size;
		
		//reconstruct mid bytes...
		size_t * mid_byte_size_a = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_a = comp_data_pos;
		comp_data_pos += *mid_byte_size_a;
		
		size_t * mid_byte_size_b = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_b = comp_data_pos;	
		comp_data_pos += *mid_byte_size_b;	
		
		size_t * mid_byte_size_c = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_c = comp_data_pos;	
		comp_data_pos += *mid_byte_size_c;			
		
		size_t * mid_byte_size_d = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* exactMidBytes_d = comp_data_pos;	
		comp_data_pos += *mid_byte_size_d;			

		//reconstruct the residualMidBits		
		size_t * resiMidBites_a_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_a = comp_data_pos;
		comp_data_pos += *resiMidBites_a_size;

		size_t * resiMidBites_b_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_b = comp_data_pos;
		comp_data_pos += *resiMidBites_b_size;
		
		size_t * resiMidBites_c_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_c = comp_data_pos;
		comp_data_pos += *resiMidBites_c_size;
		
		size_t * resiMidBites_d_size = (size_t *) comp_data_pos;
		comp_data_pos += sizeof(size_t);
		unsigned char* residualMidBits_d = comp_data_pos;
		comp_data_pos += *resiMidBites_d_size;				
		
		//perform the decompression using the reconstructed leadNum, exactMidBytes and residualMidBits....
		decompressExactDataArray_float(leadNum_a, exactMidBytes_a, residualMidBits_a, reg_count, reqLength_a, medianValue_a, &dec_a);
		decompressExactDataArray_float(leadNum_b, exactMidBytes_b, residualMidBits_b, reg_count, reqLength_b, medianValue_b, &dec_b);
		decompressExactDataArray_float(leadNum_c, exactMidBytes_c, residualMidBits_c, reg_count, reqLength_c, medianValue_c, &dec_c);
		decompressExactDataArray_float(leadNum_d, exactMidBytes_d, residualMidBits_d, reg_count, reqLength_d, medianValue_d, &dec_d);
		// for(size_t i=0; i<reg_count; i++){
		// 	reg_params_buf[4*i] = dec_a[i];
		// 	reg_params_buf[4*i + 1] = dec_b[i];
		// 	reg_params_buf[4*i + 2] = dec_c[i];
		// 	reg_params_buf[4*i + 3] = dec_d[i];
		// }
		// free(dec_a);
		// free(dec_b);
		// free(dec_c);
		// free(dec_d);
		free(leadNum_a);
		free(leadNum_b);
		free(leadNum_c);
		free(leadNum_d);				
		
	}
	// float * reg_params = reg_params_buf;
	// float * reg_params = (float *) comp_data_pos;
	// comp_data_pos += reg_count * 4 * sizeof(float);
	// reorder reg_params

	// for(size_t i=0; i<reg_count; i++){
	// 	for(size_t j=0; j<4; j++){
	// 		reg_params_buf[4*i + j] = reg_params[j*reg_count + i];
	// 	}
	// }
	// reg_params = reg_params_buf;
	// {
	// 	int status;
	// 	writeFloatData_inBytes(reg_params, reg_count * 4, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/params.dat", &status);
	// }
	size_t total_unpred;
	memcpy(&total_unpred, comp_data_pos, sizeof(size_t));
	comp_data_pos += sizeof(size_t);
	float * unpred_data = (float *) comp_data_pos;
	comp_data_pos += total_unpred * sizeof(float);
	// printf("total_unpred: %ld\n", total_unpred);
	// for(size_t i=0; i<total_unpred; i++){
	// 	printf("%.2f ", unpred_data[i]);
	// }
	// printf("\n");
	// getchar();

	int * result_type = (int *) malloc(num_elements * sizeof(int));
	decode(comp_data_pos, num_elements, root, result_type);
	int * type;
	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	size_t index = 0;
	float * data_pos = *data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t cur_unpred_count;
	// size_t type_offset = 0;
	// size_t decomp_unpred = 0;
	// printf("decompress offset to start: %ld\n", comp_data_pos - tdps->data);
	// fflush(stdout);
	float * dec_a_pos = dec_a;
	float * dec_b_pos = dec_b;
	float * dec_c_pos = dec_c;
	float * dec_d_pos = dec_d;
	unsigned char * indicator_pos = indicator;
	if(use_mean){
		type = result_type;
		// for(size_t i=0; i<num_x; i++){
		// 	for(size_t j=0; j<num_y; j++){
		// 		for(size_t k=0; k<num_z; k++){
		// 			offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
		// 			offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
		// 			offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
		// 			data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

		// 			current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
		// 			current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
		// 			current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

		// 			// type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
		// 			// type = result_type + type_offset;
		// 			size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
		// 			// index = i * num_y * num_z + j * num_z + k;

		// 			// printf("i j k: %ld %ld %ld\toffset: %ld %ld %ld\tindicator: %ld\n", i, j, k, offset_x, offset_y, offset_z, indicator[index]);
		// 			if(*indicator_pos){
		// 				// decompress by SZ
		// 				// cur_unpred_count = decompressDataSeries_float_3D_blocked_nonblock_pred(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, i, j, k, realPrecision, type, unpred_data);
		// 				float * block_data_pos = data_pos;
		// 				float pred;
		// 				size_t index = 0;
		// 				int type_;
		// 				// d111 is current data
		// 				size_t unpredictable_count = 0;
		// 				float d000, d001, d010, d011, d100, d101, d110;
		// 				for(size_t ii=0; ii<current_blockcount_x; ii++){
		// 					for(size_t jj=0; jj<current_blockcount_y; jj++){
		// 						for(size_t kk=0; kk<current_blockcount_z; kk++){
		// 							type_ = type[index];
		// 							if(type_ == intvRadius){
		// 								*block_data_pos = mean;
		// 							}
		// 							else if(type_ == 0){
		// 								*block_data_pos = unpred_data[unpredictable_count ++];
		// 							}
		// 							else{
		// 								d000 = d001 = d010 = d011 = d100 = d101 = d110 = 1;
		// 								if(i == 0 && ii == 0){
		// 									d000 = d001 = d010 = d011 = 0;
		// 								}
		// 								if(j == 0 && jj == 0){
		// 									d000 = d001 = d100 = d101 = 0;
		// 								}
		// 								if(k == 0 && kk == 0){
		// 									d000 = d010 = d100 = d110 = 0;
		// 								}
		// 								if(d000){
		// 									d000 = block_data_pos[- dim0_offset - dim1_offset - 1];
		// 								}
		// 								if(d001){
		// 									d001 = block_data_pos[- dim0_offset - dim1_offset];
		// 								}
		// 								if(d010){
		// 									d010 = block_data_pos[- dim0_offset - 1];
		// 								}
		// 								if(d011){
		// 									d011 = block_data_pos[- dim0_offset];
		// 								}
		// 								if(d100){
		// 									d100 = block_data_pos[- dim1_offset - 1];
		// 								}
		// 								if(d101){
		// 									d101 = block_data_pos[- dim1_offset];
		// 								}
		// 								if(d110){
		// 									d110 = block_data_pos[- 1];
		// 								}
		// 								if(type_ < intvRadius) type_ += 1;
		// 								pred = d110 + d101 + d011 - d100 - d010 - d001 + d000;
		// 								*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
		// 							}
		// 							index ++;
		// 							block_data_pos ++;
		// 						}
		// 						block_data_pos += dim1_offset - current_blockcount_z;
		// 					}
		// 					block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
		// 				}
		// 				cur_unpred_count = unpredictable_count;
		// 			}
		// 			else{
		// 				// decompress by regression
		// 				// cur_unpred_count = decompressDataSeries_float_3D_RA_block_all_by_regression(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, reg_params, type, unpred_data);
		// 				{
		// 					float * block_data_pos = data_pos;
		// 					float curData;
		// 					float pred;
		// 					int type_;
		// 					size_t index = 0;
		// 					size_t unpredictable_count = 0;
		// 					for(size_t ii=0; ii<current_blockcount_x; ii++){
		// 						for(size_t jj=0; jj<current_blockcount_y; jj++){
		// 							for(size_t kk=0; kk<current_blockcount_z; kk++){
		// 								type_ = type[index];
		// 								if (type_ != 0){
		// 									// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
		// 									pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
		// 									*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
		// 								}
		// 								else{
		// 									*block_data_pos = unpred_data[unpredictable_count ++];
		// 								}
		// 								// if(*block_data_pos - *(tmp_data + (block_data_pos - tmp_dec_data)) < -1.3046569824 || *block_data_pos - *(tmp_data + (block_data_pos - tmp_dec_data)) > 1.3046569824){
		// 									// printf("DEC REG PPPPPPP\n");
		// 								// }
		// 								index ++;	
		// 								block_data_pos ++;
		// 							}
		// 							block_data_pos += dim1_offset - current_blockcount_z;
		// 						}
		// 						block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
		// 					}
		// 					cur_unpred_count = unpredictable_count;
		// 				}
		// 				// reg_params += 4;
		// 				dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
		// 			}

		// 			type += current_block_elements;
		// 			indicator_pos ++;
		// 			unpred_data += cur_unpred_count;
		// 			// decomp_unpred += cur_unpred_count;
		// 			// printf("block comp done, data_offset from %ld to %ld: diff %ld\n", *data, data_pos, data_pos - *data);
		// 			// fflush(stdout);
		// 		}
		// 	}
		// }
		type = result_type;
		// i == 0
		{
			// j == 0
			{
				// k == 0
				{
					data_pos = *data;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = 0;
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;						
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim0_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				// i == 0 j == 0 k != 0
				for(size_t k=1; k<num_z; k++){
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_z;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					// type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
					// type = result_type + type_offset;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}// end j==0
			for(size_t j=1; j<num_y; j++){
				// k == 0
				{
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					data_pos = *data + offset_y * dim1_offset;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						// reg_params += 4;
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_y * dim1_offset + offset_z;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					// type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
					// type = result_type + type_offset;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}
		} // end i==0
		for(size_t i=1; i<num_x; i++){
			// j == 0
			{
				// k == 0
				{
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					data_pos = *data + offset_x * dim0_offset;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim0_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_x * dim0_offset + offset_z;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					// type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
					// type = result_type + type_offset;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						float d000, d001, d010, d011, d100, d101, d110;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}// end j = 0
			for(size_t j=1; j<num_y; j++){
				// k == 0
				{
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						// reg_params += 4;
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					// type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
					// type = result_type + type_offset;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						float d000, d001, d010, d011, d100, d101, d110;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == intvRadius){
										*block_data_pos = mean;
									}
									else if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										if(type_ < intvRadius) type_ += 1;
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}
		}
	}
	else{
		type = result_type;
		// i == 0
		{
			// j == 0
			{
				// k == 0
				{
					data_pos = *data;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = 0;
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;						
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim0_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				// i == 0 j == 0 k != 0
				for(size_t k=1; k<num_z; k++){
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_z;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					// type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
					// type = result_type + type_offset;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}// end j==0
			for(size_t j=1; j<num_y; j++){
				// k == 0
				{
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					data_pos = *data + offset_y * dim1_offset;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						// reg_params += 4;
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_y * dim1_offset + offset_z;

					current_blockcount_x = early_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					// type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
					// type = result_type + type_offset;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						// ii == 0
						{
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] - block_data_pos[- dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						for(size_t ii=1; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}
		} // end i==0
		for(size_t i=1; i<num_x; i++){
			// j == 0
			{
				// k == 0
				{
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					data_pos = *data + offset_x * dim0_offset;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim0_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_x * dim0_offset + offset_z;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = early_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					// type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
					// type = result_type + type_offset;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						float d000, d001, d010, d011, d100, d101, d110;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							// jj == 0
							{
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							for(size_t jj=1; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}// end j = 0
			for(size_t j=1; j<num_y; j++){
				// k == 0
				{
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = early_blockcount_z;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								{
									// kk == 0
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim0_offset - dim1_offset];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								for(size_t kk=1; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						// reg_params += 4;
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				} // end k == 0
				for(size_t k=1; k<num_z; k++){
					offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
					offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
					offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
					data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

					current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
					current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
					current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

					// type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
					// type = result_type + type_offset;
					size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
					if(*indicator_pos){
						// decompress by SZ
						float * block_data_pos = data_pos;
						float pred;
						size_t index = 0;
						int type_;
						size_t unpredictable_count = 0;
						float d000, d001, d010, d011, d100, d101, d110;
						for(size_t ii=0; ii<current_blockcount_x; ii++){
							for(size_t jj=0; jj<current_blockcount_y; jj++){
								for(size_t kk=0; kk<current_blockcount_z; kk++){
									type_ = type[index];
									if(type_ == 0){
										*block_data_pos = unpred_data[unpredictable_count ++];
									}
									else{
										pred = block_data_pos[- 1] + block_data_pos[- dim1_offset] + block_data_pos[- dim0_offset] - block_data_pos[- dim1_offset - 1] - block_data_pos[- dim0_offset - 1] - block_data_pos[- dim0_offset - dim1_offset] + block_data_pos[- dim0_offset - dim1_offset - 1];
										*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
									}
									index ++;
									block_data_pos ++;
								}
								block_data_pos += dim1_offset - current_blockcount_z;
							}
							block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
						}
						cur_unpred_count = unpredictable_count;
					}
					else{
						// decompress by regression
						{
							float * block_data_pos = data_pos;
							float curData;
							float pred;
							int type_;
							size_t index = 0;
							size_t unpredictable_count = 0;
							for(size_t ii=0; ii<current_blockcount_x; ii++){
								for(size_t jj=0; jj<current_blockcount_y; jj++){
									for(size_t kk=0; kk<current_blockcount_z; kk++){
										type_ = type[index];
										if (type_ != 0){
											// pred = reg_params[0] * ii + reg_params[1] * jj + reg_params[2] * kk + reg_params[3];
											pred = dec_a_pos[0] * ii + dec_b_pos[0] * jj + dec_c_pos[0] * kk + dec_d_pos[0];
											*block_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
										}
										else{
											*block_data_pos = unpred_data[unpredictable_count ++];
										}
										index ++;	
										block_data_pos ++;
									}
									block_data_pos += dim1_offset - current_blockcount_z;
								}
								block_data_pos += dim0_offset - current_blockcount_y * dim1_offset;
							}
							cur_unpred_count = unpredictable_count;
						}
						dec_a_pos ++, dec_b_pos ++, dec_c_pos ++, dec_d_pos ++;
					}

					indicator_pos ++;
					type += current_block_elements;
					unpred_data += cur_unpred_count;
				}
			}
		}
	}
	if(NULL != dec_a) free(dec_a);
	if(NULL != dec_b) free(dec_b);
	if(NULL != dec_c) free(dec_c);
	if(NULL != dec_d) free(dec_d);

	free(indicator);
	// free(reg_params_buf);
	free(result_type);
}

void decompressDataSeries_float_3D_RA_blocked_with_regression(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	size_t num_elements = r1 * r2 * r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	// calculate block dims
	size_t num_x, num_y, num_z;
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	SZ_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	SZ_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	SZ_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2;
	// intvRadius = (int)((tdps->intervals - 1)/ 2);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;
	unsigned char * indicator = (unsigned char *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(unsigned char);
	unsigned short * unpred_count = (unsigned short *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(unsigned short);
	size_t total_unpred = 0;
	for(size_t i=0; i<num_blocks; i++){
		total_unpred += unpred_count[i];
	}
	float * unpred_data = (float *) comp_data_pos;
	comp_data_pos += total_unpred * sizeof(float);
	float * reg_params = (float *) comp_data_pos;
	comp_data_pos += num_blocks * 4 * sizeof(float);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	decode(comp_data_pos, num_elements, root, result_type);
	int * type;
	size_t unpredictable_count;
	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	size_t index = 0;
	float * data_pos = *data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t cur_unpred_count;
	size_t type_offset = 0;
	// printf("decompress offset to start: %ld\n", comp_data_pos - tdps->data);
	// fflush(stdout);
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
				data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

				type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
				type = result_type + type_offset;
				size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
				unpredictable_count = * unpred_count;
				index = i * num_y * num_z + j * num_z + k;
				if(indicator[index] == 0)
					cur_unpred_count = decompressDataSeries_float_3D_RA_block_all_by_regression(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, reg_params, type, unpred_data);
				else
					cur_unpred_count = decompressDataSeries_float_3D_RA_block_no_mean(data_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, type, unpred_data);

				if(cur_unpred_count != unpredictable_count){
					printf("Check bugs, unpredictable_count is not the same: %d %d\n", unpredictable_count, cur_unpred_count);
					printf("Current index: %d %d %d\n\n", i, j, k);
					for(int i=0; i<current_block_elements; i++){
						printf("%d ", type[i]);
					}
					printf("\n");
					exit(0);
				}

				unpred_data += unpredictable_count;
				unpred_count ++;
				reg_params += 4;

				// printf("block comp done, data_offset from %ld to %ld: diff %ld\n", *data, data_pos, data_pos - *data);
				// fflush(stdout);

			}
		}
	}
	free(result_type);

}

void decompressDataSeries_float_3D_RA(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);
	double elapsed_time = 0.0;
	clock_t start, end;
	start = clock();

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	size_t num_elements = r1 * r2 * r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	// calculate block dims
	size_t num_x, num_y, num_z;
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	SZ_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	SZ_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	SZ_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2;
	// intvRadius = (int)((tdps->intervals - 1)/ 2);

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;
	unsigned int * block_pos = (unsigned int *) comp_data_pos;
	// skip block index here
	comp_data_pos += num_blocks * sizeof(unsigned int);
	unsigned short * unpred_count = (unsigned short *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(unsigned short);
	float * mean_pos = (float *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(float);

	size_t unpredictable_count;

	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	int * type = (int *) malloc(max_num_block_elements * sizeof(int));
	size_t unpred_data_max_size = max_num_block_elements;
	float * unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float));
	float mean;
	unsigned char * tmp;
	size_t unpredictableEncodeSize;
	size_t index = 0;
	float * data_pos = *data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t cur_unpred_count;
	// printf("decompress offset to start: %ld\n", comp_data_pos - tdps->data);
	// fflush(stdout);
	end = clock();
	elapsed_time = ((double)(end - start)) /CLOCKS_PER_SEC;
	printf("Read info time: %.4f\n", elapsed_time);
	start = clock();
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
				data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

				tmp = comp_data_pos;
				mean = *mean_pos;
				// mean = 0.004558146;
				// mean = *((float *) tmp);
				// tmp += 4;
				unpredictable_count = *(unpred_count);
				if(unpredictable_count > 0){
					unpredictableEncodeSize = unpredictable_count * sizeof(float);
					memcpy(unpredictable_data, tmp, unpredictableEncodeSize);
					tmp += unpredictableEncodeSize;
				}
				size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
				decode(tmp, current_block_elements, root, type);

				int cur_unpred_count = decompressDataSeries_float_3D_RA_block(data_pos, mean, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, type, unpredictable_data);
				// cur_unpred_count = decompressDataSeries_float_3D_RA_block_3D_pred(data_pos, mean, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, type, unpredictable_data);
				//int cur_unpred_count = decompressDataSeries_float_3D_RA_block_1D_pred(data_pos, mean, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, type, unpredictable_data);
				// if(i==0 && j==0 && k==1){
				// 	int status;
				// 	writeIntData_inBytes(type, current_block_elements, "/Users/LiangXin/github/SZ-develop/example/openmp/decomp001_type.dat", &status);
				// }

				if(cur_unpred_count != unpredictable_count){
					printf("Check bugs, unpredictable_count is not the same: %d %d\n", unpredictable_count, cur_unpred_count);
					printf("Current index: %d %d %d\n\n", i, j, k);
					size_t count = 0;
					for(int i=0; i<current_block_elements; i++){
						// printf("%d ", type[i]);
						if(type[i] == 0) count ++;
					}
					printf("%ld\n", count);
					exit(0);
				}

				comp_data_pos += *block_pos;
				block_pos ++;
				unpred_count ++;
				mean_pos ++;

				// printf("block comp done, data_offset from %ld to %ld: diff %ld\n", *data, data_pos, data_pos - *data);
				// fflush(stdout);

			}
		}
	}
	end = clock();
	elapsed_time = ((double)(end - start)) /CLOCKS_PER_SEC;
	printf("Decode & decompress elapsed time: %.4f\n", elapsed_time);

	free(type);
	free(unpredictable_data);

}

void decompressDataSeries_float_3D_nonblocked_multi_means(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// calculate block dims
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t num_elements = r1 * r2 * r3;
	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);
	comp_data_pos += 4 + tree_size;

	float dense_pos = *((float *) comp_data_pos);
	comp_data_pos += 4;
	unsigned mean_count = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	float * means = (float *) comp_data_pos;
	comp_data_pos += mean_count * sizeof(float);

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2*((mean_count + 1)/2);
	intvRadius = intvCapacity/2 + 2*((mean_count + 1)/2);

	printf("decompress dense_pos %.8f mean_count %d intervals %d\n", dense_pos, mean_count, intervals);
	printf("capacity %d radius %d\n", intvCapacity, intvRadius);

	size_t unpred_count = *((unsigned int *)comp_data_pos);
	comp_data_pos += 4;
	float * unpredictable_data = (float *) comp_data_pos;
	comp_data_pos += unpred_count * sizeof(float);
	
	// float mean = *((float *) comp_data_pos);
	// comp_data_pos += 4;

	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	int * type = (int *) malloc(num_elements * sizeof(int));
	float * data_pos = *data;
	decode(comp_data_pos, num_elements, root, type);
	// size_t cur_unpred_count = decompressDataSeries_float_3D_RA_block_3D_pred(data_pos, mean, r1, r2, r3, r1, r2, r3, realPrecision, type, unpredictable_data);
	size_t cur_unpred_count = decompressDataSeries_float_3D_RA_block_3D_pred_multi_means(data_pos, mean_count, means, dense_pos, r1, r2, r3, r1, r2, r3, realPrecision, type, unpredictable_data);
	if(cur_unpred_count != unpred_count){
		printf("Check bugs, unpredictable_count is not the same: %d %d\n", unpred_count, cur_unpred_count);
		for(size_t i=0; i<512; i++){
			printf("%d ", type[i]);
		}
		printf("\n");
		exit(0);
	}
	free(type);

}

void decompressDataSeries_float_3D_RA_multi_means(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	size_t num_elements = r1 * r2 * r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;
	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	// calculate block dims
	size_t num_x, num_y, num_z;
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r2, num_y, block_size);
	SZ_COMPUTE_3D_NUMBER_OF_BLOCKS(r3, num_z, block_size);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	SZ_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	SZ_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	SZ_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);

	comp_data_pos += 4 + tree_size;

	float dense_pos = *((float *) comp_data_pos);
	comp_data_pos += 4;
	unsigned mean_count = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	float * means = (float *) comp_data_pos;
	comp_data_pos += mean_count * sizeof(float);
	unsigned short * block_pos = (unsigned short *) comp_data_pos;
	// skip block index here
	comp_data_pos += num_blocks * sizeof(unsigned short);
	unsigned short * unpred_count = (unsigned short *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(unsigned short);

	size_t unpredictable_count;

	updateQuantizationInfo(intervals);
	intvCapacity = intervals - 2*((mean_count + 1)/2);
	intvRadius = intvCapacity/2 + 2*((mean_count + 1)/2);

	printf("decompress dense_pos %.8f mean_count %d intervals %d\n", dense_pos, mean_count, intervals);
	printf("capacity %d radius %d\n", intvCapacity, intvRadius);
	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	int * type = (int *) malloc(max_num_block_elements * sizeof(int));
	size_t unpred_data_max_size = max_num_block_elements;
	float * unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float));
	unsigned char * tmp;
	size_t unpredictableEncodeSize;
	size_t index = 0;
	float * data_pos = *data;
	size_t offset_x, offset_y, offset_z;
	size_t current_blockcount_x, current_blockcount_y, current_blockcount_z;
	size_t cur_unpred_count;
	// printf("decompress offset to start: %ld\n", comp_data_pos - tdps->data);
	// fflush(stdout);
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
				offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
				offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
				data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

				current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
				current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
				current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;

				tmp = comp_data_pos;
				unpredictable_count = *(unpred_count);
				if(unpredictable_count > 0){
					unpredictableEncodeSize = unpredictable_count * sizeof(float);
					memcpy(unpredictable_data, tmp, unpredictableEncodeSize);
					tmp += unpredictableEncodeSize;
				}
				size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
				decode(tmp, current_block_elements, root, type);

				cur_unpred_count = decompressDataSeries_float_3D_RA_block_3D_pred_multi_means(data_pos, mean_count, means, dense_pos, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, type, unpredictable_data);

				if(cur_unpred_count != unpredictable_count){
					printf("Check bugs, unpredictable_count is not the same: %d %d\n", unpredictable_count, cur_unpred_count);
					printf("Current index: %d %d %d\n\n", i, j, k);
					for(int i=0; i<current_block_elements; i++){
						printf("%d ", type[i]);
					}
					printf("\n");
					exit(0);
				}

				comp_data_pos += *block_pos;
				block_pos ++;
				unpred_count ++;
				// printf("block comp done, data_offset from %ld to %ld: diff %ld\n", *data, data_pos, data_pos - *data);
				// fflush(stdout);

			}
		}
	}
	free(type);
	free(unpredictable_data);

}

void decompressDataSeries_float_3D_all_by_interpolation(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	size_t num_elements = r1 * r2 * r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	size_t block_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	// calculate block dims
	size_t num_x, num_y, num_z;
	num_x = r1 / block_size;
	num_y = r2 / block_size;
	num_z = r3 / block_size;

	if(r1 != num_x * block_size + 1){
		printf("dim x not aligned\n");
		exit(0);
	}
	if(r2 != num_y * block_size + 1){
		printf("dim y not aligned\n");
		exit(0);
	}
	if(r3 != num_z * block_size + 1){
		printf("dim z not aligned\n");
		exit(0);
	}
	size_t sample_strip0 = (num_y + 1) * (num_z + 1);
	size_t sample_strip1 = num_z + 1;
	size_t downsample_size = (num_x+1) * (num_y+1) * (num_z+1);

	double realPrecision = bytesToDouble(comp_data_pos);
	comp_data_pos += 8;
	unsigned int intervals = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;

	updateQuantizationInfo(intervals);
	unsigned int tree_size = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	allNodes = bytesToInt_bigEndian(comp_data_pos);
	stateNum = allNodes/2;
	SZ_Reset(allNodes, stateNum);
	// printf("Reconstruct huffman tree with node count %ld\n", nodeCount);
	// fflush(stdout);
	node root = reconstruct_HuffTree_from_bytes_anyStates(comp_data_pos+4, allNodes);
	comp_data_pos += 4 + tree_size;

	size_t unpredictable_count = *((size_t *) comp_data_pos);
	comp_data_pos += sizeof(size_t);
	float * unpredictable_data = (float *) comp_data_pos;
	comp_data_pos += unpredictable_count * sizeof(float);
	float * downsamples = (float *) comp_data_pos;
	comp_data_pos += downsample_size * sizeof(float);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	decode(comp_data_pos, num_elements, root, result_type);
	printf("SZ_TYPE_SIZE: %d\nsize: %ld\n\n", SZ_SIZE_TYPE, comp_data_pos - comp_data + MetaDataByteLength + 4 + SZ_SIZE_TYPE);
	for(int i=0; i<10; i++){
		printf("%d ", result_type[1000+i]);
	}
	printf("\n");
	printf("decompress start, unpredictable_count %zu, intervals %d, block_size %d\n", unpredictable_count, intervals, block_size);
	// {
	// 	int status;
	// 	writeFloatData_inBytes(downsamples, downsample_size, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/decomp_sample.dat", &status);
	// 	writeIntData_inBytes(result_type, num_elements, "/Users/LiangXin/Documents/research/anl/lossy_comp/data/NYX/decomp_type.dat", &status);
	// }
	// exit(0);

	int * type = result_type;
	// printf("decompress offset to start: %ld\n", comp_data_pos - tdps->data);
	// fflush(stdout);
	size_t index = 0;
	float * block_data_pos;
	float * cur_data_pos;
	int type_;
	float c00, c01, c10, c11, c0, c1;
	float pred;
	float * coeff = (float *) malloc((block_size+1) * sizeof(float));
	for(size_t i=0; i<=block_size; i++){
		coeff[i] = i * 1.0 / block_size;
	}
	int block_size_x, block_size_y, block_size_z;
	size_t unpred_count = 0;
	for(size_t i=0; i<num_x; i++){
		for(size_t j=0; j<num_y; j++){
			for(size_t k=0; k<num_z; k++){
				// printf("i, j, k: %d %d %d, index %zu\n", i, j, k, index);
				// fflush(stdout);

				block_data_pos = *data + i*block_size*dim0_offset + j*block_size*dim1_offset + k*block_size;
				cur_data_pos = block_data_pos;
				// trilinear interpolate block
				if(i == num_x - 1) block_size_x = block_size+1;
				else block_size_x = block_size;
				if(j == num_y - 1) block_size_y = block_size+1;
				else block_size_y = block_size;
				if(k == num_z - 1) block_size_z = block_size+1;
				else block_size_z = block_size;
				for(size_t ii=0; ii<block_size_x; ii++){
					for(size_t jj=0; jj<block_size_y; jj++){
						for(size_t kk=0; kk<block_size_z; kk++){
							// decompress
							type_ = type[index];
							// if(i==0 && j==0 && k==8){
							// 	printf("ii, jj, kk: %d %d %d, type_ %d\n", ii, jj, kk, type_);
							// 	fflush(stdout);
							// }
							if (type_ != 0){
								c00 = downsamples[i*sample_strip0 + j*sample_strip1 + k] * (1 - coeff[ii]) + downsamples[i*sample_strip0 + j*sample_strip1 + k+1] * coeff[ii];
								c01 = downsamples[i*sample_strip0 + (j+1)*sample_strip1 + k] * (1 - coeff[ii]) + downsamples[i*sample_strip0 + (j+1)*sample_strip1 + k+1] * coeff[ii];
								c10 = downsamples[(i+1)*sample_strip0 + j*sample_strip1 + k] * (1 - coeff[ii]) + downsamples[(i+1)*sample_strip0 + j*sample_strip1 + k+1] * coeff[ii];
								c11 = downsamples[(i+1)*sample_strip0 + (j+1)*sample_strip1 + k] * (1 - coeff[ii]) + downsamples[(i+1)*sample_strip0 + (j+1)*sample_strip1 + k+1] * coeff[ii];
								c0 = c00 * (1 - coeff[jj]) + c10 * coeff[jj];
								c1 = c01 * (1 - coeff[jj]) + c11 * coeff[jj];
								pred = c0 * (1 - coeff[ii]) + c1 * coeff[ii];
								*cur_data_pos = pred + 2 * (type_ - intvRadius) * realPrecision;
							}
							else{
								*cur_data_pos = unpredictable_data[unpred_count ++];
							}
							if(block_size*i+ii == 0 && block_size*j+jj == 1 && block_size*k+kk == 505){
								printf("type_ %d, index %d, data %f\n", type_, index, *cur_data_pos);
							}
							index ++;
							cur_data_pos ++;
						}
						cur_data_pos += dim1_offset - block_size_z;
					}
					cur_data_pos += dim0_offset - dim1_offset * block_size_y;
				}

			}
		}
	}
	printf("Actual unpredictable_count %d\n", unpred_count);
	free(result_type);
}
