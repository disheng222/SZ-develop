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
				;
			else if(dim == 3)
				decompressDataSeries_float_3D_nonblocked(newData, r3, r2, r1, tdps->raBytes);
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
			decompressDataSeries_float_1D_pwrgroup(data, dataSeriesLength, tdps);
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
			decompressDataSeries_float_2D_pwr(data, r1, r2, tdps);
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
			decompressDataSeries_float_3D_pwr(data, r1, r2, r3, tdps);
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
			decompressDataSeries_float_3D_pwr(data, r1*r2, r3, r4, tdps);
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
			decompressDataSeries_float_3D_pwr(&decmpData, r1*r2, r3, r4, tdps);
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
	COMPUTE_2D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	COMPUTE_2D_NUMBER_OF_BLOCKS(r2, num_y, block_size);

	size_t split_index_x, split_index_y;
	size_t early_blockcount_x, early_blockcount_y;
	size_t late_blockcount_x, late_blockcount_y;
	COLL_BASE_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);

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

	COMPUTE_1D_NUMBER_OF_BLOCKS(r1, num_x, block_size);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);

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

void decompressDataSeries_float_3D_RA(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
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

				// int cur_unpred_count = decompressDataSeries_float_3D_RA_block(data_pos, mean, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, type, unpredictable_data);
				cur_unpred_count = decompressDataSeries_float_3D_RA_block_3D_pred(data_pos, mean, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, type, unpredictable_data);
				//int cur_unpred_count = decompressDataSeries_float_3D_RA_block_1D_pred(data_pos, mean, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, type, unpredictable_data);

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
				mean_pos ++;

				// printf("block comp done, data_offset from %ld to %ld: diff %ld\n", *data, data_pos, data_pos - *data);
				// fflush(stdout);

			}
		}
	}
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
