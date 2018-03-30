#include "sz_omp.h"
#include <math.h>
#include <time.h>

unsigned char * SZ_compress_float_3D_MDQ_openmp(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t * comp_size){

	double elapsed_time = 0.0;

	elapsed_time = -omp_get_wtime();
	unsigned int quantization_intervals;
	float dense_pos;
	float mean_flush_freq;
	float sz_sample_correct_freq;
	if(optQuantMode==1)
	{
		// quantization_intervals = optimize_intervals_float_3D(oriData, r1, realPrecision);
		quantization_intervals = optimize_intervals_float_3D_with_freq_and_dense_pos(oriData, r1, r2, r3, realPrecision, &dense_pos, &sz_sample_correct_freq, &mean_flush_freq);
		printf("3D number of bins: %d\nerror bound %.20f dense position %.20f\n", quantization_intervals, realPrecision, dense_pos);
		// exit(0);		
		updateQuantizationInfo(quantization_intervals);
	}	
	else{
		quantization_intervals = intvCapacity;
	}
	elapsed_time += omp_get_wtime();
	printf("opt interval time: %.4f\n", elapsed_time);

	elapsed_time = -omp_get_wtime();
	int thread_num = omp_get_max_threads();
	thread_num = 4;
	int thread_order = (int)log2(thread_num);
	size_t num_x, num_y, num_z;
	{
		int block_thread_order = thread_order / 3;
		switch(thread_order % 3){
			case 0:{
				num_x = 1 << block_thread_order;
				num_y = 1 << block_thread_order;
				num_z = 1 << block_thread_order;
				break;
			}
			case 1:{
				num_x = 1 << (block_thread_order + 1);
				num_y = 1 << block_thread_order;
				num_z = 1 << block_thread_order;
				break;
			}
			case 2:{
				num_x = 1 << (block_thread_order + 1);
				num_y = 1 << (block_thread_order + 1);
				num_z = 1 << block_thread_order;
				break;
			}
		}
		thread_num = num_x * num_y * num_z;
	}
	omp_set_num_threads(thread_num);
	// calculate block dims
	printf("number of blocks: %d %d %d\n", num_x, num_y, num_z);

	size_t split_index_x, split_index_y, split_index_z;
	size_t early_blockcount_x, early_blockcount_y, early_blockcount_z;
	size_t late_blockcount_x, late_blockcount_y, late_blockcount_z;
	COLL_BASE_COMPUTE_BLOCKCOUNT(r1, num_x, split_index_x, early_blockcount_x, late_blockcount_x);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r2, num_y, split_index_y, early_blockcount_y, late_blockcount_y);
	COLL_BASE_COMPUTE_BLOCKCOUNT(r3, num_z, split_index_z, early_blockcount_z, late_blockcount_z);

	size_t max_num_block_elements = early_blockcount_x * early_blockcount_y * early_blockcount_z;
	size_t num_blocks = num_x * num_y * num_z;
	size_t num_elements = r1 * r2 * r3;
	// printf("max_num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	
	// printf("malloc blockinfo array start\n");
	// fflush(stdout);

	size_t buffer_size = early_blockcount_y * early_blockcount_z * sizeof(float);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	size_t unpred_data_max_size = max_num_block_elements;
	float * result_unpredictable_data = (float *) malloc(unpred_data_max_size * sizeof(float) * num_blocks);

	unsigned int * unpredictable_count = (unsigned int *) malloc(num_blocks * sizeof(unsigned int));
	float * mean = malloc(num_blocks * sizeof(float));
	int num_yz = num_y * num_z;
	#pragma omp parallel for
	for(int t=0; t<thread_num; t++){
		int id = omp_get_thread_num();
		int i = id/(num_yz);
		int j = (id % num_yz) / num_y;
		int k = id % num_z;
		// printf("%d: %d %d %d\n", omp_get_thread_num(), i, j, k);
		size_t offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
		size_t offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
		size_t offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
		float * data_pos = oriData + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

		size_t current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
		size_t current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
		size_t current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
		size_t type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
		int * type = result_type + type_offset;

		float * unpredictable_data = result_unpredictable_data + id * unpred_data_max_size;
		float *P0, *P1; // buffer
		P0 = (float *) malloc(buffer_size);
		P1 = (float *) malloc(buffer_size);
		unpredictable_count[id] = SZ_compress_float_3D_MDQ_RA_block(data_pos, mean + id, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, P0, P1, type, unpredictable_data);
		free(P0);
		free(P1);
	}
	elapsed_time += omp_get_wtime();
	printf("compression and quantization time: %.4f\n", elapsed_time);
	elapsed_time = -omp_get_wtime();
	// printf("unpred count:\n");
	// for(int i=0; i<num_blocks; i++){
	// 	printf("%d ", unpredictable_count[i]);
	// }
	// printf("\n");
	// printf("total_unpred num: %d\n", total_unpred);
	// printf("Block wise compression end, num_elements %ld\n", num_elements);
	// huffman encode
	SZ_Reset(allNodes, stateNum);
	size_t nodeCount = 0;
	init(result_type, num_elements);
	elapsed_time += omp_get_wtime();
	printf("Build Huffman: %.4f\n", elapsed_time);
	elapsed_time = -omp_get_wtime();
	for (size_t i = 0; i < stateNum; i++)
		if (code[i]) nodeCount++;
	nodeCount = nodeCount*2-1;
	unsigned char *treeBytes;
	unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(nodeCount, &treeBytes);

	unsigned int meta_data_offset = 3 + 1 + MetaDataByteLength;
	size_t total_unpred = 0;
	for(int i=0; i<num_blocks; i++){
		total_unpred += unpredictable_count[i];
		// printf("%d: %d mean %.2f\n", i, unpredictable_count[i], mean[i]);
	}
	// total size 										metadata		real precision		intervals	nodeCount		huffman 	 	block index 						unpredicatable count						mean 					 	unpred size 				elements
	unsigned char * result = (unsigned char *) malloc(meta_data_offset + sizeof(double) + sizeof(int) + sizeof(int) + treeByteSize + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(unsigned short) + num_blocks * sizeof(float) + total_unpred * sizeof(float) + num_elements * sizeof(int));
	unsigned char * result_pos = result;
	initRandomAccessBytes(result_pos);
	result_pos += meta_data_offset;

	size_t enCodeSize = 0;

	intToBytes_bigEndian(result_pos, thread_num);
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

	memcpy(result_pos, unpredictable_count, num_blocks * sizeof(unsigned int));
	result_pos += num_blocks * sizeof(unsigned int);
	memcpy(result_pos, mean, num_blocks * sizeof(float));
	result_pos += num_blocks * sizeof(float);	
	// printf("unpred offset: %ld\n", result_pos - result);
	// storage unpredicable data
	float * unpred_pos = (float *) result_pos;
	for(int t=0; t<thread_num; t++){
		float * unpredictable_data = result_unpredictable_data + t * unpred_data_max_size;
		memcpy(result_pos, unpredictable_data, unpredictable_count[t] * sizeof(float));		
		result_pos += unpredictable_count[t]*sizeof(float);
	}
	elapsed_time += omp_get_wtime();
	printf("write misc time: %.4f\n", elapsed_time);
	elapsed_time = -omp_get_wtime();

	size_t * block_pos = (size_t *) result_pos;
	result_pos += num_blocks * sizeof(size_t);
	unsigned char * encoding_buffer = (unsigned char *) malloc(max_num_block_elements * sizeof(int) * num_blocks);
	#pragma omp parallel for
	for(int t=0; t<thread_num; t++){
		int id = omp_get_thread_num();
		int i = id/(num_yz);
		int j = (id % num_yz) / num_z;
		int k = id % num_z;
		unsigned char * encoding_buffer_pos = encoding_buffer + id * max_num_block_elements * sizeof(int);
		size_t enCodeSize = 0;
		size_t offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
		size_t offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
		size_t offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
		float * data_pos = oriData + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;
		size_t current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
		size_t current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
		size_t current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
		size_t current_block_elements = current_blockcount_x * current_blockcount_y * current_blockcount_z;
		size_t type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
		int * type = result_type + type_offset;
		encode(type, current_block_elements, encoding_buffer_pos, &enCodeSize);
		block_pos[id] = enCodeSize;
	}
	elapsed_time += omp_get_wtime();
	printf("Parallel Huffman encoding elapsed time: %.4f\n", elapsed_time);
	elapsed_time = -omp_get_wtime();
	for(int t=0; t<thread_num; t++){
		memcpy(result_pos, encoding_buffer + t * max_num_block_elements * sizeof(int), block_pos[t]);
		result_pos += block_pos[t];
	}
	elapsed_time += omp_get_wtime();
	printf("Final copy elapsed time: %.4f\n", elapsed_time);
	free(encoding_buffer);
	// {
	// 	int status;
	// 	writeIntData_inBytes(result_type, num_elements, "/Users/LiangXin/github/SZ-develop/example/openmp/comp001_type.dat", &status);
	// }

	// int status;
	// writeIntData_inBytes(result_type, num_elements, "/Users/LiangXin/github/SZ-develop/example/openmp/omp_type.dat", &status);
	// printf("type array size: %ld\n", enCodeSize);
	result_pos += enCodeSize;
	size_t totalEncodeSize = 0;
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

void decompressDataSeries_float_3D_openmp(float** data, size_t r1, size_t r2, size_t r3, unsigned char* comp_data){
	// printf("num_block_elements %d num_blocks %d\n", max_num_block_elements, num_blocks);
	// fflush(stdout);
	double elapsed_time = 0.0;
	elapsed_time = -omp_get_wtime();

	size_t dim0_offset = r2 * r3;
	size_t dim1_offset = r3;
	size_t num_elements = r1 * r2 * r3;

	*data = (float*)malloc(sizeof(float)*num_elements);
	
	unsigned char * comp_data_pos = comp_data;
	//int meta_data_offset = 3 + 1 + MetaDataByteLength;
	//comp_data_pos += meta_data_offset;

	int thread_num = bytesToInt_bigEndian(comp_data_pos);
	comp_data_pos += 4;
	int thread_order = (int)log2(thread_num);
	size_t num_x, num_y, num_z;
	{
		int block_thread_order = thread_order / 3;
		switch(thread_order % 3){
			case 0:{
				num_x = 1 << block_thread_order;
				num_y = 1 << block_thread_order;
				num_z = 1 << block_thread_order;
				break;
			}
			case 1:{
				num_x = 1 << (block_thread_order + 1);
				num_y = 1 << block_thread_order;
				num_z = 1 << block_thread_order;
				break;
			}
			case 2:{
				num_x = 1 << (block_thread_order + 1);
				num_y = 1 << (block_thread_order + 1);
				num_z = 1 << block_thread_order;
				break;
			}
		}
	}
	printf("number of blocks: %d %d %d, thread_num %d\n", num_x, num_y, num_z, thread_num);
	omp_set_num_threads(thread_num);
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
	unsigned int * unpred_count = (unsigned int *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(unsigned int);
	float * mean_pos = (float *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(float);
	float * result_unpredictable_data = (float *) comp_data_pos;
	size_t total_unpred = 0;
	size_t * unpred_offset = (size_t *) malloc(num_blocks * sizeof(size_t));
	for(int i=0; i<num_blocks; i++){
		unpred_offset[i] = total_unpred;
		total_unpred += unpred_count[i];
	}
	comp_data_pos += total_unpred * sizeof(float);

	// printf("unpred count:\n");
	// for(int i=0; i<num_blocks; i++){
	// 	printf("%d ", unpred_count[i]);
	// }
	// printf("\n");
	// for(int i=0; i<1000; i++){
	// 	printf("%.2f ", result_unpredictable_data[i]);
	// }
	// printf("\ntotal_unpred num: %d\n", total_unpred);
	
	// for(int i=0; i<num_blocks; i++){
	// 	printf("%d unpred offset %ld\n", i, unpred_offset[i]);
	// 	for(int tmp=0; tmp<10; tmp++){
	// 		printf("%.2f ", (result_unpredictable_data + unpred_offset[i])[tmp]);
	// 	}
	// 	printf("\n");
	// }
	// exit(0);
	// printf("Block wise decompression start: %d %d %d\n", early_blockcount_x, early_blockcount_y, early_blockcount_z);
	// fflush(stdout);
	int * result_type = (int *) malloc(num_elements * sizeof(int));
	// decode(comp_data_pos, num_elements, root, result_type);
	size_t * block_offset = (size_t *) malloc(num_blocks * sizeof(size_t));
	size_t * block_pos = (size_t *) comp_data_pos;
	comp_data_pos += num_blocks * sizeof(size_t);
	block_offset[0] = 0;
	for(int t=1; t<thread_num; t++){
		block_offset[t] = block_pos[t-1] + block_offset[t-1];
	}
	int num_yz = num_y * num_z;
	elapsed_time += omp_get_wtime();
	printf("Read data info elapsed time: %.4f\n", elapsed_time);
	elapsed_time = -omp_get_wtime();
	#pragma omp parallel for
	for(int t=0; t<thread_num; t++){
		int id = omp_get_thread_num();
		int i = id/(num_yz);
		int j = (id % num_yz) / num_z;
		int k = id % num_z;
		size_t offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
		size_t offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
		size_t offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
		size_t current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
		size_t current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
		size_t current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
		int type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
		int * type = result_type + type_offset;
		decode(comp_data_pos + block_offset[t], current_blockcount_x*current_blockcount_y*current_blockcount_z, root, type);
	}
	elapsed_time += omp_get_wtime();
	printf("Parallel Huffman decoding elapsed time: %.4f\n", elapsed_time);
	elapsed_time = -omp_get_wtime();

	#pragma omp parallel for
	for(int t=0; t<thread_num; t++){
		int id = omp_get_thread_num();
		int i = id/(num_yz);
		int j = (id % num_yz) / num_y;
		int k = id % num_z;
		// printf("%d: %d %d %d\n", omp_get_thread_num(), i, j, k);
		size_t offset_x = (i < split_index_x) ? i * early_blockcount_x : i * late_blockcount_x + split_index_x;
		size_t offset_y = (j < split_index_y) ? j * early_blockcount_y : j * late_blockcount_y + split_index_y;
		size_t offset_z = (k < split_index_z) ? k * early_blockcount_z : k * late_blockcount_z + split_index_z;
		float * data_pos = *data + offset_x * dim0_offset + offset_y * dim1_offset + offset_z;

		size_t current_blockcount_x = (i < split_index_x) ? early_blockcount_x : late_blockcount_x;
		size_t current_blockcount_y = (j < split_index_y) ? early_blockcount_y : late_blockcount_y;
		size_t current_blockcount_z = (k < split_index_z) ? early_blockcount_z : late_blockcount_z;
		int type_offset = offset_x * dim0_offset +  offset_y * current_blockcount_x * dim1_offset + offset_z * current_blockcount_x * current_blockcount_y;
		int * type = result_type + type_offset;

		float * unpredictable_data = result_unpredictable_data + unpred_offset[id];
		float mean = mean_pos[id];
		// printf("\n%d\ndata_offset: %ld\n", t, offset_x * dim0_offset + offset_y * dim1_offset + offset_z);
		// printf("mean: %.2f\n", mean);
		// for(int tmp=0; tmp<10; tmp++){
		// 	printf("%.2f ", unpredictable_data[tmp]);
		// }
		// printf("\n\n");
		int cur_unpred_count = decompressDataSeries_float_3D_RA_block(data_pos, mean, r1, r2, r3, current_blockcount_x, current_blockcount_y, current_blockcount_z, realPrecision, type, unpredictable_data);
	}	
	elapsed_time += omp_get_wtime();
	printf("Parallel decompress elapsed time: %.4f\n", elapsed_time);

	free(result_type);
	free(unpred_offset);

}





