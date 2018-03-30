#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include <math.h>

int main(){

	printf("Max thread number=%d %d\n", omp_get_max_threads(), (int)log2(omp_get_max_threads()));
	omp_set_num_threads(4);
	#pragma omp parallel
	{
		printf("Thread number=%d\n", omp_get_thread_num());
	}
}