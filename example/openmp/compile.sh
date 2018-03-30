gcc -g -O0 -fopenmp -c sz_omp.c -o sz_omp -I/Users/LiangXin/utils/sz_dev/include
gcc -g -O0 -fopenmp sz_omp_test.c sz_omp -I/Users/LiangXin/utils/sz_dev/include /Users/LiangXin/utils/sz_dev/lib/libSZ.a /Users/LiangXin/utils/sz_dev/lib/libzlib.a
