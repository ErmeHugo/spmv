# SpMV: Sparse Matrix-Vector product

This code is based on the use of GSL (GNU Scientific Library) and MKL 
(MKL Math Kernel Library) for the implementation of the baseline 
operations used for comparison:
- dense matrix-vector product: `cblas_dgemv()`
- sparse matrix-vector product: `gsl_spblas_dgemv()`, `mkl_sparse_d_mv()`


# Compilation and optimization

For this project we have a Makefile, to compile you can just write :

```bash 
 make 
 ```

If you want to add optimization for the program you can do so like this :

```bash 
#GCC Compiler
make OPT="-O2 -fno-tree-vectorize"   
make OPT="-O3 -ftree-vectorize"      
make OPT="-Ofast -ftree-vectorize" 

#ICC Compiler
make OPT=-Ofast     
make OPT="-O2 -no-vec" 
make OPT="-O3 -xAVX"   
 ```

In order to use the mkl library please put the line 13 in comments and 
decomment line 14 of the Makefile

The file spmv_mkl.c is used for the calculation of mkl reference times


# Benchmarking Results 

Benchmarking results are located in the results.ods file and an example of 
the results can be found in the log files.

