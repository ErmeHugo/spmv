# SpMV: Sparse Matrix-Vector product

This code is based on the use of GSL (GNU Scientific Library) for the
implementation of the baseline operations used for comparison:
- dense matrix-vector product: `cblas_dgemv()`
- sparse matrix-vector product: `gsl_spblas_dgemv()`

# Compilation 

For this project we have a Makefile, to compile you can just write :

```bash 
 make 
 ```
