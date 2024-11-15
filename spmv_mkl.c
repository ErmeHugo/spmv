#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>       
#include "timer.h"
#include "spmv.h"

#define DEFAULT_SIZE 16384
#define DEFAULT_DENSITY 0.1

unsigned int populate_sparse_matrix(double mat[], unsigned int n, double density, unsigned int seed)
{
    unsigned int nnz = 0;

    srand(seed);

    for (unsigned int i = 0; i < n * n; i++) {
        if ((rand() % 100) / 100.0 < density) {
            // Get a pseudorandom value between -9.99 e 9.99
            mat[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
            nnz++;
        } else {
            mat[i] = 0;
        }
    }

    return nnz;
}

unsigned int populate_vector(double vec[], unsigned int size, unsigned int seed)
{
    srand(seed);

    for (unsigned int i = 0; i < size; i++) {
        vec[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
    }

    return size;
}

void spmv_mkl_csr(double *mat, double *vec, double *result, unsigned int n, unsigned int nnz)
{
    sparse_matrix_t mkl_matrix;
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;

    MKL_INT *rows_start = (MKL_INT *)malloc((n + 1) * sizeof(MKL_INT));
    MKL_INT *columns = (MKL_INT *)malloc(nnz * sizeof(MKL_INT));
    double *values = (double *)malloc(nnz * sizeof(double));

    unsigned int k = 0;
    rows_start[0] = 0;
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            if (mat[i * n + j] != 0) {
                columns[k] = j;
                values[k] = mat[i * n + j];
                k++;
            }
        }
        rows_start[i + 1] = k;
    }

    mkl_sparse_d_create_csr(&mkl_matrix, SPARSE_INDEX_BASE_ZERO, n, n, rows_start, rows_start + 1, columns, values);

    // Perform sparse matrix-vector multiplication
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, mkl_matrix, descr, vec, 0.0, result);

    mkl_sparse_destroy(mkl_matrix);
    free(rows_start);
    free(columns);
    free(values);
}

void spmv_mkl_csc(double *mat, double *vec, double *result, unsigned int n, unsigned int nnz)
{
    sparse_matrix_t mkl_matrix;
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;

    MKL_INT *cols_start = (MKL_INT *)malloc((n + 1) * sizeof(MKL_INT));
    MKL_INT *rows = (MKL_INT *)malloc(nnz * sizeof(MKL_INT));
    double *values = (double *)malloc(nnz * sizeof(double));

    unsigned int k = 0;
    cols_start[0] = 0;
    for (unsigned int j = 0; j < n; j++) {
        for (unsigned int i = 0; i < n; i++) {
            if (mat[i * n + j] != 0) {
                rows[k] = i;
                values[k] = mat[i * n + j];
                k++;
            }
        }
        cols_start[j + 1] = k;
    }

    mkl_sparse_d_create_csc(&mkl_matrix, SPARSE_INDEX_BASE_ZERO, n, n, cols_start, cols_start + 1, rows, values);

    // Perform sparse matrix-vector multiplication
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, mkl_matrix, descr, vec, 0.0, result);

    mkl_sparse_destroy(mkl_matrix);
    free(cols_start);
    free(rows);
    free(values);
}

void spmv_mkl_coo(double *mat, double *vec, double *result, unsigned int n, unsigned int nnz)
{
    sparse_matrix_t mkl_matrix;
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;

    MKL_INT *row_indices = (MKL_INT *)malloc(nnz * sizeof(MKL_INT));
    MKL_INT *col_indices = (MKL_INT *)malloc(nnz * sizeof(MKL_INT));
    double *values = (double *)malloc(nnz * sizeof(double));

    unsigned int k = 0;
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            if (mat[i * n + j] != 0) {
                row_indices[k] = i;
                col_indices[k] = j;
                values[k] = mat[i * n + j];
                k++;
            }
        }
    }

    mkl_sparse_d_create_coo(&mkl_matrix, SPARSE_INDEX_BASE_ZERO, n, n, nnz, row_indices, col_indices, values);

    // Perform sparse matrix-vector multiplication
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, mkl_matrix, descr, vec, 0.0, result);

    mkl_sparse_destroy(mkl_matrix);
    free(row_indices);
    free(col_indices);
    free(values);
}

int main(int argc, char *argv[])
{

    long diff_milli(struct timespec *start, struct timespec *end) {
      return (end->tv_sec - start->tv_sec) * 1000 + (end->tv_nsec - start->tv_nsec) / 1000000;
    }

    unsigned int size = DEFAULT_SIZE;
    double density = DEFAULT_DENSITY;
    unsigned int seed = 42;

    double *mat = (double *)malloc(size * size * sizeof(double));
    double *vec = (double *)malloc(size * sizeof(double));
    double *result = (double *)malloc(size * sizeof(double));

    unsigned int nnz = populate_sparse_matrix(mat, size, density, seed);
    populate_vector(vec, size, seed);

    printf("Running CSR format:\n");
    timeinfo start, now;
    timestamp(&start);
    spmv_mkl_csr(mat, vec, result, size, nnz);
    timestamp(&now);
    printf("Time taken by MKL CSR computation: %ld ms\n", diff_milli(&start, &now));

    printf("\nRunning CSC format:\n");
    timestamp(&start);
    spmv_mkl_csc(mat, vec, result, size, nnz);
    timestamp(&now);
    printf("Time taken by MKL CSC computation: %ld ms\n", diff_milli(&start, &now));

    printf("\nRunning COO format:\n");
    timestamp(&start);
    spmv_mkl_coo(mat, vec, result, size, nnz);
    timestamp(&now);
    printf("Time taken by MKL COO computation: %ld ms\n", diff_milli(&start, &now));

    free(mat);
    free(vec);
    free(result);

    return 0;
}
