#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>      // CBLAS in GSL (the GNU Scientific Library)
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
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

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return fabs(x - y) <= epsilon * fabs(x);
  // see Knuth section 4.2.2 pages 217-218
}

unsigned int check_result(double ref[], double result[], unsigned int size)
{
  for(unsigned int i = 0; i < size; i++) {
    if (!is_nearly_equal(ref[i], result[i]))
      return 0;
  }

  return 1;
}

typedef struct {
    double *values;   // Non-zero values
    int *col_i;       // Column indices of non-zero values
    int *offset;      // Row pointer array
} CSR_Matrix;

typedef struct {
    double *values;   // Non-zero values
    int *col_j;       // Column indices of non-zero values
    int *row_i;      // Row indices of non-zero values
} COO_Matrix;

typedef struct {
    double *values;   // Non-zero values
    int *offset;       // Column pointer array
    int *row_j;      // Row indices of non-zero values
} CSC_Matrix;

CSR_Matrix convert_to_csr(double mat[], int size,int nnz)
{
  CSR_Matrix csr;

  csr.offset = (int *)malloc((size+1)*sizeof(int));
  csr.col_i = (int *)malloc(nnz*sizeof(int));
  csr.values = (double *)malloc(nnz*sizeof(double));
  
  int count = 0;
  for (int i=0; i<size; i++){
    csr.offset[i] = count;
    for (int j=0; j<size; j++){
      if (mat[size*i + j] != 0){
        csr.values[count] = mat[size*i + j];
        csr.col_i[count] = j;
        count++;
      }
    }
  }
  csr.offset[size] = count;
  return csr;
}

COO_Matrix convert_to_coo(double mat[], int size, int nnz) {
    COO_Matrix coo;

    coo.row_i = (int *)malloc(nnz * sizeof(int));
    coo.col_j = (int *)malloc(nnz * sizeof(int));
    coo.values = (double *)malloc(nnz * sizeof(double));

    // Fill COO structure
    int count = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (mat[i * size + j] != 0) {
                coo.row_i[count] = i;
                coo.col_j[count] = j;
                coo.values[count] = mat[i * size + j];
                count++;
            }
        }
    }
    return coo;
}

CSC_Matrix convert_to_csc(double mat[], int size, int nnz) {
    CSC_Matrix csc;


    csc.offset = (int *)malloc((size + 1) * sizeof(int));
    csc.row_j = (int *)malloc(nnz * sizeof(int));
    csc.values = (double *)malloc(nnz * sizeof(double));

    // Initialize column offsets
    for (int j = 0; j <= size; j++) {
        csc.offset[j] = 0;
    }

    // Count non-zero entries in each column for col_offset
    for (int j = 0; j < size; j++) {
        for (int i = 0; i < size; i++) {
            if (mat[i * size + j] != 0) {
                csc.offset[j + 1]++;
            }
        }
    }

    // Compute cumulative sum for col_offset
    for (int j = 1; j <= size; j++) {
        csc.offset[j] += csc.offset[j - 1];
    }

    // Fill CSC structure
    int *col_counts = (int *)calloc(size, sizeof(int));
    for (int j = 0; j < size; j++) {
        for (int i = 0; i < size; i++) {
            if (mat[i * size + j] != 0) {
                int index = csc.offset[j] + col_counts[j];
                csc.row_j[index] = i;
                csc.values[index] = mat[i * size + j];
                col_counts[j]++;
            }
        }
    }

    return csc;
}

int main(int argc, char *argv[])
{
  int size;        // number of rows and cols (size x size matrix)
  double density;  // aprox. ratio of non-zero values

  if (argc < 2) {
    size = DEFAULT_SIZE;
    density = DEFAULT_DENSITY;
  } else if (argc < 3) {
    size = atoi(argv[1]);
    density = DEFAULT_DENSITY;
  } else {
    size = atoi(argv[1]);
    density = (double) atoi(argv[2]) / 100.0;
  }

  double *mat, *vec, *refsol, *mysol;

  mat = (double *) malloc(size * size * sizeof(double));
  vec = (double *) malloc(size * sizeof(double));
  refsol = (double *) malloc(size * sizeof(double));
  mysol = (double *) malloc(size * sizeof(double));

  unsigned int nnz = populate_sparse_matrix(mat, size, density, 1);
  populate_vector(vec, size, 2);

  printf("Matriz size: %d x %d (%d elements)\n", size, size, size*size);
  printf("%d non-zero elements (%.2lf%%)\n\n", nnz, (double) nnz / (size*size) * 100.0);

  //
  // Dense computation using CBLAS (eg. GSL's CBLAS implementation)
  //
  printf("Dense computation\n----------------\n");

  timeinfo start, now;
  timestamp(&start);

  cblas_dgemv(CblasRowMajor, CblasNoTrans, size, size, 1.0, mat, size, vec, 1, 0.0, refsol, 1);

  timestamp(&now);
  printf("Time taken by CBLAS dense computation: %ld ms\n", diff_milli(&start, &now));

  //
  // Using your own dense implementation
  //
  timestamp(&start);

  my_dense(size, mat, vec, mysol);

  timestamp(&now);
  printf("Time taken by my dense matrix-vector product: %ld ms\n", diff_milli(&start, &now));

  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");


  //
  // Let's try now SpMV: Sparse Matrix - Dense Vector computation
  //

  // Use the gsl_spmatrix struct as datatype
  gsl_spmatrix *sp = gsl_spmatrix_alloc(size, size);
  gsl_spmatrix *sp_csr;
  gsl_spmatrix *sp_csc;
  

   
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (mat[i * size + j] != 0) {
        gsl_spmatrix_set(sp, i, j, mat[i * size + j]);
        }
    }
  }

  // Copy the matrix into the ccs and crs format
  sp_csc = gsl_spmatrix_ccs(sp);
  sp_csr = gsl_spmatrix_crs(sp);
  
  gsl_vector *x = gsl_vector_alloc(size);
  for (int i = 0; i < size; i++) {
      gsl_vector_set(x, i, vec[i]);
  }

  // Allocate a result vector to store the multiplication result
  gsl_vector *result = gsl_vector_alloc(size);

  //
  // Sparse computation using GSL's sparse algebra functions
  //

  
  timestamp(&start);
  gsl_spblas_dgemv(CblasNoTrans, 1, sp, x, 0, result);
  timestamp(&now);
  printf("Time taken by SPBlas_coo sparse computation: %ld ms\n", diff_milli(&start, &now));
  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");

  timestamp(&start);
  gsl_spblas_dgemv(CblasNoTrans, 1, sp_csc, x, 0, result);
  timestamp(&now);
  printf("Time taken by SPBlas_csc sparse computation: %ld ms\n", diff_milli(&start, &now));
  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");

  timestamp(&start);
  gsl_spblas_dgemv(CblasNoTrans, 1, sp_csr, x, 0, result);
  timestamp(&now);
  printf("Time taken by SPBlas_csr sparse computation: %ld ms\n", diff_milli(&start, &now));
  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");

  
  // Convert mat to a sparse format: CSR
  CSR_Matrix csr;
  csr = convert_to_csr(mat,size,nnz);

  // Your own sparse implementation
  timestamp(&start);
  my_sparse_csr(csr.values, csr.col_i, csr.offset, vec, mysol, size);
  timestamp(&now);
  printf("Time taken by my sparse_csr matrix-vector product: %ld ms\n", diff_milli(&start, &now));
  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");


  // Convert mat to a sparse format: CSC
  CSC_Matrix csc;
  csc = convert_to_csc(mat,size,nnz);

  // Your own sparse implementation
  timestamp(&start);
  my_sparse_csc(csc.values, csc.row_j, csc.offset, vec, mysol, size);
  timestamp(&now);
  printf("Time taken by my sparse_csc matrix-vector product: %ld ms\n", diff_milli(&start, &now));
  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");


  // Convert mat to a sparse format: COO
  COO_Matrix coo;
  coo = convert_to_coo(mat,size,nnz);

  // Your own sparse implementation
  timestamp(&start);
  my_sparse_coo(coo.values, coo.row_i, coo.col_j, vec, mysol, nnz, size);
  timestamp(&now);
  printf("Time taken by my sparse_coo matrix-vector product: %ld ms\n", diff_milli(&start, &now));
  if (check_result(refsol, mysol, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");


  // Free resources
  free(mat);
  free(vec);
  free(refsol);
  free(mysol);
  return 0;
}