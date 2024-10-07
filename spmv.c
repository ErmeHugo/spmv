#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>      // CBLAS in GSL (the GNU Scientific Library)
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include "timer.h"
#include "spmv.h"

#define DEFAULT_SIZE 1024
#define DEFAULT_DENSITY 0.25

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

  // Convert mat to a sparse format: CSR
  CSR_Matrix csr;
  csr = convert_to_csr(mat,size,nnz);
  /*
  // Use the gsl_spmatrix struct as datatype
  gsl_matrix_view gsl_mat_view = gsl_matrix_view_array(mat, size, size);
  gsl_vector_view gsl_vec_view = gsl_vector_view_array(vec, size);

  gsl_matrix *gsl_mat = &gsl_mat_view.matrix;
  gsl_vector *gsl_vec = &gsl_vec_view.vector;
  gsl_spmatrix *gsl_mat_csr;

  gsl_spmatrix_csr(*gsl_mat_csr, *gsl_mat);
  //
  // Sparse computation using GSL's sparse algebra functions
  //
  printf("Dense computation\n----------------\n");

  timestamp(&start);
  gsl_vector *result = gsl_vector_alloc(size);

  // Perform the matrix-vector multiplication 
  cgsl_cblas_dgemv(CblasNoTrans, 1.0, gsl_mat, gsl_vec, 0.0, result);

  timestamp(&now);
  printf("Time taken by GSL CBLAS dense computation: %ld ms\n", diff_milli(&start, &now));
  */
  //
  // Your own sparse implementation
  //
  my_sparse(csr.values, csr.col_i, csr.offset, vec, mysol, size);
  // Compare times (and computation correctness!)


  // Free resources
  free(mat);
  free(vec);
  free(refsol);
  free(mysol);

  return 0;
}
