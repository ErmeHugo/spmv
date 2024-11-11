#include "spmv.h"

// Function to multiply CSR matrix by vector
int my_sparse_csr(const double *values, const int *col_i, const int *offset, 
              const double *vec, double *result, int size)
  {
    // Initialize the result vector to zero
      for (int i = 0; i < size; ++i) {
          result[i] = 0.0;
      }

      // Perform the matrix-vector multiplication using CSR format
      for (int i = 0; i < size; ++i) {
          for (int j = offset[i]; j < offset[i + 1]; ++j) {
              result[i] += values[j] * vec[col_i[j]];
          }
      }
      return 1;
  }

// Function to multiply CSC matrix by vector
int my_sparse_csc(const double *values, const int *row_j, const int *offset, 
              const double *vec, double *result, int size)
  {
    // Initialize the result vector to zero
      for (int i = 0; i < size; ++i) {
          result[i] = 0.0;
      }

      // Perform the matrix-vector multiplication using CSC format
        for (int col = 0; col < size; ++col) {
            for (int j = offset[col]; j < offset[col + 1]; ++j) {
                result[row_j[j]] += values[j] * vec[col];
            }
        }
      
      return 1;
  }

// Function to multiply COO matrix by vector
int my_sparse_coo(const double *values, const int *row_i, const int *col_j, 
              const double *vec, double *result, int nnz, int size)
{
    // Initialize the result vector to zero
    for (int i = 0; i < size; ++i) {
        result[i] = 0.0;
    }

    // Perform the matrix-vector multiplication using COO format
    for (int i = 0; i < nnz; ++i) {
        result[row_i[i]] += values[i] * vec[col_j[i]];
    }
    
    return 1;
}