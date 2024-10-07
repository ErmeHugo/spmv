#include "spmv.h"

// Function to multiply CSR matrix by vector
int my_sparse(const double *values, const int *col_i, const int *offset, 
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
