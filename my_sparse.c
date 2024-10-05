#include "spmv.h"

int my_sparse(const unsigned int n, const double mat[], double vec[], double result[])
{
  for (int i = 0; i < n; i++){
    double sum = 0;
    for (int j = 0; j < n; j++){
      sum += mat[n*i + j] * vec[j];
    }
    result[i] = sum;
  }
  return 0;
}
