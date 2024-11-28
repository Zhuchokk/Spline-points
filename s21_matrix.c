#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int return_value = 0;
  if (rows < 1 || columns < 1) return_value = 1;
  if (return_value == 0) {
    result->rows = rows;
    result->columns = columns;
    result->matrix = malloc(result->rows * result->columns * sizeof(double) +
                            result->rows * sizeof(double *));
    double *ptr = (double *)(result->matrix + result->rows);
    for (int i = 0; i < result->rows; i++)
      result->matrix[i] = ptr + result->columns * i;
  }
  return return_value;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL && A->matrix != NULL) {
    free(A->matrix);
    A->matrix = NULL;
    A->rows = 0;
    A->columns = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int return_value = SUCCESS;
  if (A->rows == B->rows && A->columns == B->columns) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++)
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= EPS) {
          return_value = FAILURE;
          break;
        }
      if (return_value == FAILURE) break;
    }
  } else {
    return_value = FAILURE;
  }
  return return_value;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int return_value = 0;
  if (A == NULL || B == NULL || A->rows < 1 || A->columns < 1 || B->rows < 1 ||
      B->columns < 1)
    return_value = 1;
  else if (A->rows != B->rows || A->columns != B->columns)
    return_value = 2;
  if (return_value == 0) {
    if (s21_create_matrix(A->rows, A->columns, result) == 0) {
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++)
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    } else {
      return_value = 1;
    }
  }
  return return_value;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int return_value = 0;
  if (A == NULL || B == NULL || A->rows < 1 || A->columns < 1 || B->rows < 1 ||
      B->columns < 1)
    return_value = 1;
  else if (A->rows != B->rows || A->columns != B->columns)
    return_value = 2;
  if (return_value == 0) {
    if (s21_create_matrix(A->rows, A->columns, result) == 0) {
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++)
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    } else {
      return_value = 1;
    }
  }
  return return_value;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int return_value = 0;
  if (A == NULL || A->rows < 1 || A->columns < 1)
    return_value = 1;
  else
    return_value = s21_create_matrix(A->rows, A->columns, result);
  if (return_value == 0) {
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[i][j] = A->matrix[i][j] * number;
  }
  return return_value;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int return_value = 0;
  if (A == NULL || B == NULL || A->rows < 1 || A->columns < 1 || B->rows < 1 ||
      B->columns < 1)
    return_value = 1;
  else if (A->columns != B->rows)
    return_value = 2;
  if (return_value == 0) {
    if (s21_create_matrix(A->rows, B->columns, result) == 0) {
      double sum = 0;
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          sum = 0;
          for (int k = 0; k < A->columns; k++)
            sum += A->matrix[i][k] * B->matrix[k][j];
          result->matrix[i][j] = sum;
        }
      }
    } else {
      return_value = 1;
    }
  }
  return return_value;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int return_value = 0;
  if (A == NULL || A->rows < 1 || A->columns < 1)
    return_value = 1;
  else
    return_value = s21_create_matrix(A->columns, A->rows, result);
  if (return_value == 0) {
    for (int i = 0; i < A->columns; i++)
      for (int j = 0; j < A->rows; j++) result->matrix[i][j] = A->matrix[j][i];
  }
  return return_value;
}

int s21_determinant(matrix_t *A, double *result) {
  *result = 0;
  int return_value = 0;
  if (A == NULL || A->rows < 1 || A->columns < 1)
    return_value = 1;
  else if (A->rows != A->columns)
    return_value = 2;
  if (return_value == 0) {
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else if (A->rows == 2) {
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    } else {
      int temp_col = 0;
      double temp_det_num = 0;
      matrix_t temp_matrix;
      int temp_i = 0;
      int temp_j = 0;
      s21_create_matrix(A->rows - 1, A->rows - 1, &temp_matrix);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 1; j < A->rows; j++) {
          for (int k = 0; k < A->rows; k++)
            if (k != temp_col) {
              temp_matrix.matrix[temp_i][temp_j] = A->matrix[j][k];
              ++temp_j;
            }
          ++temp_i;
          temp_j = 0;
        }
        temp_i = 0;
        s21_determinant(&temp_matrix, &temp_det_num);
        *result += A->matrix[0][temp_col] * temp_det_num * pow(-1, temp_col);
        ++temp_col;
      }
      s21_remove_matrix(&temp_matrix);
    }
  }
  return return_value;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int return_value = 0;
  if (A == NULL || A->rows < 1 || A->columns < 1)
    return_value = 1;
  else if (A->rows != A->columns)
    return_value = 2;
  if (return_value == 0 && A->rows == 1) {
    s21_create_matrix(A->rows, A->rows, result);
    result->matrix[0][0] = A->matrix[0][0];
  } else if (return_value == 0) {
    s21_create_matrix(A->rows, A->rows, result);
    int res_i = 0;
    int res_j = 0;
    double temp_det_num = 0;
    matrix_t temp_matrix;
    int temp_i = 0;
    int temp_j = 0;
    s21_create_matrix(A->rows - 1, A->rows - 1, &temp_matrix);
    for (int k = 0; k < A->rows * A->rows; k++) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->rows; j++) {
          if (res_i != i && res_j != j) {
            temp_matrix.matrix[temp_i][temp_j] = A->matrix[i][j];
            ++temp_j;
          }
        }
        if (temp_j == A->rows - 1) {
          ++temp_i;
          temp_j = 0;
        }
      }
      temp_i = 0;
      s21_determinant(&temp_matrix, &temp_det_num);
      result->matrix[res_i][res_j] = temp_det_num * pow(-1, res_i + res_j);
      ++res_j;
      if (res_j == A->rows) {
        ++res_i;
        res_j = 0;
      }
    }
    s21_remove_matrix(&temp_matrix);
  }
  return return_value;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  double det = 0;
  int return_value = 0;
  if (A == NULL || A->rows < 1 || A->columns < 1) {
    return_value = 1;
  } else {
    s21_determinant(A, &det);
    if (det == 0 || A->rows != A->columns) return_value = 2;
  }
  if (return_value == 0) {
    double det_koef = 1 / det;
    matrix_t complements;
    s21_calc_complements(A, &complements);
    matrix_t transp;
    s21_transpose(&complements, &transp);
    s21_mult_number(&transp, det_koef, result);
    s21_remove_matrix(&complements);
    s21_remove_matrix(&transp);
  }
  return return_value;
}
