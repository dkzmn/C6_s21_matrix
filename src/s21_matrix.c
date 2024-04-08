#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int res = 0;
  if (rows <= 0 || columns <= 0 || !result) {
    res = 1;
  } else {
    result->matrix =
        malloc(columns * rows * sizeof(double) + rows * sizeof(double *));
    if (result->matrix) {
      result->rows = rows;
      result->columns = columns;
      double *ptr = (double *)(result->matrix + rows);
      for (int i = 0; i < rows; i++) result->matrix[i] = ptr + columns * i;
    } else {
      res = 1;
    }
  }
  return res;
}

void s21_remove_matrix(matrix_t *A) {
  if (A) {
    if (A->matrix) {
      free(A->matrix);
      A->columns = 0;
      A->rows = 0;
    }
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;
  if (s21_check_matrix(A) || s21_check_matrix(B)) {
    res = FAILURE;
  } else {
    if (A->rows != B->rows || A->columns != B->columns) {
      res = FAILURE;
    } else {
      for (int i = 0; i < A->rows && res != FAILURE; i++)
        for (int j = 0; j < A->columns && res != FAILURE; j++)
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) > ACC) res = FAILURE;
    }
  }
  return res;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  if (s21_check_matrix(A) || s21_check_matrix(B) || !result) {
    res = 1;
  } else {
    if (A->rows != B->rows || A->columns != B->columns) {
      res = 2;
    } else {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++)
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }
  return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  if (s21_check_matrix(A) || s21_check_matrix(B) || !result) {
    res = 1;
  } else {
    if (A->rows != B->rows || A->columns != B->columns) {
      res = 2;
    } else {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++)
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }
  return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int res = 0;
  if (s21_check_matrix(A) || !result) {
    res = 1;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[i][j] = A->matrix[i][j] * number;
  }
  return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  double tmp;
  if (s21_check_matrix(A) || s21_check_matrix(B) || !result) {
    res = 1;
  } else {
    if (A->columns != B->rows) {
      res = 2;
    } else {
      s21_create_matrix(A->rows, B->columns, result);
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < B->columns; j++) {
          tmp = 0;
          for (int k = 0; k < A->columns; k++)
            tmp += A->matrix[i][k] * B->matrix[k][j];
          result->matrix[i][j] = tmp;
        }
    }
  }
  return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int res = 0;
  if (s21_check_matrix(A) || !result) {
    res = 1;
  } else {
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[j][i] = A->matrix[i][j];
  }
  return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int res = 0;
  matrix_t matrix_tmp = {0};
  double det_temp = 0;
  if (s21_check_matrix(A) || !result) {
    res = 1;
  } else {
    if (A->columns != A->rows) {
      res = 2;
    } else {
      s21_create_matrix(A->columns, A->rows, result);
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++) {
          s21_minor(A, j, i, &matrix_tmp);
          s21_determinant(&matrix_tmp, &det_temp);
          s21_remove_matrix(&matrix_tmp);
          if ((i + j) % 2 == 1) det_temp = -det_temp;
          result->matrix[i][j] = det_temp;
        }
    }
  }
  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  int res = 0;
  matrix_t matrix_tmp = {0};
  double det_temp = 0;
  if (s21_check_matrix(A) || !result) {
    res = 1;
  } else {
    if (A->columns != A->rows) {
      res = 2;
    } else {
      if (A->columns == 1) {
        *result = A->matrix[0][0];
      } else if (A->columns == 2) {
        *result = A->matrix[0][0] * A->matrix[1][1] -
                  A->matrix[0][1] * A->matrix[1][0];
      } else {
        *result = 0;
        for (int i = 0; i < A->columns; i++) {
          s21_minor(A, i, 0, &matrix_tmp);
          s21_determinant(&matrix_tmp, &det_temp);
          s21_remove_matrix(&matrix_tmp);
          if (i % 2 == 1) det_temp = -det_temp;
          *result = *result + A->matrix[0][i] * det_temp;
        }
      }
    }
  }
  return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = 0, res_tmp = 0;
  matrix_t matrix_complements = {0}, matrix_complements_t = {0};
  double det = 0;
  res_tmp = s21_determinant(A, &det);
  if (res_tmp == 1 || !result) {
    res = 1;
  } else if (res_tmp == 2 || fabs(det) < ACC) {
    res = 2;
  } else {
    if (A->columns == 1) {
      s21_create_matrix(1, 1, result);
      result->matrix[0][0] = 1 / A->matrix[0][0];
    } else {
      s21_calc_complements(A, &matrix_complements);
      s21_transpose(&matrix_complements, &matrix_complements_t);
      s21_mult_number(&matrix_complements_t, 1 / det, result);
      s21_remove_matrix(&matrix_complements);
      s21_remove_matrix(&matrix_complements_t);
    }
  }
  return res;
}

void s21_minor(matrix_t *A, int col, int row, matrix_t *result) {
  int newi, newj;
  s21_create_matrix(A->columns - 1, A->rows - 1, result);
  for (int i = 0; i < A->rows - 1; i++) {
    if (i >= row)
      newi = i + 1;
    else
      newi = i;
    for (int j = 0; j < A->columns - 1; j++) {
      if (j >= col)
        newj = j + 1;
      else
        newj = j;
      result->matrix[i][j] = A->matrix[newi][newj];
    }
  }
}

int s21_check_matrix(matrix_t *A) {
  int res = 0;
  if (!A) {
    res = 1;
  } else {
    if (A->matrix == NULL || A->columns <= 0 || A->rows <= 0) res = 1;
  }
  return res;
}
