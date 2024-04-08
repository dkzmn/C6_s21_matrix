#include <check.h>
#include <stdio.h>
#include <stdlib.h>

#include "s21_matrix.h"

START_TEST(create_1) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = 0;
      ck_assert_ldouble_eq_tol(0, m.matrix[i][j], ACC);
    }
  }
  ck_assert_int_eq(m.rows, rows);
  ck_assert_int_eq(m.columns, cols);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(create_2) {
  const int rows = 0;
  const int cols = 10;
  matrix_t m = {0};
  ck_assert_int_eq(s21_create_matrix(rows, cols, &m), 1);
}
END_TEST

START_TEST(create_3) {
  const int rows = 10;
  const int cols = 0;
  matrix_t m = {0};
  ck_assert_int_eq(s21_create_matrix(rows, cols, &m), 1);
}
END_TEST

START_TEST(eq_1) {
  const int rows = (int)rand() % 10 + 1;
  const int cols = (int)rand() % 10 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = (double)rand() / rand();
      m.matrix[i][j] = rand_val;
      mtx.matrix[i][j] = rand_val;
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&m, &mtx), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
}
END_TEST

START_TEST(eq_2) {
  const int rows = (int)rand() % 10 + 1;
  const int cols = (int)rand() % 10 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = (double)rand() / rand();
      m.matrix[i][j] = rand_val + 1;
      mtx.matrix[i][j] = rand_val;
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&m, &mtx), FAILURE);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
}
END_TEST

START_TEST(eq_3) {
  const int rows = (int)rand() % 10 + 1;
  const int cols = (int)rand() % 10 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols + 1, &mtx);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = (double)rand() / rand();
      m.matrix[i][j] = rand_val + 1;
      mtx.matrix[i][j] = rand_val;
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&m, &mtx), FAILURE);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
}
END_TEST

START_TEST(eq_4) {
  const int rows = (int)rand() % 10 + 1;
  const int cols = (int)rand() % 10 + 1;
  matrix_t m = {0};
  matrix_t mtx = {0};
  m.matrix = NULL;
  mtx.matrix = NULL;
  m.columns = cols;
  mtx.columns = cols;
  m.rows = rows;
  mtx.rows = rows;
  ck_assert_int_eq(s21_eq_matrix(&m, &mtx), FAILURE);
  ck_assert_int_eq(s21_eq_matrix(NULL, &mtx), FAILURE);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
}
END_TEST

START_TEST(sum_1) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);
  matrix_t check = {0};
  s21_create_matrix(rows, cols, &check);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = (double)rand() / rand();
      mtx.matrix[i][j] = (double)rand() / rand();
      check.matrix[i][j] = m.matrix[i][j] + mtx.matrix[i][j];
    }
  }
  matrix_t res = {0};
  ck_assert_int_eq(s21_sum_matrix(&m, &mtx, &res), 0);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(sum_2) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols + 1, &mtx);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = (double)rand() / rand();
      mtx.matrix[i][j] = (double)rand() / rand();
    }
  }
  matrix_t res = {0};
  ck_assert_int_eq(s21_sum_matrix(&m, &mtx, &res), 2);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(sum_3) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  m.matrix = NULL;
  m.columns = cols;
  m.rows = rows;
  matrix_t mtx = {0};
  mtx.matrix = NULL;
  mtx.columns = cols;
  mtx.rows = rows;
  matrix_t res = {0};
  ck_assert_int_eq(s21_sum_matrix(&m, &mtx, &res), 1);
  ck_assert_int_eq(s21_sum_matrix(&m, NULL, &res), 1);
  ck_assert_int_eq(s21_sum_matrix(&m, &mtx, NULL), 1);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(sub_1) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols, &mtx);
  matrix_t check = {0};
  s21_create_matrix(rows, cols, &check);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = (double)rand() / rand();
      mtx.matrix[i][j] = (double)rand() / rand();
      check.matrix[i][j] = m.matrix[i][j] - mtx.matrix[i][j];
    }
  }
  matrix_t res = {0};
  ck_assert_int_eq(s21_sub_matrix(&m, &mtx, &res), 0);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(sub_2) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(rows, cols + 1, &mtx);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = (double)rand() / rand();
      mtx.matrix[i][j] = (double)rand() / rand();
    }
  }
  matrix_t res = {0};
  ck_assert_int_eq(s21_sub_matrix(&m, &mtx, &res), 2);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(sub_3) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  m.matrix = NULL;
  m.columns = cols;
  m.rows = rows;
  matrix_t mtx = {0};
  mtx.matrix = NULL;
  mtx.columns = cols;
  mtx.rows = rows;
  matrix_t res = {0};
  ck_assert_int_eq(s21_sub_matrix(&m, &mtx, &res), 1);
  ck_assert_int_eq(s21_sub_matrix(&m, NULL, &res), 1);
  ck_assert_int_eq(s21_sub_matrix(&m, &mtx, NULL), 1);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(mult_number_1) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t check = {0};
  s21_create_matrix(rows, cols, &check);
  double mult_number = (double)rand() / rand();
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      m.matrix[i][j] = (double)rand() / rand();
      check.matrix[i][j] = m.matrix[i][j] * mult_number;
    }
  }
  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_number(&m, mult_number, &res), 0);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(mult_number_2) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  m.matrix = NULL;
  m.columns = cols;
  m.rows = rows;
  matrix_t mtx = {0};
  mtx.matrix = NULL;
  mtx.columns = cols;
  mtx.rows = rows;
  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_number(&m, 1, &res), 1);
  ck_assert_int_eq(s21_mult_number(NULL, 1, &res), 1);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(mult_matrix_1) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(cols, rows, &mtx);
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++) m.matrix[i][j] = (double)rand() / rand();
  for (int i = 0; i < cols; i++)
    for (int j = 0; j < rows; j++) mtx.matrix[i][j] = (double)rand() / rand();
  matrix_t check = {0};
  s21_create_matrix(m.rows, mtx.columns, &check);
  for (int i = 0; i < m.rows; i++) {
    for (int j = 0; j < mtx.columns; j++) {
      check.matrix[i][j] = 0;
      for (int k = 0; k < m.columns; k++)
        check.matrix[i][j] += m.matrix[i][k] * mtx.matrix[k][j];
    }
  }
  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_matrix(&m, &mtx, &res), 0);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(mult_matrix_2) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(cols + 1, rows, &mtx);
  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_matrix(&m, &mtx, &res), 2);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(mult_matrix_3) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  m.matrix = NULL;
  m.columns = cols;
  m.rows = rows;
  matrix_t mtx = {0};
  mtx.matrix = NULL;
  mtx.columns = cols;
  mtx.rows = rows;
  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_matrix(&m, &mtx, &res), 1);
  ck_assert_int_eq(s21_mult_matrix(NULL, &mtx, &res), 1);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(mult_matrix_4) {
  const int rows = 2;
  const int cols = 3;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t mtx = {0};
  s21_create_matrix(cols, rows, &mtx);
  for (int i = 0, c = 1; i < rows; i++)
    for (int j = 0; j < cols; j++) m.matrix[i][j] = c++;
  for (int i = 0, c = 7; i < cols; i++)
    for (int j = 0; j < rows; j++) mtx.matrix[i][j] = c++;
  matrix_t check = {0};
  s21_create_matrix(m.rows, mtx.columns, &check);
  check.matrix[0][0] = 58;
  check.matrix[0][1] = 64;
  check.matrix[1][0] = 139;
  check.matrix[1][1] = 154;
  matrix_t res = {0};
  ck_assert_int_eq(s21_mult_matrix(&m, &mtx, &res), 0);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&mtx);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(transpose_1) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  s21_create_matrix(rows, cols, &m);
  matrix_t check = {0};
  s21_create_matrix(cols, rows, &check);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = (double)rand() / rand();
      m.matrix[i][j] = rand_val;
      check.matrix[j][i] = rand_val;
    }
  }
  matrix_t res = {0};
  ck_assert_int_eq(s21_transpose(&m, &res), 0);
  ck_assert_int_eq(s21_eq_matrix(&check, &res), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&res);
  s21_remove_matrix(&check);
}
END_TEST

START_TEST(transpose_2) {
  const int rows = (int)rand() % 100 + 1;
  const int cols = (int)rand() % 100 + 1;
  matrix_t m = {0};
  m.matrix = NULL;
  m.columns = cols;
  m.rows = rows;
  matrix_t res = {0};
  ck_assert_int_eq(s21_transpose(&m, &res), 1);
  ck_assert_int_eq(s21_transpose(NULL, &res), 1);
  ck_assert_int_eq(s21_transpose(&m, NULL), 1);
  s21_remove_matrix(&m);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(determinant_1) {
  const int size = 5;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) m.matrix[i][j] = j;
  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 0, ACC);
  ck_assert_int_eq(code, 0);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant_2) {
  const int size = 4;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) m.matrix[i][j] = j + i;
  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 0, ACC);
  ck_assert_int_eq(code, 0);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant_3) {
  const int size = 5;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = 0;
  m.matrix[0][1] = 6;
  m.matrix[0][2] = -2;
  m.matrix[0][3] = -1;
  m.matrix[0][4] = 5;
  m.matrix[1][0] = 0;
  m.matrix[1][1] = 0;
  m.matrix[1][2] = 0;
  m.matrix[1][3] = -9;
  m.matrix[1][4] = -7;
  m.matrix[2][0] = 0;
  m.matrix[2][1] = 15;
  m.matrix[2][2] = 35;
  m.matrix[2][3] = 0;
  m.matrix[2][4] = 0;
  m.matrix[3][0] = 0;
  m.matrix[3][1] = -1;
  m.matrix[3][2] = -11;
  m.matrix[3][3] = -2;
  m.matrix[3][4] = 1;
  m.matrix[4][0] = -2;
  m.matrix[4][1] = -2;
  m.matrix[4][2] = 3;
  m.matrix[4][3] = 0;
  m.matrix[4][4] = -2;
  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 2480, ACC);
  ck_assert_int_eq(code, 0);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant_4) {
  const int size = 3;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = 2;
  m.matrix[0][1] = 3;
  m.matrix[0][2] = 1;
  m.matrix[1][0] = 7;
  m.matrix[1][1] = 4;
  m.matrix[1][2] = 1;
  m.matrix[2][0] = 9;
  m.matrix[2][1] = -2;
  m.matrix[2][2] = 1;
  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, -32, ACC);
  ck_assert_int_eq(code, 0);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant_5) {
  const int size = 2;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = -5;
  m.matrix[0][1] = -4;
  m.matrix[1][0] = -2;
  m.matrix[1][1] = -3;
  double res = 0;
  int code = s21_determinant(&m, &res);
  ck_assert_double_eq_tol(res, 7, ACC);
  ck_assert_int_eq(code, 0);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(determinant_6) {
  matrix_t m = {0};
  m.matrix = NULL;
  double res = 0;
  ck_assert_int_eq(s21_determinant(&m, &res), 1);
  ck_assert_int_eq(s21_determinant(NULL, &res), 1);
  ck_assert_int_eq(s21_determinant(&m, NULL), 1);
}
END_TEST

START_TEST(complements_1) {
  matrix_t a = {0};
  s21_create_matrix(3, 3, &a);
  matrix_t res = {0};
  s21_create_matrix(3, 3, &res);
  matrix_t m;
  res.matrix[0][0] = 0;
  res.matrix[0][1] = 10;
  res.matrix[0][2] = -20;
  res.matrix[1][0] = 4;
  res.matrix[1][1] = -14;
  res.matrix[1][2] = 8;
  res.matrix[2][0] = -8;
  res.matrix[2][1] = -2;
  res.matrix[2][2] = 4;
  a.matrix[0][0] = 1;
  a.matrix[0][1] = 2;
  a.matrix[0][2] = 3;
  a.matrix[1][0] = 0;
  a.matrix[1][1] = 4;
  a.matrix[1][2] = 2;
  a.matrix[2][0] = 5;
  a.matrix[2][1] = 2;
  a.matrix[2][2] = 1;
  ck_assert_int_eq(s21_calc_complements(&a, &m), 0);
  ck_assert_int_eq(s21_eq_matrix(&res, &m), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&a);
  s21_remove_matrix(&res);
}
END_TEST

START_TEST(complements_2) {
  matrix_t A, B, C;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &C);
  A.matrix[0][0] = 1.0;
  A.matrix[0][1] = 2.0;
  A.matrix[0][2] = 3.0;
  A.matrix[1][0] = 0.0;
  A.matrix[1][1] = 4.0;
  A.matrix[1][2] = 2.0;
  A.matrix[2][0] = 5.0;
  A.matrix[2][1] = 2.0;
  A.matrix[2][2] = 1.0;
  ck_assert_int_eq(s21_calc_complements(&A, &B), 0);
  C.matrix[0][0] = 0.0;
  C.matrix[0][1] = 10.0;
  C.matrix[0][2] = -20.0;
  C.matrix[1][0] = 4.0;
  C.matrix[1][1] = -14.0;
  C.matrix[1][2] = 8.0;
  C.matrix[2][0] = -8.0;
  C.matrix[2][1] = -2.0;
  C.matrix[2][2] = 4.0;
  int res = s21_eq_matrix(&B, &C);
  ck_assert_int_eq(res, SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(complements_3) {
  matrix_t a = {0};
  matrix_t b = {0};
  ck_assert_int_eq(s21_calc_complements(&a, &b), 1);
}
END_TEST

START_TEST(complements_4) {
  matrix_t a = {0};
  matrix_t b = {0};
  a.matrix = NULL;
  ck_assert_int_eq(s21_calc_complements(&a, &b), 1);
  ck_assert_int_eq(s21_calc_complements(NULL, &b), 1);
  ck_assert_int_eq(s21_calc_complements(&a, NULL), 1);
  s21_create_matrix(1, 2, &a);
  ck_assert_int_eq(s21_calc_complements(&a, &b), 2);
  s21_remove_matrix(&a);
}
END_TEST

START_TEST(complements_5) {
  const int size = 3;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = 1;
  m.matrix[0][1] = 2;
  m.matrix[0][2] = 3;
  m.matrix[1][0] = 0;
  m.matrix[1][1] = 4;
  m.matrix[1][2] = 2;
  m.matrix[2][0] = 5;
  m.matrix[2][1] = 2;
  m.matrix[2][2] = 1;
  matrix_t expected = {0};
  s21_create_matrix(size, size, &expected);
  expected.matrix[0][0] = 0;
  expected.matrix[0][1] = 10;
  expected.matrix[0][2] = -20;
  expected.matrix[1][0] = 4;
  expected.matrix[1][1] = -14;
  expected.matrix[1][2] = 8;
  expected.matrix[2][0] = -8;
  expected.matrix[2][1] = -2;
  expected.matrix[2][2] = 4;
  matrix_t res = {0};
  ck_assert_int_eq(s21_calc_complements(&m, &res), 0);
  ck_assert_int_eq(s21_eq_matrix(&expected, &res), SUCCESS);
  s21_remove_matrix(&m);
  s21_remove_matrix(&res);
  s21_remove_matrix(&expected);
}
END_TEST

START_TEST(inverse_1) {
  matrix_t A, B, C;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &C);
  C.matrix[0][0] = 44300.0 / 367429.0;
  C.matrix[0][1] = -236300.0 / 367429.0;
  C.matrix[0][2] = 200360.0 / 367429.0;
  C.matrix[1][0] = 20600.0 / 367429.0;
  C.matrix[1][1] = 56000.0 / 367429.0;
  C.matrix[1][2] = -156483.0 / 367429.0;
  C.matrix[2][0] = 30900.0 / 367429.0;
  C.matrix[2][1] = 84000.0 / 367429.0;
  C.matrix[2][2] = -51010.0 / 367429.0;
  A.matrix[0][0] = 2.8;
  A.matrix[0][1] = 1.3;
  A.matrix[0][2] = 7.01;
  A.matrix[1][0] = -1.03;
  A.matrix[1][1] = -2.3;
  A.matrix[1][2] = 3.01;
  A.matrix[2][0] = 0;
  A.matrix[2][1] = -3;
  A.matrix[2][2] = 2;
  s21_inverse_matrix(&A, &B);
  int res = s21_eq_matrix(&B, &C);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}

END_TEST

START_TEST(inverse_2) {
  matrix_t A, C;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &C);
  C.matrix[0][0] = 1.0;
  C.matrix[0][1] = -1.0;
  C.matrix[0][2] = 1.0;
  C.matrix[1][0] = -38.0;
  C.matrix[1][1] = 41.0;
  C.matrix[1][2] = -34.0;
  C.matrix[2][0] = 27.0;
  C.matrix[2][1] = -29.0;
  C.matrix[2][2] = 24.0;
  A.matrix[0][0] = 2.0;
  A.matrix[0][1] = 5.0;
  A.matrix[0][2] = 7.0;
  A.matrix[1][0] = 6.0;
  A.matrix[1][1] = 3.0;
  A.matrix[1][2] = 4.0;
  A.matrix[2][0] = 5.0;
  A.matrix[2][1] = -2.0;
  A.matrix[2][2] = -3.0;
  matrix_t B;
  s21_inverse_matrix(&A, &B);
  int res = s21_eq_matrix(&B, &C);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&C);
}
END_TEST

START_TEST(inverse_3) {
  matrix_t m = {0};
  matrix_t result = {0};
  int codec = s21_create_matrix(1, 1, &m);
  if (!codec) {
    m.matrix[0][0] = 1431.12312331;
    int code = s21_inverse_matrix(&m, &result);
    ck_assert_int_eq(result.matrix[0][0] == (1.0 / m.matrix[0][0]), 1);
    ck_assert_int_eq(code, 0);
    s21_remove_matrix(&m);
    s21_remove_matrix(&result);
  }
}
END_TEST

START_TEST(inverse_4) {
  matrix_t m = {0};
  matrix_t result = {0};
  int codec = s21_create_matrix(1, 1, &m);
  if (!codec) {
    m.matrix[0][0] = 0;
    int code = s21_inverse_matrix(&m, &result);
    ck_assert_int_eq(code, 2);
    s21_remove_matrix(&m);
  }
}
END_TEST

START_TEST(inverse_5) {
  matrix_t m = {0};
  matrix_t result = {0};
  ck_assert_int_eq(s21_inverse_matrix(&m, &result), 1);
  ck_assert_int_eq(s21_inverse_matrix(NULL, &result), 1);
  ck_assert_int_eq(s21_inverse_matrix(&m, NULL), 1);
}
END_TEST

START_TEST(inverse_6) {
  matrix_t m = {0};
  matrix_t result = {0};
  int codec = s21_create_matrix(1, 4, &m);
  if (!codec) {
    int code = s21_inverse_matrix(&m, &result);
    ck_assert_int_eq(code, 2);
    s21_remove_matrix(&m);
  }
}
END_TEST

START_TEST(inverse_7) {
  matrix_t m = {0};
  matrix_t expected = {0};
  int codec1, codec2;
  codec1 = s21_create_matrix(3, 3, &m);
  if (!codec1) codec2 = s21_create_matrix(3, 3, &expected);
  if (!codec1 && !codec2) {
    m.matrix[0][0] = 2;
    m.matrix[0][1] = 5;
    m.matrix[0][2] = 7;
    m.matrix[1][0] = 6;
    m.matrix[1][1] = 3;
    m.matrix[1][2] = 4;
    m.matrix[2][0] = 5;
    m.matrix[2][1] = -2;
    m.matrix[2][2] = -3;
    expected.matrix[0][0] = 1;
    expected.matrix[0][1] = -1;
    expected.matrix[0][2] = 1;
    expected.matrix[1][0] = -38;
    expected.matrix[1][1] = 41;
    expected.matrix[1][2] = -34;
    expected.matrix[2][0] = 27;
    expected.matrix[2][1] = -29;
    expected.matrix[2][2] = 24;
    matrix_t result = {0};
    int code = s21_inverse_matrix(&m, &result);
    ck_assert_int_eq(s21_eq_matrix(&result, &expected), SUCCESS);
    ck_assert_int_eq(code, 0);
    s21_remove_matrix(&m);
    s21_remove_matrix(&result);
    s21_remove_matrix(&expected);
  }
}
END_TEST

START_TEST(inverse_8) {
  const int size = 3;
  matrix_t m = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = 2;
  m.matrix[0][1] = 5;
  m.matrix[0][2] = 7;
  m.matrix[1][0] = 6;
  m.matrix[1][1] = 3;
  m.matrix[1][2] = 4;
  m.matrix[2][0] = 5;
  m.matrix[2][1] = -2;
  m.matrix[2][2] = -3;
  matrix_t res = {0};
  s21_inverse_matrix(&m, &res);
  matrix_t expected = {0};
  s21_create_matrix(size, size, &expected);
  expected.matrix[0][0] = 1;
  expected.matrix[0][1] = -1;
  expected.matrix[0][2] = 1;
  expected.matrix[1][0] = -38;
  expected.matrix[1][1] = 41;
  expected.matrix[1][2] = -34;
  expected.matrix[2][0] = 27;
  expected.matrix[2][1] = -29;
  expected.matrix[2][2] = 24;
  ck_assert_int_eq(s21_eq_matrix(&expected, &res), SUCCESS);
  s21_remove_matrix(&expected);
  s21_remove_matrix(&res);
  s21_remove_matrix(&m);
}
END_TEST

START_TEST(complex_1) {
  const int size = 3;
  matrix_t m = {0}, inv = {0}, res = {0}, expected = {0};
  s21_create_matrix(size, size, &m);
  m.matrix[0][0] = 1;
  m.matrix[0][1] = 2;
  m.matrix[0][2] = 3;
  m.matrix[1][0] = 4;
  m.matrix[1][1] = 5;
  m.matrix[1][2] = 6;
  m.matrix[2][0] = 7;
  m.matrix[2][1] = 8;
  m.matrix[2][2] = 10;
  s21_inverse_matrix(&m, &inv);
  s21_mult_matrix(&m, &inv, &res);
  s21_create_matrix(size, size, &expected);
  expected.matrix[0][0] = 1;
  expected.matrix[0][1] = 0;
  expected.matrix[0][2] = 0;
  expected.matrix[1][0] = 0;
  expected.matrix[1][1] = 1;
  expected.matrix[1][2] = 0;
  expected.matrix[2][0] = 0;
  expected.matrix[2][1] = 0;
  expected.matrix[2][2] = 1;
  ck_assert_int_eq(s21_eq_matrix(&expected, &res), SUCCESS);
  s21_remove_matrix(&expected);
  s21_remove_matrix(&res);
  s21_remove_matrix(&inv);
  s21_remove_matrix(&m);
}
END_TEST

Suite *s21_matrix_suite(void) {
  Suite *suite = suite_create("s21_matrix");
  TCase *tcase_core_create = tcase_create("create");
  tcase_add_loop_test(tcase_core_create, create_1, 0, 100);
  tcase_add_test(tcase_core_create, create_2);
  tcase_add_test(tcase_core_create, create_3);
  suite_add_tcase(suite, tcase_core_create);

  TCase *tcase_core_eq = tcase_create("eq");
  tcase_add_loop_test(tcase_core_eq, eq_1, 0, 100);
  tcase_add_test(tcase_core_eq, eq_2);
  tcase_add_test(tcase_core_eq, eq_3);
  tcase_add_test(tcase_core_eq, eq_4);
  suite_add_tcase(suite, tcase_core_eq);

  TCase *tcase_core_sum = tcase_create("sum");
  tcase_add_loop_test(tcase_core_sum, sum_1, 0, 100);
  tcase_add_test(tcase_core_sum, sum_2);
  tcase_add_test(tcase_core_sum, sum_3);
  suite_add_tcase(suite, tcase_core_sum);

  TCase *tcase_core_sub = tcase_create("sub");
  tcase_add_loop_test(tcase_core_sub, sub_1, 0, 100);
  tcase_add_test(tcase_core_sub, sub_2);
  tcase_add_test(tcase_core_sub, sub_3);
  suite_add_tcase(suite, tcase_core_sub);

  TCase *tcase_core_mult_number = tcase_create("mult_number");
  tcase_add_loop_test(tcase_core_mult_number, mult_number_1, 0, 100);
  tcase_add_test(tcase_core_mult_number, mult_number_2);
  suite_add_tcase(suite, tcase_core_mult_number);

  TCase *tcase_core_mult = tcase_create("mult");
  tcase_add_loop_test(tcase_core_mult, mult_matrix_1, 0, 100);
  tcase_add_test(tcase_core_mult, mult_matrix_2);
  tcase_add_test(tcase_core_mult, mult_matrix_3);
  tcase_add_test(tcase_core_mult, mult_matrix_4);
  suite_add_tcase(suite, tcase_core_mult);

  TCase *tcase_core_trans = tcase_create("transpose_matrix");
  tcase_add_loop_test(tcase_core_trans, transpose_1, 0, 100);
  tcase_add_test(tcase_core_trans, transpose_2);
  suite_add_tcase(suite, tcase_core_trans);

  TCase *tcase_core_determinant = tcase_create("determinant");
  tcase_add_test(tcase_core_determinant, determinant_1);
  tcase_add_test(tcase_core_determinant, determinant_2);
  tcase_add_test(tcase_core_determinant, determinant_3);
  tcase_add_test(tcase_core_determinant, determinant_4);
  tcase_add_test(tcase_core_determinant, determinant_5);
  tcase_add_test(tcase_core_determinant, determinant_6);
  suite_add_tcase(suite, tcase_core_determinant);

  TCase *tcase_core_complements = tcase_create("complements");
  tcase_add_test(tcase_core_complements, complements_1);
  tcase_add_test(tcase_core_complements, complements_2);
  tcase_add_test(tcase_core_complements, complements_3);
  tcase_add_test(tcase_core_complements, complements_4);
  tcase_add_test(tcase_core_complements, complements_5);
  suite_add_tcase(suite, tcase_core_complements);

  TCase *tcase_core_inverse = tcase_create("inverse");
  tcase_add_test(tcase_core_inverse, inverse_1);
  tcase_add_test(tcase_core_inverse, inverse_2);
  tcase_add_test(tcase_core_inverse, inverse_3);
  tcase_add_test(tcase_core_inverse, inverse_4);
  tcase_add_test(tcase_core_inverse, inverse_5);
  tcase_add_test(tcase_core_inverse, inverse_6);
  tcase_add_test(tcase_core_inverse, inverse_7);
  tcase_add_test(tcase_core_inverse, inverse_8);
  suite_add_tcase(suite, tcase_core_inverse);

  TCase *tcase_core_complex = tcase_create("complex");
  tcase_add_test(tcase_core_complex, complex_1);
  suite_add_tcase(suite, tcase_core_complex);

  return suite;
}

int main(void) {
  int failed_count;
  Suite *suite = s21_matrix_suite();
  SRunner *suite_runner = srunner_create(suite);
  srunner_set_fork_status(suite_runner, CK_NOFORK);
  srunner_run_all(suite_runner, CK_NORMAL);
  failed_count = srunner_ntests_failed(suite_runner);
  srunner_free(suite_runner);
  return (failed_count == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
