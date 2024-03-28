#ifndef SRC_UNIT_TESTS_S21_MATRIX_H_
#define SRC_UNIT_TESTS_S21_MATRIX_H_

#include <check.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "../s21_matrix.h"

double get_rand(double min, double max);

Suite *suite_create_matrix(void);
Suite *suite_transpose_matrix(void);
Suite *suite_eq_matrix(void);
Suite *suite_sum_matrix(void);
Suite *suite_sub_matrix(void);
Suite *suite_mult_matrix(void);
Suite *suite_mult_number_matrix(void);
Suite *suite_calc_complements(void);
Suite *suite_determinant(void);
Suite *suite_inverse_matrix(void);

#endif  // SRC_UNIT_TESTS_S21_MATRIX_H_