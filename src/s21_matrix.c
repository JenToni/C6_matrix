#include "s21_matrix.h"

// Создание матрицы
int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error = OK;
  if (rows < 1 || columns < 1)
    error = INCORRECT_MATRIX;
  else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    for (int i = 0; i < rows; i++)
      result->matrix[i] = (double *)calloc(columns, sizeof(double));
  }
  return error;
}

// Удаление матрицы
void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
    A->columns = 0;
    A->rows = 0;
    A->matrix = NULL;
  }
}

// Сравнение матриц до 7 знака
int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;
  if (A->columns == B->columns && A->rows == B->rows &&
      (s21_is_Emty(A) == OK && s21_is_Emty(B) == OK)) {
    for (int i = 0; i < A->rows && res; i++)
      for (int j = 0; j < A->columns && res; j++)
        if (round(A->matrix[i][j] * pow(10, 7)) !=
            round(B->matrix[i][j] * pow(10, 7)))
          res = FAILURE;
  } else
    res = FAILURE;
  return res;
}

// Сумма матриц
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (s21_is_Emty(A) == OK && s21_is_Emty(B) == OK) {
    if (A->columns == B->columns && A->rows == B->rows) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++)
          result->matrix[i][j] = (A->matrix[i][j] + B->matrix[i][j]);
    } else
      res = CALCULATION_ERROR;
  } else
    res = INCORRECT_MATRIX;
  return res;
}

// Разность матриц
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (s21_is_Emty(A) == OK && s21_is_Emty(B) == OK) {
    if (A->columns == B->columns && A->rows == B->rows) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++)
          result->matrix[i][j] = (A->matrix[i][j] - B->matrix[i][j]);
    } else
      res = CALCULATION_ERROR;
  } else
    res = INCORRECT_MATRIX;
  return res;
}

// Произведение матрицы и числа
int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int res = OK;
  if (s21_is_Emty(A) == OK) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[i][j] = (A->matrix[i][j] * number);
  } else
    res = INCORRECT_MATRIX;
  return res;
}

// Произведение матриц
// C(i,j) = A(i,1) × B(1,j) + A(i,2) × B(2,j) + … + A(i,k) × B(k,j).
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (s21_is_Emty(A) == OK && s21_is_Emty(B) == OK) {
    if (A->columns == B->rows) {
      s21_create_matrix(A->rows, B->columns, result);
      for (int i = 0; i < result->rows; i++)
        for (int j = 0; j < result->columns; j++)
          for (int k = 0; k < B->rows; k++)
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
    } else
      res = CALCULATION_ERROR;
  } else
    res = INCORRECT_MATRIX;
  return res;
}

// проверка матрицы на корректность
int s21_is_Emty(matrix_t *matrix) {
  int res = OK;
  if (matrix == NULL || matrix->matrix == NULL || matrix->rows <= 0 ||
      matrix->columns <= 0) {
    res = INCORRECT_MATRIX;
  } else {
    res = OK;
  }
  return res;
}

// Транспортировние матрицы
int s21_transpose(matrix_t *A, matrix_t *result) {
  int res = OK;
  if (s21_is_Emty(A) == OK) {
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[j][i] = A->matrix[i][j];
  } else
    res = INCORRECT_MATRIX;
  return res;
}

// Минор матрицы и матрица алгебраических дополнений
int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int res = OK;
  if (s21_is_Emty(A) == OK) {
    if (A->rows == A->columns && A->rows > 1 && A->columns > 1) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++) {
          matrix_t minor = {0};
          double determinant = 0;
          s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
          s21_copy_matrix(i, j, A, &minor);
          determinant = s21_get_determinant(&minor);
          result->matrix[i][j] = pow(-1, (i + j)) * determinant;
          s21_remove_matrix(&minor);
        }
    } else
      res = CALCULATION_ERROR;
  } else
    res = INCORRECT_MATRIX;
  return res;
}

// Определитель матрицы
int s21_determinant(matrix_t *A, double *result) {
  *result = 0.0;
  int res = OK;
  if (s21_is_Emty(A) == 0) {
    if (A->rows == A->columns) {
      *result = s21_get_determinant(A);
    } else {
      res = CALCULATION_ERROR;
    }
  } else {
    res = INCORRECT_MATRIX;
  }
  return res;
}

void s21_copy_matrix(int row, int col, matrix_t *A, matrix_t *result) {
  int r_row = 0;
  int r_col = 0;
  for (int i = 0; i < A->rows; i++) {
    if (i == row) continue;
    r_col = 0;
    for (int j = 0; j < A->columns; j++) {
      if (j == col) continue;
      result->matrix[r_row][r_col] = A->matrix[i][j];
      r_col++;
    }
    r_row++;
  }
}

double s21_get_determinant(matrix_t *matrix) {
  double determ = 0.0;
  if (matrix->rows == 1)
    determ = matrix->matrix[0][0];
  else {
    matrix_t minor = {0};
    s21_create_matrix(matrix->rows - 1, matrix->columns - 1, &minor);
    for (int i = 0; i < matrix->columns; i++) {
      s21_copy_matrix(0, i, matrix, &minor);
      if (i % 2)
        determ -= matrix->matrix[0][i] * s21_get_determinant(&minor);
      else
        determ += matrix->matrix[0][i] * s21_get_determinant(&minor);
    }
    s21_remove_matrix(&minor);
  }
  return determ;
}

// Обратная матрица
int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = INCORRECT_MATRIX;
  if (s21_is_Emty(A) == OK) {
    res = CALCULATION_ERROR;
    if (A->rows == 1 && A->columns == 1 && A->matrix[0][0] != 0) {
      res = s21_create_matrix(A->rows, A->rows, result);
      if (res == 0) result->matrix[0][0] = 1.0 / A->matrix[0][0];
      return res;
    }
    double det = 0.0;
    s21_determinant(A, &det);
    if (det != 0) {
      matrix_t calc = {0};
      res = s21_calc_complements(A, &calc);
      if (res == OK) {
        matrix_t trans = {0};
        res = s21_transpose(&calc, &trans);
        if (res == OK) {
          s21_mult_number(&trans, 1 / det, result);
        }
        s21_remove_matrix(&trans);
      }
      s21_remove_matrix(&calc);
    }
  }
  return res;
}
