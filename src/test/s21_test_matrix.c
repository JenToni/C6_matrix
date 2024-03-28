#include "s21_test_matrix.h"

int main() {
  int failed = 0;
  int passed = 0;
  Suite *matrix_tests[] = {suite_create_matrix(),
                           suite_eq_matrix(),
                           suite_sub_matrix(),
                           suite_sum_matrix(),
                           suite_mult_matrix(),
                           suite_mult_number_matrix(),
                           suite_transpose_matrix(),
                           suite_calc_complements(),
                           suite_determinant(),
                           suite_inverse_matrix(),
                           NULL};

  for (int i = 0; matrix_tests[i] != NULL; i++) {
    SRunner *sr = srunner_create(matrix_tests[i]);

    srunner_set_fork_status(sr, CK_NOFORK);
    srunner_run_all(sr, CK_NORMAL);

    failed += srunner_ntests_failed(sr);
    passed += srunner_ntests_run(sr) - srunner_ntests_failed(sr);
    srunner_free(sr);
  }

  printf("========= PASSED: %d =========\n", passed);
  printf("========= FAILED: %d =========\n", failed);

  return failed == 0 ? 0 : 1;
}
