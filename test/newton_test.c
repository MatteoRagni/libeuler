#include <stdio.h>
#include <math.h>
#include "libnewton.h"

void function(double *f, const double *x, const double *u, const double **p, void *data) {
  f[0] = 2 * x[0] - x[1] - exp(-x[0]);
  f[1] = -x[0] + 2 * x[1] - exp(-x[1]);
}

void gradient(double *df, const double *x, const double *u, const double **p, void *data) {
  df[0] = 2 + exp(-x[0]);
  df[1] = -1;
  df[2] = -1;
  df[3] = 2 + exp(-x[1]);
}


newton_options options = {
  LAPACK_COL_MAJOR,
  2,
  2,
  1e-12,
  1e-12,
  100,
  function,
  gradient
};

int main() {

  double x[2] = {10, 10};
  double f[2] = {0, 0};

  newton_ret ret = newton_solve(&options, x, NULL, NULL, NULL);
  function(f, x, NULL, NULL, NULL);

  printf("EXIT = %d\n", ret);
  printf("  f(% 5.10f, % 5.10f) = (% 5.10f, % 5.10f)\n", x[0], x[1], f[0], f[1]);
  printf("  |f| = % 5.20f\n", options.f_tol);
  printf("  |x| = % 5.20f\n", options.x_tol);
  printf(" iter = %d\n", options.max_iter);

  return 0;
}
