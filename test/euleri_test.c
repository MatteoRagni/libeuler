#include <stdio.h>
#include <math.h>
#include "libeuler.h"

void f(double *f, double t, const double *x, const double *u, const double **p, void *data)
{
  double A1 = 0.180;
  double k = 0.003;
  double a1 = 0.006;
  double g = 9.810;
  double A2 = 0.080;
  double a2 = 0.008;

  f[0] = 1.0 / A1 * (k * u[0] - a1 * sqrt(2 * g * x[0]));
  f[1] = 1.0 / A2 * (a1 * sqrt(2 * g * x[0]) - a2 * sqrt(2 * g * x[1]));
}

void df(double *df, double t, const double *x, const double *u, const double **p, void *data)
{
  double A1 = 0.180;
  double k = 0.003;
  double a1 = 0.006;
  double g = 9.810;
  double A2 = 0.080;
  double a2 = 0.008;

  df[0] = -(a1 * sqrt(g)) / ( A1 * sqrt(2 * x[0]));
  df[1] = (a1 * sqrt(g)) / (A2 * sqrt(2 * x[0]));
  df[2] = 0;
  df[3] = -(a2 * sqrt(g)) / ( A2 * sqrt(2 * x[1]));
}

euler_options opt = {
    .ts = 1e-2,
    .alpha = 0.5,
    .x_size = 2,
    .u_offset = 0,
    .ordering = LAPACK_COL_MAJOR,
    .s_tol = 1e-12,
    .x_tol = 1e-12,
    .max_iter = 100,
    .f = f,
    .df = df,
    .data = NULL};

double input(double t)
{
  if (t < 251)
    return 10.0;
  if (t < 451)
    return 5.0;
  return 8.0;
}

int main()
{
  double t = 0;
  double x[2] = {1e-6, 0.1};
  double xp[2] = {0, 0.1};
  double u = input(t);

  while (t < 500)
  {
    printf("% 5.3f, % 5.6f, % 5.6f, % 5.6f\n", t, u, x[0], x[1]);

    euler(&opt, xp, t, x, &u, NULL, NULL);

    t += opt.ts;
    u = input(t);
    x[0] = xp[0];
    x[1] = xp[1];
  }
}