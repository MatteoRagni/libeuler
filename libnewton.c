/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright 2018 - Matteo Ragni, Matteo Cocetti - University of Trento
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include "libnewton.h"


newton_ret newton_solve(newton_options *opt, const double t, double *x, const double *u, const double **p, void *data) {
  double set_f_tol = opt->f_tol;
  double set_x_tol = opt->x_tol;
  lapack_int counts = 0;

  newton_ret ret = NEWTON_GENERIC_ERROR;
  double *df, *f;
  
  lapack_int ldb = opt->f_size > opt->x_size ? opt->f_size : opt->x_size;
  lapack_int lda = opt->f_size;

  f = (double*)calloc(ldb, sizeof(double));
  if (!f)
    return NEWTON_MALLOC_ERROR;
  df = (double*)calloc(opt->f_size * opt->x_size, sizeof(double));
  if (!df) {
    free(f);
    return NEWTON_MALLOC_ERROR;
  }

  while (counts <= opt->max_iter) {
    opt->f(f, t, x, u, p, data);                                                  /* FUNCTION EVALUATION */
    opt->f_tol = cblas_dnrm2(opt->f_size, f, 1);

    /* Function Tollerance condition */
    if (opt->f_tol < set_f_tol) {
      ret = NEWTON_F_TOL;
      break;
    }
    cblas_dscal(opt->f_size, -1.0, f, 1);
  
    opt->df(df, t, x, u, p, data);                                                /* JACOBIAN EVALUATION */
    
    lapack_int sol_ret = -1;
    sol_ret = LAPACKE_dgels(opt->ordering, 'N', opt->f_size, opt->x_size, 1, df, lda, f, ldb);
    if (sol_ret != 0) {
      if (sol_ret > 0)
        ret = NEWTON_SINGULAR_JACOBIAN;
      else
        ret = NEWTON_ILLEGAL_JACOBIAN;
      break;
    }
    
    /* Update Step condition */
    opt->x_tol = cblas_dnrm2(opt->x_size, f, 1);
    if (opt->x_tol < set_x_tol) {
      ret = NEWTON_X_TOL;
      break;
    }
    
    cblas_daxpy(opt->x_size, 1, f, 1, x, 1);
    counts++;
  }
  
  ret = opt->max_iter <= counts ? NEWTON_MAX_ITER : ret;
  opt->max_iter = counts;  

  free(f);
  free(df);
  return ret;
}
