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

#include <cblas.h>
#include "libeuler.h"

/**
 * @brief Internal: Strut that will be passed to euler functions callbacks
 */
typedef struct euler_passthrough
{
  double ts;              /**< Integration steps. Taken from options struct */
  double alpha;           /**< Tustin coefficient. Taken from options struct */
  lapack_int u_offset;    /**< Input offset for \f$t+h\f$ callbacks. Taken from options struct */
  euler_ode_function f;   /**< Vector field to integrate. Taken from input struct */
  euler_ode_jacobian df;  /**< Jacobian of the vector field. Taken from options struct */
  const double *xk;       /**< Current state. Taken from parameters */
  lapack_int x_size;      /**< Ode dimension, taken from the input struct */
  double *work_f;         /**< Working space. Allocated in integration step */
  double *work_df;        /**< Working space. Allocated in integration step */
  void *data;             /**< User supplied data. Taken from parameters */
} euler_passtrough;

/**
 * @brief Euler implicit step wrapper
 * 
 * The vector field wrapper implements the following:
 * \f{
 *   g(\cdot) = -x(t+h) + x(t) + (1-\alpha) h f(x(t), u_{1..u_{off}, p} + 
 *            \alpha h f(x(t+h), u_{u_{off}..dim(u)}, p)
 * \f}
 * by using user supplied vector field callback. The implicit step search the root
 * of \f$g(\cdot)\f$.
 */
void euler_function_wrapper(double *f, const double t, const double *x, const double *u, const double **p, void *data);

/**
 * @brief Euler implicit step wrapper jacobian
 * 
 * The vector field wrapper implements the following jacobian matrix:
 * \f{
 *   \nabla_{x(t+h)} g(\cdot) = -I + \alpha h \nabla_{x(t+h)} f(x(t+h), u_{u_{off}..dim(u)}, p) 
 * \f}
 * by using user supplied vector field callback
 */
void euler_jacobian_wrapper(double *f, const double t, const double *x, const double *u, const double **p, void *data);

    euler_ret euler(const euler_options *opt, double *xp, const double t, const double *x, const double *u, const double **p, void *data)
{

  /* EXPLICIT IMPLEMENTATION */
/* Performin a very simple step if alpha == 0 */
if (opt->alpha == 0)
{
  opt->f(xp, t, x, u, p, opt->data);
  cblas_dscal(opt->x_size, opt->ts, xp, 1);
  cblas_daxpy(opt->x_size, 1, x, 1, xp, 1);
  return EULER_SUCCESS;
  }

  /* IMPLICIT IMPLEMENTTION */
  /* Allocating working memory */
  double *work_f = (double *)calloc(2 * opt->x_size, sizeof(double));
  if (!work_f)
    return EULER_EMALLOC;
  double *work_df = (double *)calloc(2 * opt->x_size * opt->x_size, sizeof(double));
  if(!work_df) {
    free(work_f);
    return EULER_EMALLOC;
  }
  /* Setting up the identity matrix in work_df second space, once forever */
  for (lapack_int i = 0; i < opt->x_size; i++)
    work_df[(opt->x_size * opt->x_size) + i + i * opt->x_size] = -1.0;

  /* Setting up options for Euler step */
  newton_options newton_opts = {
    opt->ordering, opt->x_size, opt->x_size, 
    opt->s_tol, opt->x_tol, opt->max_iter,
    euler_function_wrapper,
    euler_jacobian_wrapper
  };

  euler_passtrough pt = {
    opt->ts, opt->alpha, opt->u_offset,
    opt->f, opt->df, x, opt->x_size,
    work_f, work_df,
    opt->data
  };

  cblas_dcopy(opt->x_size, x, 1, xp, 1);
  newton_ret nwt = newton_solve(&newton_opts, t, xp, u, p, ((void *)&pt)); 

  /* Freeing space */
  free(work_f);
  free(work_df);
  if (nwt > NEWTON_MAX_ITER)
    return EULER_GENERIC;
  return EULER_SUCCESS;
}

void euler_function_wrapper(double *f, const double t, const double *x, const double *u, const double **p, void *data) {
  euler_passtrough *_data = ((euler_passtrough *)data);

  /* Evaluating f(x(k), u(k)) and f(x(k+1), u(k+1)), storing result in work_f */
  _data->f(_data->work_f, t, _data->xk, u, p, _data->data);
  _data->f(_data->work_f + _data->x_size, t, x, u + _data->u_offset, p, _data->data);

  /* Computing: x(k) - x(k+1) + (1-alpha) ts f(x(k), u(k)) + alpha ts f(x(k+1), u(k+1)) */
  cblas_dscal(_data->x_size, (1 - _data->alpha) * _data->ts, _data->work_f, 1);
  cblas_dscal(_data->x_size, _data->alpha * _data->ts, _data->work_f + _data->x_size, 1);
  cblas_daxpy(_data->x_size, 1, _data->work_f + _data->x_size, 1, _data->work_f, 1);
  cblas_daxpy(_data->x_size, -1, x, 1, _data->work_f, 1);
  cblas_daxpy(_data->x_size, 1, _data->xk, 1, _data->work_f, 1);

  /* Copying result in output */
  cblas_dcopy(_data->x_size, _data->work_f, 1, f, 1);
}

void euler_jacobian_wrapper(double *df, const double t, const double *x, const double *u, const double **p, void *data) {
  euler_passtrough *_data = ((euler_passtrough *)data);

  /* Evaluating JAC(f)(x(k+1), u(k+1)) */
  _data->df(_data->work_df, t, x, u + _data->u_offset, p, _data->data);

  /* Computing: -I + alpha ts JAC(f)(x(k+1), u(k+1)) (still using level 1 Blas) */
  cblas_dscal(_data->x_size * _data->x_size, _data->alpha * _data->ts, _data->work_df, 1);
  cblas_daxpy(_data->x_size * _data->x_size, 1, _data->work_df + _data->x_size * _data->x_size, 1, _data->work_df, 1);

  /* Copying result in output */
  cblas_dcopy(_data->x_size * _data->x_size, _data->work_df, 1, df, 1);
}