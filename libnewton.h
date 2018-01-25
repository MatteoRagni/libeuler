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

#ifndef LIBNEWTON_H_
#define LIBNEWTON_H_

#include <lapacke.h>
#include <cblas.h>

/**
 * @brief Jacobian Callback for the Newton algorithm
 * 
 * The callback receivs the current input, the current control, and the list of
 * parameter vectors. It should store the first argument the Jacobian matrix (stored as an array).
 * To select the ordering in the matrix, you must select the correct ordering in newton_options 
 * structure, at the ordering input.
 * @param df output vector (vectorized Jacobian matrix)
 * @param t current time for evaluation
 * @param x current point for evaluation
 * @param u current control for evaluation
 * @param p array of parameter vectors
 * @param data user space input (simply use it as a pointer casted to void)
 */
typedef void (*newton_jacobian)(
    double *df,
    const double t,
    const double *x,
    const double *u,
    const double **p,
    void *data);

/**
 * @brief Vector Field Callback for the Newton algorithm
 * 
 * The callback receivs the current input, the current control, and the list of
 * parameter vectors. It should store the first argument the function output vector.
 * @param df output vector (vectorized Jacobian matrix)
 * @param t current time for evaluation
 * @param x current point for evaluation
 * @param u current control for evaluation
 * @param p array of parameter vectors
 * @param data user space input (simply use it as a pointer casted to void)
 */
typedef void (*newton_function)(
    double *f,
    const double t,
    const double *x,
    const double *u,
    const double **p,
    void *data);

/**
 * @brief Options for the Newton Algorithm
 * 
 * This structure contains all the options for the Newton algorithm, alongside the callbacks.
 * The structure will be modified by the algorithm with some debug information, such as number
 * of iterations, tolerances and ordering for the jacobian matrix.
 */
typedef struct newton_options {
  lapack_int ordering; /**< Should be LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR */
  lapack_int f_size;   /**< Vector field size */
  lapack_int x_size;   /**< Variable vector size */
  double f_tol;        /**< Stopping tolerance for the zero. 
                              At the end will contain the 2 norm of the 
                              vector field in the solution */
  double x_tol;        /**< Stopping tolerance for the x vector step.
                              At the end will contain the 2 norm of the last update step */
  lapack_int max_iter; /**< Maximum number of iteration,
                              At the end will contain the number of step executed  */
  newton_function f;   /**< Pointer to vector field callback */
  newton_jacobian df;  /**<  Pointer to Jacobian callback */
} newton_options;

/**
 * @brief Error code returned by the algorithm
 */
typedef enum newton_ret {
  NEWTON_F_TOL = 0,         /**< (0) The solver reached the required tolerance limit for the vector field */
  NEWTON_X_TOL,             /**< (1) The last step for solution update was less than the minimum */
  NEWTON_MAX_ITER,          /**< (2) Maximum number of iterations reached */
  NEWTON_SINGULAR_JACOBIAN, /**< (3 LAPACKE) The jacobian is singular */
  NEWTON_ILLEGAL_JACOBIAN,  /**< (4 LAPACKE) Illegal jacobian */
  NEWTON_MALLOC_ERROR,      /**< (5) Cannot allocate memory */
  NEWTON_GENERIC_ERROR      /**< (6) Generic error in the execution of the algorithm */
} newton_ret;

/**
 * @brief Boolean implementation
 */
typedef enum newton_bool {
  NEWTON_FALSE = 0, /**< False */
  NEWTON_TRUE       /**< True */
} newton_bool;

/**
 * @brief Executes the Newton algorithm for root finding
 * 
 * Executes the Newton algorithm for root finding. The step Performed is:
 * \f{
 *   x_{k+1} - x_{k} = -\nabla F^{-1}(x_k, u, p) F(x_k, u, p)
 * \f}
 * and the solution is found by using DGELS defined in LAPACK library.
 * The stopping conditions are:
 *  * Number of iterations bigger than maximum allowed (specified in newton_options)
 *  * \f$ |x_{k+1} - x_k| \leq x_{tol}\f$
 *  * \f$ |f(x_k, u, p)| \leq y_{tol}\f$
 * and the function and jacobian are evaluated through callbacks. On exit, the values 
 * inside the option structure are update for debuggin purposes (this is why Newton 
 * options is not a const pointer). It returns a status enum.
 * @param opt option structure
 * @param x root position and initial condition. Will be modified
 * @param u control action input. It can be NULL.
 * @param p parameter array of vectors. It can be NULL.
 * @param data space for user data. Will be passed to callbacks. It can be NULL
 * @return a status exit code as described in newton_ret enum.
 */
newton_ret newton_solve(
    newton_options *opt,
    const double t,
    double *x,
    const double *u,
    const double **p,
    void *data);

#endif
