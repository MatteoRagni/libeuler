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

#ifndef LIBEULER_H_
#define LIBEULER_H_

#include <cblas.h>
#include "libnewton.h" /**< Custom solver for implicit step */

/**
 * @brief Callback for the ode vector field
 * The callback for the ODE model. The callback stores the output in
 * the first pointer (f). The function does not need to allocate nor free
 * the input vector, but may overflow if exceedes the dimension that is
 * declared in the euler options structure.
 * When accessing an element of the parameter, please remember that the
 * implementation keeps in mind a MATLAB-like interface. If the parameters
 * (in MATLAB) are passed like:
 * @code
 * f(t, x, u, p1, [p2, p3])
 * @endcode
 * then, in the C code, p1, p2 and p3 may be accessed in p as:
 * @code
 * rk4_float p1, p2, p3;
 * p1 = p[0][0];
 * p2 = p[1][0];
 * p3 = p[1][1];
 * @endcode
 * If you need the time step **inside** the callback, you should directly define
 * it inside or pass it as a parameter (and remeber to update it in the euler_options
 * structure).
 * @param f otput vector for the ODE
 * @param t current time for the function
 * @param x state for the ODE evaluation
 * @param u external input for the ODE evaluation
 * @param p pointer to arrays of parameters
 * @param data auxiliary data pointer to void for user data
 * @returns nothing
 * @warning The callback may produce a buffer overflow error if tries to
 *          write more than x_size elements in xdot.
 */
typedef void (*euler_ode_function)(
    double *f,
    const double t,
    const double *x,
    const double *u,
    const double **p,
    void *data);

/**
 * @brief Callback for the ode vector field Jacobian
 * The callback for the ODE Jacobian. The callback stores the output in
 * the first pointer (df). The function does not need to allocate nor free
 * the input vector, but may overflow if exceedes the dimension that is
 * declared in the euler options structure. The matrix is stored continuosly
 * inan array, with ordering as specified in options structure.
 * @param f otput vector for the ODE
 * @param t current time for the function
 * @param x state for the ODE evaluation
 * @param u external input for the ODE evaluation
 * @param p pointer to arrays of parameters
 * @param data auxiliary data pointer to void for user data
 * @returns nothing
 * @warning The callback may produce a buffer overflow error if tries to
 *          write more than x_size elements in xdot.
 */
typedef void (*euler_ode_jacobian)(
    double *f,
    const double t,
    const double *x,
    const double *u,
    const double **p,
    void *data);

/**
 * @brief Returning value for the integrator
 */
typedef enum euler_ret {
  EULER_SUCCESS = 0,  /**< Correct execution */
  EULER_EMALLOC,      /**< Memory allocation error */
  EULER_NULLPTR,      /**< Received a null pointer */
  EULER_GENERIC       /**< Generic error raised. @todo: should be divided for Newton error */
} euler_ret;

typedef struct euler_options {
  double ts;               /**< Integration step */
  double alpha;            /**< Tustin transform coefficient. \f$\alpha \in [0,1]\f$. 
                                For \f$\alpha = 0\f$, no optimization is performed */
  lapack_int x_size;       /**< State dimensions. Will also be vector field size */
  lapack_int u_offset;     /**< Input offset for implicit steps. It can also be 0 */
  lapack_int ordering;     /**< Ordering: LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR */
  double s_tol;            /**< Tolerance for Newton solution */
  double x_tol;            /**< Step tolerance for Newton solution */
  lapack_int max_iter;     /**< Maximum number of iterations for Newton solver */
  euler_ode_function f;    /**< Actual ODE vector field */
  euler_ode_jacobian df;   /**< Jacobian of the vector field */
  void *data;              /**< Empty space for user data */
} euler_options;

/**
 * @brief Euler step (explicit or implicit)
 * 
 * The function performs the euler step.
 * If \f$ \alpha = 0\f$ the step performed is:
 * \f{
 *   x(t+h) = x(t) + f(x(t), u_{1..dim(u)}, p)
 * \f}
 * while for \f$\alpha \in (0,1]\f$ the step perfomed is the solution of
 * the following non linear problem:
 * \f{
 *   x(t+h) = x(t) + (1-\alpha) h f(x(t), u_{1..u_{off}, p} + 
 *            \alpha h f(x(t+h), u_{u_{off}..dim(u)}, p)
 * \f}
 * which is an implicit integration step. The initial guess for the solution of the 
 * non linear problem is the state \f$x(t)\f$.
 * @param opt pointer to scruct with options
 * @param xp next integration step
 * @param t current integration time
 * @param x current state
 * @param u control vector. May contain both current and next control input,
 *          in sequence: \f$u = [u(t), u(t+h)]\f$. The offset allows to get as 
 *          input in the callback the next input.
 * @param p  pointer to arrays of parameters
 * @param data void pointer to userspace data
 * @return an exit code to check if integration step succeeded
 */
euler_ret euler(
  const euler_options *opt,
  double *xp,
  const double t,
  const double *x,
  const double *u,
  const double **p,
  void *data);

#endif /* LIBEULER_H_ */