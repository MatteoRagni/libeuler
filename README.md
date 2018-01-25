# libeuler - Explicit or Implicit Euler integration step

The library implements an integration step for an Euler integrator, with Tustin coefficient.

If the Tustin coefficient is null, the next step is evaluated as:

<img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;x(t&plus;h)&space;=&space;x(t)&space;&plus;&space;h&space;f(t,&space;x(t),&space;u_{1..dim(u)},&space;p)" title="x(t+h) = x(t) + h f(t, x(t), u_{1..dim(u)}, p)" />

if the Tustin is different than zero, the next step is the solution of the following system
of equation:

<img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;\begin{align*}&space;g(\cdot)&space;=&space;&&space;x(t)&space;-&space;x(t&plus;h)&space;&plus;&space;\\&space;&&space;(1-\alpha)\,&space;h\,&space;f(t,&space;x(t),&space;u_{1..q},&space;p)&space;&plus;&space;\\&space;&&space;\alpha\,&space;h\,&space;f(t&space;&plus;&space;h,&space;x(t&space;&plus;&space;h),&space;u_{q..dim(u)},&space;p)&space;\end{align*}" title="\begin{align*} g(\cdot) = & x(t) - x(t+h) + \\ & (1-\alpha)\, h\, f(t, x(t), u_{1..q}, p) + \\ & \alpha\, h\, f(t + h, x(t + h), u_{q..dim(u)}, p) \end{align*}" />

The solution is found by using the Newton algorithm, a multistep method in which
the solution is iteratively found as solution of:

<img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;x^&plus;(t&plus;h)&space;=&space;x(t&plus;h)&space;&plus;&space;\nabla^{-1}g(\cdot)\,g(\cdot)" title="x^+(t+h) = x(t+h) + \nabla^{-1}g(\cdot)\,g(\cdot)" />

where the Jacobian is defined as:

<img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;\nabla&space;g(\cdot)&space;=&space;I&space;&plus;&space;\alpha&space;\nabla&space;f(t&plus;h,&space;x(t&plus;h),&space;u_{q..dim(u)},&space;p)" title="\nabla g(\cdot) = I + \alpha \nabla f(t+h, x(t+h), u_{q..dim(u)}, p)" />

## Dependencies

 * Lapack (in particular Lapacke)
 * Cblas

The two dependencies are available in Linux as packages in all major distribution. The two
dependencies are also distributed with MATLAB. To see how to compile a model with this integrator
check the Makefile

## Usage Example

Let's make an usage example and a comparison with the output of the equivalent Simulink model. 
The ODE to be solved is:

<img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;\begin{align*}&space;\dot{x}_1(t)&space;=&&space;a_1&space;(a_2&space;u(t)&space;-&space;a_3&space;\sqrt{\gamma&space;x_1(t)})&space;\\&space;\dot{x}_2(t)&space;=&&space;b_1&space;(a_3&space;\sqrt{\gamma&space;x_1(t)}&space;-&space;b_3&space;\sqrt{\gamma&space;x_2(t)})&space;\\&space;\end{align*}" title="\begin{align*} \dot{x}_1(t) =& a_1 (a_2 u(t) - a_3 \sqrt{\gamma x_1(t)}) \\ \dot{x}_2(t) =& b_1 (a_3 \sqrt{\gamma x_1(t)} - b_3 \sqrt{\gamma x_2(t)}) \\ \end{align*}" />

and it has the following Jacobian:

<img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;\nabla&space;f(\cdot)&space;=&space;\left[&space;\begin{align*}&space;-2&space;a_1&space;a_3&space;\sqrt{\gamma&space;x_1(t)^{-1}}&space;&&space;&&space;0&space;\\&space;b_1&space;a_3&space;\sqrt{\gamma&space;x_1(t)^{-1}}&space;&&space;&&space;b_1&space;b_3&space;\sqrt{\gamma&space;x_2(t)^{-1}}&space;\end{align*}&space;\right&space;]" title="\nabla f(\cdot) = \left[ \begin{align*} -2 a_1 a_3 \sqrt{\gamma x_1(t)^{-1}} & & 0 \\ b_1 a_3 \sqrt{\gamma x_1(t)^{-1}} & & b_1 b_3 \sqrt{\gamma x_2(t)^{-1}} \end{align*} \right ]" />

 * the implicit implementation is in `test/test_euleri.c`
 * the explicit implementation is in `test/test_eulere.c`

Those two immages summarizes the results:

![Result of dynamic](.fig/dyn.png) 

![Result of dynamic, integration error](.fig/err.png) 

