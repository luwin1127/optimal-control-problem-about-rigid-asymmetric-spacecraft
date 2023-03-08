# optimal-control-problem-about-rigid-asymmetric-spacecraft

*Direct solution of nonlinear optimal control problems using quasilinearization and Chebyshev polynomials*（DOI：10.1016/S0016-0032(02)00028-5） 

Section 5.2. Example 2: Rigid asymmetric spacecraft

The optimal control problem is given as follows.

The angular velocitires $\omega_1$, $\omega_2$ and $\omega_3$ of the spacecraft are given by

$$
\left \{ \begin{matrix}
    \dot{\omega}_1 = -\frac{I_3-I_2}{I_1} \omega_2 \omega_3 + \frac{u_1}{I_1},  \\
    \dot{\omega}_2 = -\frac{I_1-I_3}{I_2} \omega_1 \omega_3 + \frac{u_2}{I_2},  \\
    \dot{\omega}_3 = -\frac{I_2-I_1}{I_3} \omega_1 \omega_2 + \frac{u_3}{I_3},
\end{matrix}
\right.
$$

where $u_1$, $u_2$ and $u_3$ are the control torques.

The objective function is given by

$$
J = 0.5 \int_{0}^{100}(u_1^2+u_2^2+u_3^2) \text{d}t.
$$

The following initial and terminal state constraints have to be satified

$$
\begin{aligned}
    &\omega_1(0)=0.01 \ \text{r/s},\ \omega_2(0)=0.005 \ \text{r/s},\ \omega_3(0)=0.001 \ \text{r/s},    \\
    &\omega_1(t_f)=\omega_2(t_f)=\omega_3(t_f)=0 \ \text{r/s},
\end{aligned}
$$

and the requirable parameters are given by

$$
I_1=86.24 \ \text{kg m}^2,\ I_2=85.07 \ \text{kg m}^2,\ I_3=113.59 \ \text{kg m}^2.
$$

The paper considers two cases, that is,

1) case 1: no inequility constraints on the states or on the controls.
2) case 2: inequaility constraints on $\omega_1(t)$.

In case 2, $\omega_1(t)$ must satisfy the following constraint.

$$
\omega_1(t) - (5\times10^{-6}\ t^2 - 5\times10^{-4}\ t + 0.015 ) \le 0.
$$

Results are given in the matlab codes.

Run `section5_example2.m`.
