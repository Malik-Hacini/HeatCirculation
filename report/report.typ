#import "@preview/hei-synd-report:0.1.1": *
#import "metadata.typ": *
#import "extra.typ": *
//#show:make-glossary
//#register-glossary(entry-list)

//-------------------------------------
// Template config
//
#show: report.with(
  option: option,
  doc: doc,
  date: date,
  tableof: tableof,
)

//-------------------------------------
// Content
//
#counter(page).update(1)

#let project-heading-numbering(..nums) = {
  let ns = nums.pos()
  if ns.len() == 1 {
    numbering("1.", ns.at(0))
  } else if ns.len() == 2 {
    numbering("a)", ns.at(1))
  } else {
    numbering("a.1", ns.at(1), ns.at(2))
  }
}

#set heading(numbering: project-heading-numbering)

#heading(level: 1, numbering: none)[Introduction]

= Direct problem 
= Optimization problem

We use the same notation as in the statement. The oven is
$Omega = ] -1.5, 1.5 [ times ] -1, 1 [$, the cooking region is
$S = ] -0.75, 0.75 [ times ] -0.5, 0.5 [$, and the six sources are the disks
$C_i$ of radius $0.05$ centered at the prescribed points. We write
$Gamma_t$ for the top wall, $Gamma_b$ for the bottom door, and $Gamma_l$,
$Gamma_r$ for the left and right walls. The coefficient is
$kappa = 1$ in $S$ and $kappa = 10$ in $Omega minus S$.
The constants are $T_c = 343 K$ and $T_e = 293 K$: $T_c$ is the imposed
temperature on the top wall, while $T_e$ is the exterior air temperature in the
Robin law on the bottom door.

The top boundary condition is Dirichlet, so the natural affine space for the
temperature is
$ K = { w in H^1(Omega) | w = T_c " on " Gamma_t } $
and the test space is
$ V = { v in H^1(Omega) | v = 0 " on " Gamma_t } . $
For $u in H^1(Omega)$ and $v in V$, define
$
a(u, v)
  = integral_Omega kappa nabla u dot nabla v dif x
    + integral_(Gamma_b) kappa h u v dif s.
$
The bottom condition is the one written in the worksheet,
$partial_n T + h(T - T_e) = 0$. Therefore, after applying Green's formula to
$-op("div")(kappa nabla T)$, the Robin contribution is
$integral_(Gamma_b) kappa h T v dif s$ on the left-hand side and
$integral_(Gamma_b) kappa h T_e v dif s$ on the right-hand side. The lateral
insulation conditions are homogeneous Neumann conditions and therefore produce
no boundary term.

== Linear decomposition of the solution

For a coefficient vector $alpha = (alpha_1, dots, alpha_6)$, the weak direct
problem is: find $T(alpha) in K$ such that, for all $v in V$,
$
a(T(alpha), v)
  = integral_Omega (sum_(i=1)^6 alpha_i 1_(C_i)) v dif x
    + integral_(Gamma_b) kappa h T_e v dif s.
$
First, define $T_i$ for $i = 0,...6$ as the solutions of the following weak problems:
$
a(T_0, v) = integral_(Gamma_b) kappa h T_e v dif s
quad "and "
a(T_i, v) = integral_(C_i) v dif x
quad ,"for all" v in V.
$

$T_0$ is the temperature when all resistors are switched
off, depending only on the effect of the external temperature $T_e$ on the boundary, and $T_i$ is the temperature when only the resistor $i$ is on.
#pagebreak()
Let:
$ T^*(alpha) = T_0 + sum_(i=1)^6 alpha_i T_i. 
$
Since $T_0 in K$ and each $T_i in V$, we have $T^*(alpha) in K$. The bilinearity of $a$ ensures that:
$
a(T^*(alpha), v)
  &= integral_(Gamma_b) kappa h T_e v dif s
     + integral_Omega (sum_(i=1)^6 alpha_i 1_(C_i)) v dif x.
$
Thus $T^*(alpha)$ satisfies exactly the variational formulation of the direct
problem. By uniqueness of the solution $T(alpha, x)$ of the direct problem established in [LABEL PART 1],
we conclude that
$
T(alpha, x) = T_0 (x) + sum_(i=1)^6 alpha_i T_i (x).
$

== Optimization problem
Set $e_0(x) = T_0(x) - T_s$.  We have
$
J(alpha) = 1 / 2 integral_S abs(T(alpha, x) - T_s)^2 dif x = 1 / 2 integral_S abs(e_0 + sum_(j=1)^6 alpha_j T_j)^2 dif x.
$
For each
$i = 1, dots, 6$, differentiating under the integral gives
$
partial_(alpha_i) J(alpha)
  = integral_S (e_0 + sum_(j=1)^6 alpha_j T_j) T_i dif x.
$
Setting the gradient to 0 yields
$
sum_(j=1)^6 (integral_S T_i T_j dif x) alpha_j^*
  = integral_S (T_s - T_0) T_i dif x,
quad i = 1, dots, 6.
$
Therefore the optimization problem is reduced to the linear system
$A alpha^* = b,
$
where
$
A_(i j) = integral_S T_i T_j dif x,
quad
b_i = integral_S (T_s - T_0) T_i dif x.
$

The matrix $A$ is symmetric and positive semidefinite because, for every
$beta in RR^6$,
$beta^T A beta
  = integral_S abs(sum_(i=1)^6 beta_i T_i)^2 dif x >= 0.
$
To conclude that $A$ is positive definite, we need the linear independence of
the restrictions $T_1|_S, dots, T_6|_S$ in $L^2(S)$. In the present geometry
this can be proved by using the transmission conditions and unique continuation
for harmonic functions. Since this argument is outside the main scope of the
course, we give it in the #link(<app-linear-independence>)[appendix].
It shows that $A$ is definite thus invertible. We can conclude that the linear system, thus the optimization problem, has
a unique solution.

== Simulation results

The script `Project/optimization.edp` implements the previous construction with
$P_1$ finite elements. It computes $T_0$,
then the six functions $T_i$, assembles the discrete normal system
$A^h alpha^h = b^h$, solves it, and finally performs a direct simulation with
the resulting source term.

#figure(
  image("figures/opti_initial.png", width: 100%),
  caption: [Optimized temperature distribution, with the cooking region $S$ outlined. Warmer colors indicate higher temperature.],
)

For $T_s = 400 K$, the direct simulations with the optimized coefficients give
$
overline(T)_S &= 399.9 K quad "and "
sqrt(1 / abs(S) integral_S abs(T - overline(T)_S)^2 dif x) &= 2.49 K.
$
The mean temperature in $S$ is therefore very close to the target $400 K$, and
the $L^2$ error in $S$ is about $2.49 K$. 

The result is symmetric under the reflection $x -> -x$, as the domain, the
coefficient $kappa$, the region $S$, the source layout, and the boundary
conditions are invariant under this reflection. There is no analogous symmetry
in the vertical direction: the top wall is fixed at $T_c = 343 K$, whereas the
bottom wall exchanges heat with an exterior temperature $T_e = 293 K$, so the
bottom sources remain cooler than the top ones.


== Source optimization
= Inverse problem 

#heading(level: 1, numbering: none)[Conclusion]

#heading(level: 1, numbering: none)[Appendix] <app-linear-independence>
== Linear independence of the temperature responses

We prove #footnote[Please note that this proof was not AI-generated. One of the authors has studied complex and harmonic analysis before and was more than happy to apply it in this context :)]  that $T_1|_S, dots, T_6|_S$ are linearly independent in $L^2(S)$ for
the project geometry. The source disks are outside $S$, their closures are
pairwise disjoint, and none of them touches $partial S$. The coefficient
$kappa$ is constant equal to $1$ in $S$ and constant equal to $10$ outside $S$.

Assume that $sum_(i=1)^6 beta_i T_i = 0$ in $S$, and set
$W = sum_(i=1)^6 beta_i T_i$. Then
$-op("div")(kappa nabla W) = sum_(i=1)^6 beta_i 1_(C_i)$ in $Omega$, with
homogeneous boundary conditions. Since the sources are away from $partial S$,
$W$ is harmonic on both sides of each open side segment of $partial S$.

On the side of $S$, $W = 0$, hence both its trace and normal derivative vanish
there. The standard transmission conditions for a weak solution of
$-op("div")(kappa nabla W) = F$ give continuity of $W$ and of the normal flux
$kappa partial_n W$ across $partial S$. Thus, on the exterior side of any open
side segment of $partial S$, we also have $W = 0$ and $partial_n W = 0$.

Near such a segment, extend the exterior harmonic function by zero across the
segment. The zero Cauchy data make this extension weakly harmonic; by elliptic
regularity it is harmonic. It vanishes on a non-empty open set, so by analyticity
of harmonic functions it vanishes in a whole neighbourhood outside $S$. The
exterior region obtained by removing $overline(S)$ and the closed disks
$overline(C_i)$ is connected, and $W$ is harmonic there. Analytic continuation
therefore gives $W = 0$ throughout this exterior region.

Now fix one disk $C_i$. Since $W = 0$ outside the disk, the trace of $W$ and the
exterior normal derivative vanish on $partial C_i$. There is no boundary measure
in the equation, so the interior normal derivative also vanishes. Inside $C_i$,
we have $-10 Delta W = beta_i$. Integrating over $C_i$ gives
$
beta_i abs(C_i)
  = -10 integral_(C_i) Delta W dif x
  = -10 integral_(partial C_i) partial_n W dif s
  = 0.
$
Since $abs(C_i) > 0$, we get $beta_i = 0$. Thus the restrictions $T_i|_S$ are linearly
independent in $L^2(S)$.

Please note that the circular shape is not essential. This proofs works for any finite family
of pairwise disjoint Lipschitz open source regions of positive measure,
compactly contained in $Omega minus overline(S)$, provided
$Omega minus (overline(S) union union_i overline(C_i))$ is connected.


#bibliography("references.bib", style: "ieee", title: [References])
