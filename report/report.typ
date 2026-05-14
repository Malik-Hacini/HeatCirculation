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


= Direct problem

== Variational formulation

Let $Gamma_t$ be the top wall, $Gamma_b$ the bottom door, and $Gamma_l$,
$Gamma_r$ the lateral walls. We fix the source intensities
$alpha = (alpha_1, dots, alpha_6)$ and set $f = sum_(i=1)^6 alpha_i 1_(C_i)$.
Since the $C_i$ have finite measure, $f in L^2(Omega)$. The strong problem is
$
cases(
  -op("div")(kappa nabla T) = f & "in" Omega,
  T = T_c & "on" Gamma_t,
  partial_n T = 0 & "on" Gamma_l union Gamma_r,
  partial_n T + h(T - T_e) = 0 & "on" Gamma_b,
)
$

Following the course treatment of mixed boundary value problems, we encode the
Dirichlet condition in the function space. Let
$V = { v in H^1(Omega) | gamma_0 v = 0 " on " Gamma_t }$ and
$K = { w in H^1(Omega) | gamma_0 w = T_c " on " Gamma_t }$.

Assume first that $T$ is regular. For $v in V$, multiply the equation by $v$ and
integrate over $Omega$:
$
integral_Omega f v dif x = integral_Omega -op("div")(kappa nabla T) v dif x.
$
Green's formula gives
$
integral_Omega -op("div")(kappa nabla T) v dif x
  = integral_Omega kappa nabla T dot nabla v dif x
    - integral_(partial Omega) kappa partial_n T v dif s.
$
We now split the boundary term. On $Gamma_t$, $v = 0$ by definition of $V$. On
$Gamma_l union Gamma_r$, the walls are insulated, so $partial_n T = 0$. On
$Gamma_b$, the Robin condition gives $partial_n T = -h(T - T_e)$, hence
$-kappa partial_n T v = kappa h T v - kappa h T_e v$. Therefore
$
integral_Omega kappa nabla T dot nabla v dif x
  + integral_(Gamma_b) kappa h T v dif s
  = integral_Omega f v dif x + integral_(Gamma_b) kappa h T_e v dif s.
$
This identity only involves first derivatives of $T$ and is the weak formulation.
For $u, v in H^1(Omega)$, define
$
a(u, v) = integral_Omega kappa nabla u dot nabla v dif x
  + integral_(Gamma_b) kappa h u v dif s.
$
The weak problem is thus: find $T in K$ such that, for all $v in V$,
$
a(T, v) = integral_Omega f v dif x + integral_(Gamma_b) kappa h T_e v dif s.
$

We cannot directly apply Lax-Milgram on $K$ because $K$ is affine, not a vector
space.
Choose any lifting $G in H^1(Omega)$ such that $gamma_0 G = T_c$ on $Gamma_t$;
for instance the constant function $G = T_c$. We write $T = U + G$, where
$U in V$. Then $T$ solves the weak problem if and only if $U in V$ satisfies
$a(U, v) = ell(v)$ for all $v in V$, with
$
ell(v) = integral_Omega f v dif x + integral_(Gamma_b) kappa h T_e v dif s - a(G, v).
$
The term $a(G, v)$ is simply the contribution of the lifting of the nonhomogeneous
Dirichlet condition. It appears when replacing $T$ by $U + G$ in $a(T, v)$.

The space $V$ is a closed subspace of $H^1(Omega)$, hence a Hilbert space for the
$H^1$ norm.
We now check the three assumptions of the Lax-Milgram theorem.

First, $a$ is continuous on $V times V$. Indeed, $kappa in L^infinity(Omega)$ and the
trace theorem gives $norm(gamma_0 v)_(L^2(Gamma_b)) <= C norm(v)_(H^1(Omega))$.
Thus there is $M > 0$ such that
$abs(a(u, v)) <= M norm(u)_(H^1(Omega)) norm(v)_(H^1(Omega))$.

Second, $a$ is coercive on $V$. Since $kappa >= 1$,
$a(v, v) >= norm(nabla v)_(L^2(Omega))^2$. By the Poincare inequality on the
space of functions with zero trace on $Gamma_t$, there exists $C_P > 0$ such that
$norm(v)_(H^1(Omega)) <= C_P norm(nabla v)_(L^2(Omega))$. Hence
$a(v, v) >= C_P^(-2) norm(v)_(H^1(Omega))^2$.

Third, $ell$ is continuous on $V$. The volume term is controlled by Cauchy-Schwarz,
the boundary term by the trace theorem, and the lifting term by the continuity of
$a$.

By Lax-Milgram, there is a unique $U in V$, hence a unique weak solution
$T = U + G in K$.

== P1 finite element approximation

Let $cal(T)_h$ be a triangulation of $Omega$. The $P_1$ finite element space is
$
V_h = { v_h in C^0(overline(Omega)) | v_h|_K " is affine for all " K in cal(T)_h,
  v_h = 0 " on " Gamma_t }.
$
With $G_h = T_c$, we look for $T_h = U_h + G_h$ with $U_h in V_h$ such that
$a(U_h, v_h) = ell(v_h)$ for all $v_h in V_h$.

Let $(phi_1, dots, phi_N)$ be the usual nodal basis of $V_h$. Writing
$U_h = sum_(j=1)^N U_j phi_j$ and testing with each $phi_i$ gives the linear
system $A_h U = F_h$, with $(A_h)_(i j) = a(phi_j, phi_i)$ and
$(F_h)_i = ell(phi_i)$. The matrix is sparse because each hat function has local
support. By applying Lax-Milgram on the finite-dimensional Hilbert space $V_h$,
this discrete problem has a unique solution.

#pagebreak()
== Stationary simulation

The script `direct.edp` implements the $P_1$ method for the stationary direct
problem. For the test run, we used the prescribed six source locations and the
known intensities $alpha_i = 250000$ for all $i$. The direct simulation gives
$min T_h = 343 K$, $max T_h = 507.137 K$, and $overline(T_h)_S = 422.015 K$.

#figure(
  image("figures/direct.png", width: 80%),
  caption: [Stationary temperature distribution for the direct problem, with the cooking region $S$ outlined. Warmer colors indicate higher temperature.],
)

The figure shows hot spots around the resistors and a non-uniform temperature in
$S$, which motivates the optimization problem of Part 2.

== Warm-up phase

To describe the warm-up phase, we add the time derivative and consider
$partial_t T - op("div")(kappa nabla T) = f$ in $Omega times ]0, T_f[$, with the
same boundary conditions as before. We use the following implicit Euler scheme:
given $T_h^n$, find $T_h^(n+1) in G_h + V_h$ such that, for all
$v_h in V_h$,
$
integral_Omega (T_h^(n+1) v_h) / delta t dif x + a(T_h^(n+1), v_h)
  = integral_Omega (T_h^n v_h) / delta t dif x
    + integral_Omega f v_h dif x
    + integral_(Gamma_b) kappa h T_e v_h dif s.
$
The script `warmup.edp` implements this scheme with $T_h^0 = T_e$, using $200$
time steps on the interval $[0, 1.25]$. With the same source intensities as in
`direct.edp`, the mean temperature in $S$ goes from $295.491 K$ at
$t = 0.00625$ to $421.986 K$ at $t = 1$, and reaches $422.012 K$ at $t = 1.25$.
This is consistent with the convergence of the parabolic problem toward the
stationary direct solution.

#pagebreak()
= Optimization problem

== Linear decomposition of the solution

For a coefficient vector $alpha = (alpha_1, dots, alpha_6)$, set
$f_alpha = sum_(i=1)^6 alpha_i 1_(C_i)$. The direct strong problem is obtained
from Part 1 by replacing $f$ with $f_alpha$.

Define $T_0$ as the solution with no internal source,
$
cases(
  -op("div")(kappa nabla T_0) = 0 & "in" Omega,
  T_0 = T_c & "on" Gamma_t,
  partial_n T_0 = 0 & "on" Gamma_l union Gamma_r,
  partial_n T_0 + h(T_0 - T_e) = 0 & "on" Gamma_b,
)
$
and, for $i = 1, dots, 6$, define $T_i$ as the unit response of source $C_i$,
$
cases(
  -op("div")(kappa nabla T_i) = 1_(C_i) & "in" Omega,
  T_i = 0 & "on" Gamma_t,
  partial_n T_i = 0 & "on" Gamma_l union Gamma_r,
  partial_n T_i + h T_i = 0 & "on" Gamma_b,
)
$

$T_0$ is the temperature when all resistors are switched off, and $T_i$ is the
temperature response generated by source $i$ alone with homogeneous boundary
data. Let
$ T^*(alpha) = T_0 + sum_(i=1)^6 alpha_i T_i. $
By linearity of the differential operator and of the boundary operators,
$T^*(alpha)$ satisfies the same strong problem as $T(alpha)$. By uniqueness of
the solution, we conclude that
$
T(alpha, x) = T_0 (x) + sum_(i=1)^6 alpha_i T_i (x).
$

== Linear formulation of the optimization problem
Set $e_0(x) = T_0(x) - T_s$. We have
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
$A alpha^* = b$,
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
It shows that $A$ is positive definite and therefore invertible. We conclude that
the linear system, and thus the optimization problem, has
a unique solution.

== Simulation results

The script `initial_optimized.edp` implements the previous construction with
$P_1$ finite elements. It computes $T_0$,
then the six functions $T_i$, assembles the discrete normal system
$A^h alpha^h = b^h$, solves it, and finally performs a direct simulation with
the resulting source term.

#figure(
  image("figures/opti_initial.png", width: 80%),
  caption: [Optimized temperature distribution.],
)

For $T_s = 400 K$, the direct simulation with the optimized coefficients gives
$
overline(T)_S &= 399.901 K quad "and "
sqrt(1 / abs(S) integral_S abs(T - T_s)^2 dif x) &= 2.49 K.
$
The mean temperature in $S$ is therefore very close to the target $400 K$, and
the normalized error in $S$ is about $2.49 K$.

The result is symmetric under the reflection $x -> -x$, as the domain, the
coefficient $kappa$, the region $S$, the source layout, and the boundary
conditions are invariant under this reflection. There is no analogous symmetry
in the vertical direction: the top wall is fixed at $T_c = 343 K$, whereas the
bottom wall exchanges heat with an exterior temperature $T_e = 293 K$, so the
bottom sources remain cooler than the top ones.


== Source layout optimization

We then tested different source layouts while keeping the same radius
$r_C = 0.05$ and the same target temperature $T_s = 400 K$. The best tested
configuration uses eight sources: five on the upper row and three on the lower
row,
$
(-1, 0.75), (-0.5, 0.75), (0, 0.75), (0.5, 0.75), (1, 0.75)
$
and
$
(-1, -0.75), (0, -0.75), (1, -0.75).
$
We chose this distribution because the original optimization already required
larger coefficients on the upper sources. Adding more sources on the upper side
therefore gives finer control where it is most useful. We also keep symmetry
with respect to $x = 0$, since the oven geometry, the cooking region, and the
target temperature are symmetric in the horizontal direction.

#figure(
  image("figures/opti_source.png", width: 80%),
  caption: [Optimized temperature distribution for the proposed eight-source layout.],
)

The script `5t_3b_optimized.edp` gives
$
overline(T)_S &= 399.994 K quad "and "
sqrt(1 / abs(S) integral_S abs(T - T_s)^2 dif x) &= 0.804 K.
$
For comparison, the initial six-source layout gave a normalized error of about
$2.49 K$. The proposed layout therefore reduces the error by about $68$ percent,
while keeping all optimized source intensities positive and thus physically
feasible. We also observe the expected smoother temperature distribution in $S$.

#pagebreak()
= Inverse problem

== Detection method

We assume that the normal operating state of the oven is known. Thus, from a
measured temperature $T_m$, we first subtract the no-anomaly background $T_0$.
The anomaly $U = T_m - T_0$ satisfies the same homogeneous boundary conditions
as the responses $T_i$ from the optimization part, but with one unknown source
$alpha 1_(C(z))$, where $z$ is the unknown center.

For a candidate center $z$, let $u_z$ be the solution with unit intensity:
$
a(u_z, v) = integral_(C(z)) v dif x quad "for all" v in V.
$
If $B_M$ denotes the measurement operator, the predicted measurement is
$alpha B_M u_z$. We use two choices: $B_M w = w$ on $Gamma_l union Gamma_r$,
which corresponds to measuring the temperature on the vertical walls, and
$B_M w = partial_n w$ on $Gamma_t$, which corresponds to measuring the normal
derivative on the top wall.

For each candidate $z$, the best intensity is obtained explicitly by a
one-dimensional least-squares fit:
$
alpha(z) = ((d_M, B_M u_z)_M) / ((B_M u_z, B_M u_z)_M),
$
where $d_M = B_M (T_m - T_0)$ and $(dot, dot)_M$ is the $L^2$ product on the
measured boundary. The remaining relative residual is
$
r(z) = sqrt(((d_M - alpha(z) B_M u_z, d_M - alpha(z) B_M u_z)_M) / ((d_M, d_M)_M)).
$
The numerical inverse problem is then reduced to a finite search over candidate
centers: choose the $z$ minimizing $r(z)$, and take $alpha(z)$ as the recovered
intensity.

#pagebreak()
== Numerical results

The script `inverse.edp` implements this dictionary search. We generated a
synthetic anomaly with center $z^* = (0.4, -0.3)$, radius $0.05$, and intensity
$alpha^* = 5 times 10^5$. The candidate grid is
$29 times 19$ points in $[-1.4, 1.4] times [-0.9, 0.9]$; this grid contains the
true center.

#figure(
  image("figures/inverse.png", width: 80%),
  caption: [Temperature distribution with unknown source and recovered center drawn in blue.],
)
Using the temperature on the left and right walls, the algorithm returns
$z = (0.4, -0.3)$, $alpha = 5 times 10^5$, with relative residual $0$. Using
$partial_n T$ on the top wall, it also returns $z = (0.4, -0.3)$ and
$alpha = 5 times 10^5$, with relative residual $1.39 times 10^(-8)$. The small
non-zero residual in the second case comes from the finite element evaluation of
the boundary normal derivative.

These results confirm that, in the noiseless synthetic setting, both proposed
measurement types are sufficient to detect the anomaly when the candidate grid
contains the true source. In practice, the inverse problem remains sensitive to
noise and to the grid resolution because the heat equation smooths the source
before it is observed at the boundary. If the true source is not on the grid,
the same method identifies the nearest candidate response; refining the grid
then improves the location estimate.

#pagebreak()

#heading(level: 1, numbering: none)[Appendix] <app-linear-independence>
#heading(level: 2, numbering: none)[Linear independence of the temperature responses]

We prove #footnote[Please note that this proof was not AI-generated. One of the authors has previously studied complex and harmonic analysis and was happy to apply it in this context.] that $T_1|_S, dots, T_6|_S$ are linearly independent in $L^2(S)$ for
the project geometry. The source disks are outside $S$, their closures are
pairwise disjoint, and none of them touches $partial S$. The coefficient
$kappa$ is constant equal to $1$ in $S$ and constant equal to $10$ outside $S$.

Assume that $sum_(i=1)^6 beta_i T_i = 0$ in $S$, and set
$W = sum_(i=1)^6 beta_i T_i$. Then
$-op("div")(kappa nabla W) = sum_(i=1)^6 beta_i 1_(C_i)$ in $Omega$, with
homogeneous boundary conditions. Since the sources are away from $partial S$,
$W$ is harmonic on both sides of each open side segment of $partial S$.

On the interior side of $partial S$, $W = 0$, hence both its trace and normal
derivative vanish there. The standard transmission conditions for a weak solution of
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
term in the equation, so the interior normal derivative also vanishes. Inside $C_i$,
we have $-10 Delta W = beta_i$. Integrating over $C_i$ gives
$
beta_i abs(C_i)
  = -10 integral_(C_i) Delta W dif x
  = -10 integral_(partial C_i) partial_n W dif s
  = 0.
$
Since $abs(C_i) > 0$, we get $beta_i = 0$. Thus the restrictions $T_i|_S$ are
linearly independent in $L^2(S)$.

Note that the circular shape is not essential. This proof works for any finite
family of pairwise disjoint Lipschitz open source regions of positive measure,
compactly contained in $Omega minus overline(S)$, provided
$Omega minus (overline(S) union union_i overline(C_i))$ is connected.


// #bibliography("references.bib", style: "ieee", title: [References])
