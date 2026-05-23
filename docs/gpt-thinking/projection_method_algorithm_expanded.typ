#set page(margin: 22mm)
#set par(justify: true)

= Projection Method Algorithm For The Cavity Solver

This note describes the projection-method solver used by the current
finite-volume lid-driven cavity code. It is split into three layers:

1. the mathematical theory behind projection methods;
2. a solver-independent algorithm description;
3. an implementation-oriented description that maps each step to the current
   code path.

== Why A Time Step Appears In Projection Methods

Projection methods are often introduced first as transient solvers for
incompressible Navier-Stokes equations. In that setting the unknown flow field is
not only a spatial field, but a sequence of states in time:

$
(bold(u)^n, p^n) -> (bold(u)^(n+1), p^(n+1))
$

A time step is therefore needed to state how far the solution advances between
two neighboring physical states. The momentum equation contains the physical
acceleration term $rho partial_t bold(u)$. A first-order time discretization
turns it into the inertial term:

$
(rho (bold(u)^(n+1) - bold(u)^n)) / (Delta t)
$

This term is what lets a transient solver remember the previous velocity and
advance the solution by one physical time increment $Delta t$. The projection
idea then splits this time advance into two easier subproblems:

1. predict an intermediate velocity from the momentum equation;
2. correct that velocity with a pressure increment so the new velocity satisfies
   incompressibility.

For a true transient simulation, $Delta t$ is a physical time step. Reducing it
changes the temporal resolution and the physical trajectory being computed.

The current code uses the same algebraic projection structure for a steady
lid-driven cavity solve. Here `projection_dt` should be read as a pseudo-time
step: it controls the numerical march toward a steady state, but the saved
result is not a physical snapshot at a particular time. The solver stops when
the pseudo-time updates, momentum residuals, and continuity residual are small
enough.


== Theoretical Foundation Of Projection Methods

The projection method is best understood as a way to enforce the incompressible
constraint after a momentum update. In a compressible flow, pressure can be tied
to density and temperature through an equation of state. In a constant-density
incompressible flow, the density is already fixed, so pressure no longer plays
that thermodynamic role. Instead, pressure is the Lagrange multiplier that keeps
the velocity field inside the constraint set

$
bold(u) in cal(V), quad nabla dot bold(u) = 0
$

with the appropriate wall-normal boundary condition. The projection step takes a
velocity field that has been advanced by convection, diffusion, and forcing, and
then removes exactly the irrotational part that violates the divergence-free
constraint.

=== Pressure As A Constraint Force

For constant density, the incompressible equations can be written as

$
partial_t bold(u)
= - bold(u) dot nabla bold(u) + nu nabla^2 bold(u) - (1 / rho) nabla p
$

$
nabla dot bold(u) = 0
$

where $nu = mu / rho$. The pressure gradient is not an optional force added to a
compressible equation. It is the force required to prevent the momentum update
from leaving the divergence-free velocity space. This is why no independent
transport equation for pressure appears in the incompressible formulation.

A useful way to see this is to apply divergence to the momentum equation. Since
$partial_t (nabla dot bold(u)) = 0$, pressure must satisfy a Poisson-type
constraint equation:

$
nabla^2 p
= rho nabla dot ( - bold(u) dot nabla bold(u) + nu nabla^2 bold(u) )
$

This equation says: once the non-pressure terms try to create compression or
expansion, pressure adjusts globally so that the final velocity remains
incompressible. Projection methods discretize this idea in a split form.

=== Helmholtz-Hodge Decomposition

The mathematical basis of projection methods is the Helmholtz-Hodge
decomposition. Under suitable boundary conditions, a vector field can be split
into a divergence-free component and a pure gradient component:

$
bold(w) = bold(w)_perp + nabla q
$

$
nabla dot bold(w)_perp = 0
$

The projection operator $cal(P)$ maps a general vector field to its
divergence-free component:

$
cal(P) bold(w) = bold(w) - nabla q
$

where $q$ is chosen so that

$
nabla dot (bold(w) - nabla q) = 0
$

which gives

$
nabla^2 q = nabla dot bold(w)
$

In this language, the projection method computes an intermediate velocity
$bold(u)^*$ and then applies a discrete version of $cal(P)$:

$
bold(u)^(n+1) = cal(P) bold(u)^*
$

The pressure-increment Poisson equation is therefore not an arbitrary numerical
trick. It is the scalar equation needed to construct the gradient field that
must be subtracted from $bold(u)^*$.

=== Projection As A Constrained Minimization

The same correction can also be derived from a constrained minimization problem.
Among all velocities satisfying incompressibility, choose the corrected velocity
closest to the predicted velocity in the kinetic-energy norm:

$
min_(bold(u))
(1 / 2) integral_Omega (rho / Delta t) abs(bold(u) - bold(u)^*)^2 d Omega
$

subject to

$
nabla dot bold(u) = 0
$

Introduce a Lagrange multiplier $phi$ for the divergence constraint. The
stationarity condition gives

$
(rho / Delta t) (bold(u) - bold(u)^*) + nabla phi = 0
$

so

$
bold(u) = bold(u)^* - (Delta t / rho) nabla phi
$

Enforcing $nabla dot bold(u) = 0$ gives the familiar pressure-increment Poisson
equation:

$
nabla^2 phi = (rho / Delta t) nabla dot bold(u)^*
$

This interpretation is important: the projection step changes the predicted
velocity as little as possible while restoring the incompressibility constraint.
The part removed from $bold(u)^*$ is a gradient field, and gradient fields are
orthogonal to divergence-free fields in the continuous $L^2$ inner product when
the boundary terms are compatible.

=== The Role Of The Pressure Increment

There are two common ways to formulate a projection method.

1. *Non-incremental projection.* The predictor ignores pressure or uses an old
   approximation only weakly, and the Poisson solve directly produces the new
   pressure $p^(n+1)$.
2. *Incremental projection.* The predictor uses the old pressure $p^n$, and the
   Poisson solve produces a pressure correction
   $phi = p^(n+1) - p^n$.

The current solver uses the second form. The predictor already includes the
known pressure force $-nabla p^k$, so the correction only needs to supply the
missing increment $-nabla phi$. This usually gives a better pressure update than
the original non-incremental Chorin projection, especially for flows where
pressure changes gradually between iterations.

The pressure increment has one important nullspace property: adding a constant
to $phi$ does not change $nabla phi$, and therefore does not change the corrected
velocity. A reference pressure or pressure-increment cell must be fixed to make
the algebraic system nonsingular.

=== Boundary Conditions In The Projection Step

Boundary conditions are the most delicate part of projection methods. Taking the
normal component of the velocity correction gives

$
bold(u)^(n+1) dot bold(n)
= bold(u)^* dot bold(n) - (Delta t / rho) partial_n phi
$

If the corrected normal velocity is prescribed as $g_n$, then the compatible
Neumann condition for the pressure increment is

$
partial_n phi = (rho / Delta t) (bold(u)^* dot bold(n) - g_n)
$

For an impermeable wall, $g_n = 0$. If the predictor already satisfies
$bold(u)^* dot bold(n) = 0$ at the wall, this reduces to the homogeneous
condition

$
partial_n phi = 0
$

This is the condition used by the cavity solver for the pressure increment. The
moving lid imposes tangential velocity, but the wall-normal velocity is still
zero, so the normal pressure-increment condition remains homogeneous after the
predictor has been made boundary-consistent.

Tangential velocity boundary conditions are more subtle. The projection
correction naturally controls the normal divergence constraint; it does not, by
itself, guarantee that the tangential no-slip value remains exact unless the
velocity boundary conditions are reapplied after correction. This is why practical
finite-volume codes usually enforce cavity ghost cells both before and after the
projection step.

=== Solvability Of The Poisson Equation

For a pure Neumann Poisson problem, the right-hand side must satisfy a
compatibility condition. In continuous form,

$
integral_Omega nabla^2 phi d Omega
= integral_(partial Omega) partial_n phi d Gamma
$

Using the projection equation,

$
integral_Omega (rho / Delta t) nabla dot bold(u)^* d Omega
= integral_(partial Omega) (rho / Delta t) bold(u)^* dot bold(n) d Gamma
$

so the volume integral of the predicted divergence must match the boundary flux.
For a closed cavity with no normal flow, the total predicted mass imbalance
should sum to zero. In finite precision and on collocated grids, the assembled
right-hand side can contain small inconsistencies. Pinning one pressure degree
of freedom removes the constant nullspace, but it does not by itself fix a large
incompatibility in the mass balance; the face-flux construction must be
consistent with the pressure equation.

=== Discrete Projection On A Finite-Volume Grid

At the discrete level, the projection should mimic the continuous identity

$
bold(u)^(n+1) = bold(u)^* - (Delta t / rho) nabla phi
$

and the constraint

$
nabla dot bold(u)^(n+1) = 0
$

A compact matrix form is

$
bold(U)^(n+1) = bold(U)^* - M^(-1) G Phi
$

$
D bold(U)^(n+1) = 0
$

where:

- $bold(U)$ is the discrete velocity vector or face-flux vector;
- $M$ is the velocity mass or response operator;
- $G$ is the discrete pressure-gradient operator;
- $D$ is the discrete divergence operator;
- $Phi$ is the pressure-increment vector.

Substitution gives the discrete pressure equation

$
D M^(-1) G Phi = D bold(U)^*
$

For the present first-order pseudo-time projection, $M$ is essentially the
scalar factor $rho / Delta t$ times the control-volume measure, so
$M^(-1)$ corresponds to the response factor $Delta t / rho$. In a more strongly
implicit formulation, $M$ may contain the full momentum matrix or an
approximation to its diagonal.

A good finite-volume projection depends on the compatibility of $D$ and $G$.
If the divergence operator used in the Poisson equation is not the same one used
in the continuity residual, the code can solve the pressure equation accurately
but still report a nonzero mass residual. On collocated grids, pressure and
velocity live at the same cell centers, so naive interpolation can also create
checkerboard pressure modes. Rhie-Chow-type face reconstruction is used to build
pressure-sensitive face velocities and suppress this decoupling.

=== Flux Projection Versus Cell-Centered Projection

The continuous formula corrects the velocity field. A finite-volume continuity
equation, however, is written in terms of face fluxes. Therefore the most
important practical object is the corrected face-normal velocity or mass flux:

$
F_f^(n+1) = F_f^* - rho A_f d_f (phi_N - phi_P) / delta_f
$

where $A_f$ is the face area, $delta_f$ is the distance between neighboring cell
centers, and $d_f$ is the velocity response at the face. The pressure-increment
coefficients are built from exactly these face responses. When the same face
correction is used for assembly and for the final residual check, each pressure
control volume satisfies a discrete mass balance.

The cell-centered updates

$
u_P^(n+1) = u_P^* - d_P partial_x phi
$

and

$
v_P^(n+1) = v_P^* - d_P partial_y phi
$

are still needed for the next momentum solve, but the finite-volume continuity
constraint is most directly enforced on the faces.

=== Splitting Error And Projection Variants

Projection methods are fractional-step methods, so they introduce a splitting
error because the momentum and pressure-constraint parts are not solved as one
fully coupled system. The basic first-order structure is

$
"momentum predictor" -> "Poisson projection"
$

The error is usually acceptable for a steady pseudo-transient solve when the
iteration is driven to convergence, but it matters for accurate transient
simulations. Three variants are common:

#figure(
  table(
    columns: 4,
    align: left,
    inset: 6pt,
    stroke: none,
    table.hline(),
    [*Variant*], [*Predictor pressure*], [*Poisson unknown*], [*Main feature*],
    table.hline(stroke: 0.6pt),
    [Non-incremental],
    [None or weak old pressure],
    [$p^(n+1)$],
    [Simple, but pressure boundary accuracy is relatively poor.],

    [Incremental],
    [$p^n$],
    [$phi = p^(n+1) - p^n$],
    [Better pressure update; this is the form used by the current solver.],

    [Rotational incremental],
    [$p^n$],
    [$phi$ plus divergence/viscous correction],
    [Reduces pressure splitting error near boundaries for transient flows.],
    table.hline(),
  ),
  caption: [Common projection-method variants.],
)

The current implementation should therefore be interpreted as an incremental,
pseudo-transient pressure-correction method. It is close in spirit to SIMPLE
because both methods predict velocity, solve a pressure-correction equation, and
correct velocity and pressure. The difference is that SIMPLE often uses momentum
diagonal coefficients as under-relaxed velocity responses, whereas the current
projection method uses the pseudo-time response $Delta t / rho$.

=== Meaning Of Convergence For A Pseudo-Transient Projection

When the method is used for a steady cavity problem, the pseudo-time derivative
is only a continuation device. During iteration, the predictor equation contains

$
(rho / Delta t) (bold(u)^* - bold(u)^k)
$

which stabilizes the nonlinear march. At convergence,

$
bold(u)^(k+1) approx bold(u)^k
$

and the pseudo-time contribution becomes small. The remaining equations are then
the steady momentum balance and the discrete continuity equation.

The choice of $Delta t$ controls the numerical path to the steady solution:

- a smaller $Delta t$ gives a stronger diagonal pseudo-time term, usually making
  each linear solve more stable but the global march slower;
- a larger $Delta t$ weakens the pseudo-time damping and can accelerate
  convergence, but it can also make the nonlinear iteration less robust;
- for a steady solve, changing $Delta t$ should not change the final converged
  solution if the same discrete equations are converged tightly enough.


== Deriving The Equations To Solve

Start from the constant-density incompressible Navier-Stokes equations:

$
nabla dot bold(u) = 0
$

$
rho partial_t bold(u) + C(bold(u), bold(u)) - D(bold(u)) = -nabla p
$

Here $C$ represents convection and $D$ represents viscous diffusion. A direct
implicit step from time level $n$ to $n+1$ would ask for both
$bold(u)^(n+1)$ and $p^(n+1)$ at the same time:

$
(rho (bold(u)^(n+1) - bold(u)^n)) / (Delta t)
+ C(bold(u)^n, bold(u)^(n+1))
- D(bold(u)^(n+1))
= -nabla p^(n+1)
$

$
nabla dot bold(u)^(n+1) = 0
$

This is the coupled pressure-velocity problem. The pressure is not determined
by an independent equation of state; instead, it is whatever scalar field makes
the new velocity divergence-free.

Projection methods introduce a pressure increment:

$
phi = p^(n+1) - p^n
$

and split the coupled step. First, replace the unknown new pressure by the known
old pressure and solve for an intermediate velocity:

$
(rho (bold(u)^* - bold(u)^n)) / (Delta t)
+ C(bold(u)^n, bold(u)^*)
- D(bold(u)^*)
= -nabla p^n
$

This is the momentum-prediction equation. It gives a velocity that has felt
convection, diffusion, boundary forcing, and the old pressure gradient, but it
does not yet have the pressure increment needed to enforce continuity.

Next, subtract the predictor equation from the coupled target equation. The
projection approximation places the convection and diffusion changes in the
predictor step, leaving the correction step to carry the pressure increment:

$
(rho (bold(u)^(n+1) - bold(u)^*)) / (Delta t) = -nabla phi
$

Solving this expression for the corrected velocity gives:

$
bold(u)^(n+1) = bold(u)^* - (Delta t) / (rho) nabla phi
$

Now impose the incompressibility constraint on the corrected velocity:

$
0 = nabla dot bold(u)^(n+1)
= nabla dot bold(u)^* - (Delta t) / (rho) nabla^2 phi
$

Therefore the pressure increment must satisfy:

$
nabla^2 phi = (rho) / (Delta t) nabla dot bold(u)^*
$

This is the pressure-increment Poisson equation. In the code, `phi` is stored as
`pressure_correction`. Once it is solved, the pressure update is:

$
p^(n+1) = p^n + phi
$

On a finite-volume grid, the divergence equation is written as a face-flux
balance. For one pressure control volume $P$, the predicted mass imbalance is:

$
m_P^* =
rho Delta y (u_w^* - u_e^*)
+ rho Delta x (v_s^* - v_n^*)
$

The pressure-increment equation is then assembled in stencil form:

$
a_P phi_P
- a_E phi_E
- a_W phi_W
- a_N phi_N
- a_S phi_S
= m_P^*
$

For the projection solver, the velocity response factors are:

$
d_u = d_v = (Delta t) / (rho)
$

so the neighbor coefficients have the usual face-response form, for example:

$
a_E = (rho Delta y d_e) / (Delta x)
$

This is the discrete equation actually solved by the pressure-increment step.
The first pressure-increment cell is pinned because adding a constant to
$phi$ does not change its gradient.

== Part 1: Projection Method As An Algorithm

The projection method solves incompressible Navier-Stokes flow by separating a
velocity prediction step from a pressure projection step. The idea is to first
advance the velocity with the momentum equation, producing an intermediate
velocity that is not necessarily divergence-free, and then project that velocity
onto a divergence-free space by solving a pressure-increment equation.

For a constant-density incompressible flow,

$
nabla dot bold(u) = 0
$

$
rho (partial_t bold(u) + bold(u) dot nabla bold(u))
= -nabla p + mu nabla^2 bold(u)
$

The implementation here uses the method as a pseudo-transient iteration toward a
steady solution. The pseudo-time step is written as $Delta t$, and the known
state at iteration $k$ is $(bold(u)^k, p^k)$.

=== Step 1: Momentum Prediction

First compute an intermediate velocity $bold(u)^*$ from the momentum equation
using the current pressure $p^k$:

$
(rho (bold(u)^* - bold(u)^k)) / (Delta t)
+ C(bold(u)^k, bold(u)^*)
- D(bold(u)^*) = -nabla p^k
$

Here:

- $C$ is the discretized convection operator.
- $D$ is the discretized viscous diffusion operator.
- The old velocity $bold(u)^k$ appears in the pseudo-time source term.
- The predicted velocity $bold(u)^*$ does not yet satisfy exact continuity.

For the two velocity components this corresponds to solving one linear system
for $u^*$ and one for $v^*$.

=== Step 2: Pressure-Increment Poisson Equation

The projection correction has the form:

$
bold(u)^(k+1) = bold(u)^* - (Delta t) / (rho) nabla phi
$

where $phi$ is the pressure increment. Enforcing incompressibility gives:

$
nabla dot bold(u)^(k+1) = 0
$

Substitute the correction formula:

$
nabla dot bold(u)^* - (Delta t) / (rho) nabla^2 phi = 0
$

so the pressure increment satisfies:

$
nabla^2 phi = (rho) / (Delta t) nabla dot bold(u)^*
$

In finite-volume form, this is assembled from the predicted face-flux imbalance
of each pressure control volume. One pressure-increment degree of freedom must
be pinned because only pressure gradients affect the velocity.

=== Step 3: Velocity Projection

After solving for $phi$, update the cell-centered velocity:

$
u^(k+1) = u^* - (Delta t) / (rho) partial_x phi
$

$
v^(k+1) = v^* - (Delta t) / (rho) partial_y phi
$

This is the projection step: it removes the divergence contained in the
predicted velocity field.

=== Step 4: Pressure Update

The pressure is updated incrementally:

$
p^(k+1) = p^k + phi
$

The absolute pressure level is arbitrary for incompressible flow, so the
corrected pressure is shifted by subtracting one reference cell value after the
update.

=== Step 5: Boundary Conditions And Convergence

Velocity uses Dirichlet wall conditions for the lid-driven cavity. Pressure and
pressure increment use zero-normal-gradient wall conditions:

$
partial_n p = 0, quad partial_n phi = 0
$

The iteration continues until:

- the momentum linear solves reach their tolerance;
- the corrected face-flux continuity residual is small;
- the maximum velocity update between two pseudo-time states is small;
- the minimum iteration count has been reached.

== Part 2: Mapping The Algorithm To This Codebase

The current implementation lives primarily in `ProjectionSolver` and reuses the
same structured-grid and finite-volume helpers as the SIMPLE solver.

=== State Variables

The main state is stored in `FlowFields`:

- `u`, `v`: corrected cell-centered velocity, corresponding to
  $u^k, v^k$ or $u^(k+1), v^(k+1)$ depending on iteration stage.
- `u_star`, `v_star`: predicted velocity, corresponding to $u^*, v^*$.
- `pressure`: corrected pressure, corresponding to $p^k$.
- `pressure_correction`: pressure increment, corresponding to $phi$.
- `d_u`, `d_v`: velocity response factors. For projection they are set to
  `projection_dt / density`, matching $(Delta t) / (rho)$.
- `u_face`, `v_face`: Rhie-Chow reconstructed face velocities used for
  convection, continuity residuals, and pressure-increment assembly.

=== Initialization

`ProjectionSolver::ProjectionSolver()`:

1. builds the `StructuredGrid`;
2. allocates `FlowFields`;
3. checks that `projection_dt` is positive;
4. calls `set_projection_correction_factors()`;
5. applies cavity boundary conditions;
6. reconstructs initial face velocities.

This prepares a boundary-consistent zero-flow initial state.

=== Per-Iteration Code Flow

The high-level loop is `ProjectionSolver::run()`.

#figure(
  table(
    columns: (24%, 34%, 42%),
    align: (left, left, left),
    inset: 6pt,
    stroke: none,
    table.hline(),
    [*Algorithm stage*], [*Code path*], [*Role*],
    table.hline(stroke: 0.6pt),
    [Prepare state],
    [#text(size: 7.5pt)[`prepare_iteration()`]],
    [Set `d_u`, `d_v`; apply cavity ghost cells; refresh corrected face velocities.],

    [Predict $u^*$],
    [#text(size: 7.5pt)[`solve_u_predictor()`#linebreak()
      `-> assemble_u_projection_momentum()`]],
    [Solve the pseudo-transient x-momentum predictor with current pressure and old `u`.],

    [Stage face fluxes],
    [#text(size: 7.5pt)[`refresh_face_velocities(` #linebreak()
      `fields_.u_star, fields_.v)`]],
    [Let the v-predictor see the updated x-face fluxes.],

    [Predict $v^*$],
    [#text(size: 7.5pt)[`solve_v_predictor()`#linebreak()
      `-> assemble_v_projection_momentum()`]],
    [Solve the pseudo-transient y-momentum predictor with current pressure and old `v`.],

    [Solve pressure increment],
    [#text(size: 7.5pt)[`solve_pressure_increment_step()`#linebreak()
      `-> assemble_pressure_correction()`]],
    [Build the Poisson-like pressure-increment equation from predicted face imbalance.],

    [Project velocity],
    [#text(size: 7.5pt)[`project_pressure_and_velocity()`]],
    [Apply `u = u_star - d_u grad(phi)` and `v = v_star - d_v grad(phi)`.],

    [Refresh diagnostics],
    [#text(size: 7.5pt)[`project_fields_and_`#linebreak()
      `collect_metrics()`]],
    [Apply boundaries, rebuild corrected faces, compute residual metrics.],
    table.hline(),
  ),
  caption: [Projection algorithm stages and their corresponding implementation paths.],
)

=== Momentum Predictor Assembly

`assemble_u_projection_momentum()` and `assemble_v_projection_momentum()` reuse
the same upwind-convection and central-diffusion stencil as the SIMPLE
momentum assembly. The projection-specific change is the pseudo-transient term:

$
a_P^("proj") = a_P^("base") + (rho Delta x Delta y) / (Delta t)
$

The right-hand side includes the current pressure force and the old-velocity
source:

$
b_P^("u,proj") = ((p_W - p_E) Delta y) / (2)
+ (rho Delta x Delta y) / (Delta t) u_P^k
+ S_U^("u")
$

$
b_P^("v,proj") = ((p_S - p_N) Delta x) / (2)
+ (rho Delta x Delta y) / (Delta t) v_P^k
+ S_U^("v")
$

The moving-lid source contribution remains in the u-equation north boundary
treatment through $S_U^("u") = 2 a_N U_"lid"$.

=== Pressure-Increment Assembly

`assemble_pressure_correction()` is shared with SIMPLE, but in projection it is
called after:

$
d_u = d_v = (Delta t) / (rho)
$

The function reconstructs predictor face velocities from `u_star`, `v_star`,
`pressure`, `d_u`, and `d_v`, then forms the pressure-control-volume mass
imbalance:

$
m_P^* =
rho Delta y (u_w^* - u_e^*)
+ rho Delta x (v_s^* - v_n^*)
$

The pressure-increment system uses this imbalance as the right-hand side:

$
a_P phi_P
- a_E phi_E
- a_W phi_W
- a_N phi_N
- a_S phi_S
= m_P^*
$

The first pressure cell is pinned with:

$
phi_(1,1) = 0
$

Immediately after the pressure increment is loaded, the solver reapplies cavity
boundary conditions so the `pressure_correction` ghost cells satisfy the
zero-normal-gradient condition before any velocity correction gradient is read.

=== Projection Correction In Code

`project_pressure_and_velocity()` computes central pressure-increment gradients:

$
g_("pc,x") = (phi_(i+1,j) - phi_(i-1,j)) / (2 Delta x)
$

$
g_("pc,y") = (phi_(i,j+1) - phi_(i,j-1)) / (2 Delta y)
$

Then it updates:

$
u_(i,j) = u^*_(i,j) - d_u_(i,j) g_("pc,x")
$

$
v_(i,j) = v^*_(i,j) - d_v_(i,j) g_("pc,y")
$

$
p_(i,j) = p_(i,j) + phi_(i,j)
$

Finally, the pressure reference level is removed by subtracting
`pressure(1, 1)` from every interior pressure cell.

=== Verification Path

The implementation is guarded by:

- `test_projection_coefficients`: checks the pseudo-transient momentum
  coefficients and pressure-increment coefficient scaling.
- `test_projection_smoke`: checks convergence, finite fields, cavity wall
  conditions, and pressure-increment ghost consistency on an 8x8 case.
- `test_cli_smoke`: checks default projection dispatch, explicit SIMPLE
  dispatch, explicit projection dispatch, output files, and invalid solver
  handling.
- `test_projection_validation_run` plus `test_projection_validation`: runs a
  32x32 Re=100 projection case and compares centerline results with the Ghia
  benchmark through `scripts/validate_cavity.py`.
