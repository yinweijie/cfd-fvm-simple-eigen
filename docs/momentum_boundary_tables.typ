#set page(margin: 22mm)
#set par(justify: true)

= Discretization Coefficient Expressions

This note summarizes the coefficient expressions used by the main
discretization operators in the current SIMPLE solver.

When a control volume touches a cavity wall, the missing neighbor is not assembled
as an off-diagonal matrix entry. Instead, the corresponding coefficient is folded
back into the diagonal coefficient $a_P$. For the north wall in the $u$ equation,
that diagonal fold-back is paired with the moving-lid source term.

The tables below list the *positive* link-coefficient magnitudes. In the sparse
matrix, assembled neighbor entries are written as $-a_W$, $-a_E$, $-a_S$, and $-a_N$.

== Notation

The notation follows standard finite-volume conventions:

- $P$ is the current control volume center.
- $W$, $E$, $S$, and $N$ are the west, east, south, and north neighboring cells.
- $w$, $e$, $s$, and $n$ denote the corresponding cell faces.
- $a_W$, $a_E$, $a_S$, and $a_N$ are stored below as positive link-coefficient magnitudes.

== Boundary Handling Through Ghost Cells

All cell-centered fields are stored with one ghost-cell layer, but only interior
cells are unknowns in the assembled linear systems. Boundary conditions enter the
algebra by substituting ghost-cell values into the stencil.

For stationary velocity walls,

$
u_G = -u_P, quad v_G = -v_P
$

so a boundary-side contribution of the form $a_B phi_G$ becomes

$
a_B phi_G = a_B (-phi_P) = -a_B phi_P
$

which is moved to the left-hand side and therefore adds $+a_B$ to the diagonal.
This is why stationary walls do not create an explicit constant source term in the
momentum equations in this implementation.

For the moving lid in the `u` equation,

$
u_G = 2 U_"lid" - u_P
$

and the north-side contribution becomes

$
a_N u_G = a_N (2 U_"lid" - u_P) = -a_N u_P + 2 a_N U_"lid"
$

The $-a_N u_P$ part is absorbed into the diagonal, while the remaining constant
part becomes the source contribution $2 a_N U_"lid"$.

Pressure and pressure-correction use homogeneous Neumann wall conditions,

$
p_G = p_P, quad p'_G = p'_P
$

so they are handled differently from the velocity Dirichlet walls.
Their wall treatment does not inject a moving-wall-type source term.

== Definitions

$
a_E = D_x + max(-F_e, 0), quad
a_W = D_x + max(F_w, 0), quad
a_N = D_y + max(-F_n, 0), quad
a_S = D_y + max(F_s, 0)
$

$
D_x = mu (Delta y) / (Delta x), quad
D_y = mu (Delta x) / (Delta y)
$

$
F_e = rho u_e Delta y, quad
F_w = rho u_w Delta y, quad
F_n = rho v_n Delta x, quad
F_s = rho v_s Delta x
$

$
a_P^("base") = a_E + a_W + a_N + a_S + (F_e - F_w + F_n - F_s)
$

== Momentum Right-Hand Sides

The assembled momentum equations can be written in the standard linearized form

$
a_P u_P = a_W u_W + a_E u_E + a_S u_S + a_N u_N + b_P^("u")
$

$
a_P v_P = a_W v_W + a_E v_E + a_S v_S + a_N v_N + b_P^("v")
$

with complete right-hand sides

$
b_P^("u") = (p_W - p_E) Delta y / 2 + ((1 - alpha_u) / alpha_u) a_P^("base") u_P + S_U^("u")
$

$
b_P^("v") = (p_S - p_N) Delta x / 2 + ((1 - alpha_v) / alpha_v) a_P^("base") v_P + S_U^("v")
$

Here $S_U^("u")$ and $S_U^("v")$ are the boundary source terms introduced by
the $S_P / S_U$ linearization given below.

== U-Momentum

#figure(
  table(
    columns: 7,
    align: center,
    stroke: none,
    inset: 8pt,
    table.hline(),
    [*Region*], [*$a_W^("link")$*], [*$a_E^("link")$*], [*$a_S^("link")$*], [*$a_N^("link")$*], [*$a_P$*], [*$S_U$*],
    table.hline(stroke: 0.6pt),

    [Interior],      [$a_W$], [$a_E$], [$a_S$], [$a_N$], [$a_P^("base") / alpha_u$], [$0$],
    [West boundary], [$0$],   [$a_E$], [$a_S$], [$a_N$], [$a_P^("base") / alpha_u + a_W$], [$0$],
    [East boundary], [$a_W$], [$0$],   [$a_S$], [$a_N$], [$a_P^("base") / alpha_u + a_E$], [$0$],
    [South boundary],[$a_W$], [$a_E$], [$0$],   [$a_N$], [$a_P^("base") / alpha_u + a_S$], [$0$],
    [North boundary],[$a_W$], [$a_E$], [$a_S$], [$0$],   [$a_P^("base") / alpha_u + a_N$], [$2 a_N U_"lid"$],
    [South-West],    [$0$],   [$a_E$], [$0$],   [$a_N$], [$a_P^("base") / alpha_u + a_W + a_S$], [$0$],
    [South-East],    [$a_W$], [$0$],   [$0$],   [$a_N$], [$a_P^("base") / alpha_u + a_E + a_S$], [$0$],
    [North-West],    [$0$],   [$a_E$], [$a_S$], [$0$],   [$a_P^("base") / alpha_u + a_W + a_N$], [$2 a_N U_"lid"$],
    [North-East],    [$a_W$], [$0$],   [$a_S$], [$0$],   [$a_P^("base") / alpha_u + a_E + a_N$], [$2 a_N U_"lid"$],

    table.hline(),
  ),
  caption: [*u*-momentum equation link coefficients and boundary source terms.],
)

== V-Momentum

#figure(
  table(
    columns: 7,
    align: center,
    stroke: none,
    inset: 8pt,
    table.hline(),
    [*Region*], [*$a_W^("link")$*], [*$a_E^("link")$*], [*$a_S^("link")$*], [*$a_N^("link")$*], [*$a_P$*], [*$S_U$*],
    table.hline(stroke: 0.6pt),

    [Interior],      [$a_W$], [$a_E$], [$a_S$], [$a_N$], [$a_P^("base") / alpha_v$], [$0$],
    [West boundary], [$0$],   [$a_E$], [$a_S$], [$a_N$], [$a_P^("base") / alpha_v + a_W$], [$0$],
    [East boundary], [$a_W$], [$0$],   [$a_S$], [$a_N$], [$a_P^("base") / alpha_v + a_E$], [$0$],
    [South boundary],[$a_W$], [$a_E$], [$0$],   [$a_N$], [$a_P^("base") / alpha_v + a_S$], [$0$],
    [North boundary],[$a_W$], [$a_E$], [$a_S$], [$0$],   [$a_P^("base") / alpha_v + a_N$], [$0$],
    [South-West],    [$0$],   [$a_E$], [$0$],   [$a_N$], [$a_P^("base") / alpha_v + a_W + a_S$], [$0$],
    [South-East],    [$a_W$], [$0$],   [$0$],   [$a_N$], [$a_P^("base") / alpha_v + a_E + a_S$], [$0$],
    [North-West],    [$0$],   [$a_E$], [$a_S$], [$0$],   [$a_P^("base") / alpha_v + a_W + a_N$], [$0$],
    [North-East],    [$a_W$], [$0$],   [$a_S$], [$0$],   [$a_P^("base") / alpha_v + a_E + a_N$], [$0$],

    table.hline(),
  ),
  caption: [*v*-momentum equation link coefficients and boundary source terms.],
)

== Equivalent Boundary Linearization in $S_P / S_U$ Form

The same wall treatment can be written as

$
a_P = a_P^("base") / alpha - S_P
$

where $S_P <= 0$ folds the missing boundary-side link back into the diagonal and
$S_U$ carries any explicit boundary source term.

#figure(
  table(
    columns: 3,
    align: center,
    stroke: none,
    inset: 8pt,
    table.hline(),
    [*Region*], [*$S_P^("u")$*], [*$S_U^("u")$*],
    table.hline(stroke: 0.6pt),

    [Interior],      [$0$], [$0$],
    [West boundary], [$-a_W$], [$0$],
    [East boundary], [$-a_E$], [$0$],
    [South boundary],[$-a_S$], [$0$],
    [North boundary],[$-a_N$], [$2 a_N U_"lid"$],
    [South-West],    [$-(a_W + a_S)$], [$0$],
    [South-East],    [$-(a_E + a_S)$], [$0$],
    [North-West],    [$-(a_W + a_N)$], [$2 a_N U_"lid"$],
    [North-East],    [$-(a_E + a_N)$], [$2 a_N U_"lid"$],

    table.hline(),
  ),
  caption: [Equivalent boundary linearization for the *u*-momentum equation.],
)

#figure(
  table(
    columns: 3,
    align: center,
    stroke: none,
    inset: 8pt,
    table.hline(),
    [*Region*], [*$S_P^("v")$*], [*$S_U^("v")$*],
    table.hline(stroke: 0.6pt),

    [Interior],      [$0$], [$0$],
    [West boundary], [$-a_W$], [$0$],
    [East boundary], [$-a_E$], [$0$],
    [South boundary],[$-a_S$], [$0$],
    [North boundary],[$-a_N$], [$0$],
    [South-West],    [$-(a_W + a_S)$], [$0$],
    [South-East],    [$-(a_E + a_S)$], [$0$],
    [North-West],    [$-(a_W + a_N)$], [$0$],
    [North-East],    [$-(a_E + a_N)$], [$0$],

    table.hline(),
  ),
  caption: [Equivalent boundary linearization for the *v*-momentum equation.],
)

== Rhie-Chow Face Velocities

Interior cavity faces are reconstructed from neighboring cell-centered velocities
plus a pressure-correction term. Boundary faces on the physical walls are set to
zero normal velocity instead of using these formulas.

$ 
u_e = (u_P + u_E)/2 - d_e ((p_E - p_P) / Delta x - (g_(x,P) + g_(x,E)) / 2)
$

$ 
v_n = (v_P + v_N)/2 - d_n ((p_N - p_P) / Delta y - (g_(y,P) + g_(y,N)) / 2)
$

$
d_e = (d_(u,P) + d_(u,E)) / 2, quad
d_n = (d_(v,P) + d_(v,N)) / 2
$

The west and south face formulas are analogous, with the neighboring cell labels
changed from $E$ or $N$ to $W$ or $S$.

== Pressure-Correction Coefficients

The pressure-correction equation uses the predictor face velocities to form the
mass-imbalance right-hand side and uses neighboring $d_u$ and $d_v$ values to form the four-point stencil.

$
a_E^("pc") = rho Delta y ((d_(u,P) + d_(u,E)) / 2) / Delta x
$

$
a_W^("pc") = rho Delta y ((d_(u,P) + d_(u,W)) / 2) / Delta x
$

$
a_N^("pc") = rho Delta x ((d_(v,P) + d_(v,N)) / 2) / Delta y
$

$
a_S^("pc") = rho Delta x ((d_(v,P) + d_(v,S)) / 2) / Delta y
$

$
a_P^("pc") = a_E^("pc") + a_W^("pc") + a_N^("pc") + a_S^("pc")
$

$
m_P = rho Delta y (u_w^* - u_e^*) + rho Delta x (v_s^* - v_n^*)
$

The code replaces the south-west corner row $(1, 1)$ with a pressure reference:
$a_P = 1$ and $m_P = 0$.

#figure(
  table(
    columns: 7,
    align: center,
    stroke: none,
    inset: 8pt,
    table.hline(),
    [*Region*], [*$a_W^("link")$*], [*$a_E^("link")$*], [*$a_S^("link")$*], [*$a_N^("link")$*], [*$a_P$*], [*RHS*],
    table.hline(stroke: 0.6pt),

    [Interior],                [$a_W^("pc")$], [$a_E^("pc")$], [$a_S^("pc")$], [$a_N^("pc")$], [$a_P^("pc")$], [$m_P$],
    [West boundary],           [$0$],          [$a_E^("pc")$], [$a_S^("pc")$], [$a_N^("pc")$], [$a_E^("pc") + a_S^("pc") + a_N^("pc")$], [$m_P$],
    [East boundary],           [$a_W^("pc")$], [$0$],          [$a_S^("pc")$], [$a_N^("pc")$], [$a_W^("pc") + a_S^("pc") + a_N^("pc")$], [$m_P$],
    [South boundary],          [$a_W^("pc")$], [$a_E^("pc")$], [$0$],          [$a_N^("pc")$], [$a_W^("pc") + a_E^("pc") + a_N^("pc")$], [$m_P$],
    [North boundary],          [$a_W^("pc")$], [$a_E^("pc")$], [$a_S^("pc")$], [$0$],          [$a_W^("pc") + a_E^("pc") + a_S^("pc")$], [$m_P$],
    [South-East],             [$a_W^("pc")$], [$0$],          [$0$],          [$a_N^("pc")$], [$a_W^("pc") + a_N^("pc")$], [$m_P$],
    [North-West],             [$0$],          [$a_E^("pc")$], [$a_S^("pc")$], [$0$],          [$a_E^("pc") + a_S^("pc")$], [$m_P$],
    [North-East],             [$a_W^("pc")$], [$0$],          [$a_S^("pc")$], [$0$],          [$a_W^("pc") + a_S^("pc")$], [$m_P$],
    [Reference cell $(1,1)$], [$0$],          [$0$],          [$0$],          [$0$],          [$1$], [$0$],

    table.hline(),
  ),
  caption: [Pressure-correction equation link coefficients and right-hand side.],
)
