#set page(margin: 22mm)
#set par(justify: true)

= 从 SIMPLE 到投影法：面向有限体积矩阵实现的教学讲义

本文不是从抽象数学开始介绍投影法，而是从你已经熟悉的 SIMPLE 算法出发。
假设你已经知道以下内容：

- 如何对动量方程做有限体积离散，得到类似 $A_u u = b_u - G_x p$、$A_v v = b_v - G_y p$ 的线性系统；
- 如何先用一个已有压力场预测速度；
- 如何根据连续性方程构造压力修正方程；
- 如何用压力修正 $p'$ 修正速度和压力。

在这个基础上，投影法可以被理解为：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*投影法也是一种压力-速度耦合处理方法。它先求一个不一定满足连续性的中间速度 $bold(u)^*$，然后通过求解一个压力增量 Poisson 方程，把 $bold(u)^*$ 修正到满足 $nabla dot bold(u) = 0$ 的速度空间中。*
]

从实现角度看，它和 SIMPLE 长得很像：

$ "预测速度" -> "解压力修正/压力增量方程" -> "修正速度" -> "修正压力" $

但从理论角度看，投影法多了一层更清晰的解释：

$ "把速度场投影到无散度空间" $

本文的目标就是把这两种语言接起来：

- SIMPLE 语言：压力修正、速度修正、通量平衡、矩阵系数；
- 投影法语言：中间速度、压力增量、Poisson 方程、无散度投影；
- 代码语言：`u_star`、`v_star`、`pressure_correction`、`d_u`、`d_v`、face flux residual。


== 0. 先回忆你已经会的 SIMPLE 思路

对不可压缩流体，核心方程是：

$ nabla dot bold(u) = 0 $

$ rho (bold(u) dot nabla) bold(u) = - nabla p + mu nabla^2 bold(u) $

稳态 SIMPLE 的困难在于：

- 动量方程里需要压力梯度；
- 连续性方程里没有显式压力；
- 但压力必须调节速度，使速度最终满足质量守恒。

所以 SIMPLE 通常这样做。

=== 0.1 用当前压力预测速度

假设第 $k$ 次迭代已经有压力 $p^k$。先把压力当成已知量，解动量方程：

$ A_u u^* = b_u - G_x p^k $

$ A_v v^* = b_v - G_y p^k $

这里的 $u^*$、$v^*$ 是预测速度。它们满足动量方程的近似形式，但通常不满足连续性方程。

也就是说，对某个控制体 $P$，预测速度产生的面通量可能满足：

$ F_e^* - F_w^* + F_n^* - F_s^* != 0 $

这表示该控制体的质量流入和质量流出不平衡。

=== 0.2 用压力修正消除质量不平衡

SIMPLE 的核心想法是：速度误差主要来自压力误差。

令：

$ p = p^k + p' $

$ u = u^* + u' $

$ v = v^* + v' $

然后近似认为速度修正主要由压力修正的梯度决定：

$ u' approx - d_u partial_x p' $

$ v' approx - d_v partial_y p' $

这里的 $d_u$、$d_v$ 来自动量方程矩阵的对角系数。例如在很多有限体积 SIMPLE 写法里：

$ d_P approx V_P / a_P $

或者在面上写成 $d_f$。它表达的是：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*如果压力修正有一个梯度，速度会响应多少。*
]

把修正后的速度代入连续性方程，就得到压力修正方程：

$ "divergence of corrected flux" = 0 $

也就是：

$ "predicted mass imbalance" - "pressure correction induced flux change" = 0 $

最后解出 $p'$，更新压力和速度。

这就是你熟悉的 SIMPLE 图像。


== 1. 投影法和 SIMPLE 的第一层对应关系

现在先不要急着看 Helmholtz 分解，也不要急着看 Poisson 方程。先记住一个最重要的对应关系：

#figure(
  table(
    columns: 3,
    align: left,
    inset: 6pt,
    stroke: none,
    table.hline(),
    [*SIMPLE 里的概念*], [*投影法里的概念*], [*直观解释*],
    table.hline(stroke: 0.6pt),
    [预测速度 $bold(u)^*$],
    [中间速度 $bold(u)^*$],
    [先按动量方程走一步，但还没严格满足连续性。],

    [压力修正 $p'$],
    [压力增量 $phi$],
    [用来修正速度，使速度满足 $nabla dot bold(u)=0$。],

    [速度修正 $bold(u)' approx - d nabla p'$],
    [$bold(u)^(n+1) = bold(u)^* - (Delta t/rho) nabla phi$],
    [形式几乎一样，只是响应系数的来源不同。],

    [压力修正方程],
    [压力增量 Poisson 方程],
    [本质上都是由连续性方程推出来的。],

    [通量守恒],
    [无散度投影],
    [有限体积上就是每个压力控制体质量守恒。],
    table.hline(),
  ),
  caption: [SIMPLE 与投影法的概念对应。],
)

所以，对已有 SIMPLE 基础的人来说，投影法并不是完全陌生的算法。它和 SIMPLE 最像的一句话是：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*SIMPLE 用压力修正让预测速度满足连续性；投影法用压力增量把预测速度投影到无散度空间。*
]

区别主要有两个：

1. SIMPLE 通常面向稳态迭代，压力修正方程里的速度响应 $d$ 多来自松弛后的动量矩阵对角项。
2. 投影法通常从瞬态不可压缩 Navier-Stokes 推导出来，速度响应常写成 $Delta t / rho$，或者更一般地写成动量算子的逆。


== 2. 为什么需要“投影”？

先看一个有限体积控制体。

对一个二维网格单元 $P$，连续性方程离散后就是：

$ rho u_e A_e - rho u_w A_w + rho v_n A_n - rho v_s A_s = 0 $

如果是均匀二维网格，可以写成：

$ rho Delta y (u_e - u_w) + rho Delta x (v_n - v_s) = 0 $

这句话的意思很朴素：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*流入多少，必须流出多少。不可压缩流体不允许一个控制体凭空积累质量，也不允许凭空损失质量。*
]

但是动量方程预测出来的 $u^*$、$v^*$ 不一定满足这个条件。也就是说，可能出现：

$ rho Delta y (u_e^* - u_w^*) + rho Delta x (v_n^* - v_s^*) = r_P^m $

其中 $r_P^m$ 是预测速度造成的质量残差。

如果 $r_P^m > 0$，说明该控制体净流出偏多；如果 $r_P^m < 0$，说明净流入偏多。

投影法要做的事情就是：

$ bold(u)^* -> bold(u)^(n+1) $

并且要求：

$ nabla dot bold(u)^(n+1) = 0 $

也就是说，原来的 $bold(u)^*$ 可能“不合法”，因为它不满足不可压缩约束。投影法把它修正为一个“合法”的速度场。

这就是“投影”这个词的来源。它类似于线性代数中的投影：

- 原始向量不在某个子空间里；
- 找到它在这个子空间上的投影；
- 投影后的向量满足子空间约束。

在这里：

- 原始对象：预测速度 $bold(u)^*$；
- 目标子空间：所有满足 $nabla dot bold(u)=0$ 的速度场；
- 投影结果：修正速度 $bold(u)^(n+1)$。


== 3. 从 SIMPLE 公式自然过渡到投影法公式

你熟悉的 SIMPLE 速度修正可以写成：

$ bold(u) = bold(u)^* - d nabla p' $

投影法的速度修正写成：

$ bold(u)^(n+1) = bold(u)^* - (Delta t / rho) nabla phi $

这两个式子结构完全一样。

关键差别是：

#figure(
  table(
    columns: 3,
    align: left,
    inset: 6pt,
    stroke: none,
    table.hline(),
    [*项目*], [*SIMPLE*], [*投影法*],
    table.hline(stroke: 0.6pt),
    [压力修正变量], [$p'$], [$phi$],
    [速度响应系数], [$d approx V_P / a_P$], [$Delta t / rho$，或更一般的 $M^(-1)$],
    [主要来源], [稳态动量方程的压力-速度耦合近似], [瞬态动量方程的时间推进分裂],
    [压力更新], [$p^(k+1)=p^k + alpha_p p'$], [$p^(n+1)=p^n + phi$，有些格式还会加旋转修正],
    [核心解释], [迭代修正压力-速度耦合], [把速度投影到无散度空间],
    table.hline(),
  ),
  caption: [SIMPLE 与投影法公式结构对比。],
)

因此，如果你已经会写 SIMPLE 的压力修正方程，那么理解投影法最重要的不是“会不会构造矩阵”，而是理解：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*投影法中的压力增量 $phi$ 不是额外的经验修正，而是由“修正后的速度必须无散度”这个约束严格推出来的。*
]


== 4. 投影法的连续方程推导

现在从不可压缩 Navier-Stokes 方程出发。

$ partial_t bold(u) + bold(u) dot nabla bold(u) = - (1/rho) nabla p + nu nabla^2 bold(u) $

$ nabla dot bold(u) = 0 $

其中：

$ nu = mu / rho $

为了从时间 $n$ 推进到 $n+1$，最直接的隐式写法是：

$ (bold(u)^(n+1) - bold(u)^n) / Delta t
+ "convection/diffusion terms"
= - (1/rho) nabla p^(n+1) $

并且要求：

$ nabla dot bold(u)^(n+1) = 0 $

这个系统是耦合的，因为 $bold(u)^(n+1)$ 和 $p^(n+1)$ 同时未知。

投影法把它拆成两步。

=== 4.1 第一步：动量预测

先不用未知的新压力 $p^(n+1)$，而是用已知的旧压力 $p^n$，求一个中间速度 $bold(u)^*$：

$ (bold(u)^* - bold(u)^n) / Delta t
+ "convection/diffusion terms"
= - (1/rho) nabla p^n $

这个式子的意思是：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*先让速度按照惯性、对流、扩散、旧压力梯度走一步。走完之后，先不保证它完全满足不可压缩连续性。*
]

所以 $bold(u)^*$ 是预测速度，也叫中间速度。

=== 4.2 第二步：用压力增量修正速度

令压力增量为：

$ phi = p^(n+1) - p^n $

真实的新压力是：

$ p^(n+1) = p^n + phi $

由于预测步骤只用了 $p^n$，还少了 $phi$ 的压力梯度作用。这个缺少的压力作用会造成速度修正：

$ (bold(u)^(n+1) - bold(u)^*) / Delta t = - (1/rho) nabla phi $

移项得到：

$ bold(u)^(n+1) = bold(u)^* - (Delta t / rho) nabla phi $

这就是投影法的速度修正式。

对照 SIMPLE：

$ bold(u) = bold(u)^* - d nabla p' $

你会发现它们几乎一样。只是投影法里这个 $d$ 由时间离散自然给出：

$ d = Delta t / rho $

如果使用更一般的隐式投影，也可以把 $Delta t / rho$ 替换为某个动量响应算子 $M^(-1)$。

=== 4.3 第三步：由连续性推出压力增量方程

修正后的速度必须满足：

$ nabla dot bold(u)^(n+1) = 0 $

把速度修正式代入：

$ nabla dot (bold(u)^* - (Delta t / rho) nabla phi) = 0 $

展开：

$ nabla dot bold(u)^* - (Delta t / rho) nabla^2 phi = 0 $

于是得到：

$ nabla^2 phi = (rho / Delta t) nabla dot bold(u)^* $

这就是压力增量 Poisson 方程。

它和 SIMPLE 压力修正方程的来源完全一致：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*压力增量方程不是为了让压力本身更准，而是为了让压力增量的梯度修正速度后，速度满足连续性。*
]


== 5. 用有限体积语言重新推一遍

上面的推导是连续形式。现在把它翻译成你熟悉的有限体积形式。

=== 5.1 预测通量的质量残差

对控制体 $P$，预测速度给出预测面通量：

$ F_e^*, F_w^*, F_n^*, F_s^* $

连续性要求：

$ F_e - F_w + F_n - F_s = 0 $

但预测通量通常满足：

$ F_e^* - F_w^* + F_n^* - F_s^* = r_P^m $

其中 $r_P^m$ 是预测质量残差。

=== 5.2 压力增量如何改变面通量

投影法速度修正是：

$ bold(u) = bold(u)^* - d nabla phi $

其中：

$ d = Delta t / rho $

在东西方向某个面 $e$ 上，压力增量梯度近似为：

$ (partial phi / partial x)_e approx (phi_E - phi_P) / delta x_e $

所以面速度修正为：

$ u_e = u_e^* - d_e (phi_E - phi_P) / delta x_e $

对应的质量通量修正是：

$ F_e = F_e^* - rho A_e d_e (phi_E - phi_P) / delta x_e $

同理，西、北、南各个面都有类似修正。

=== 5.3 代入连续性，得到压力增量矩阵

把修正后的面通量代入：

$ F_e - F_w + F_n - F_s = 0 $

就得到一个关于 $phi_P$ 及其邻居 $phi_E, phi_W, phi_N, phi_S$ 的线性方程：

$ a_P phi_P - a_E phi_E - a_W phi_W - a_N phi_N - a_S phi_S = b_P $

其中右端项 $b_P$ 就来自预测质量不平衡，例如可以写成：

$ b_P = F_w^* - F_e^* + F_s^* - F_n^* $

或者根据你代码里的符号约定写成：

$ m_P^* = rho Delta y (u_w^* - u_e^*) + rho Delta x (v_s^* - v_n^*) $

邻居系数则来自压力增量对面通量的影响。例如东侧：

$ a_E = rho A_e d_e / delta x_e $

均匀网格二维情形下，如果 $A_e = Delta y$，$delta x_e = Delta x$，则：

$ a_E = rho Delta y d_e / Delta x $

由于投影法中：

$ d_e = Delta t / rho $

所以：

$ a_E = Delta t Delta y / Delta x $

如果把方程两边整体乘以某些体积或时间系数，形式可能略有不同，但物理意义一样：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*压力增量矩阵的系数，表示“相邻两个控制体之间的压力增量差，会造成多少穿过该面的通量修正”。*
]

这和 SIMPLE 中压力修正矩阵的含义完全一致。


== 6. 从矩阵角度看投影法

如果你习惯构造矩阵，投影法可以写得非常清楚。

定义：

- $bold(U)^*$：预测速度向量，或者预测面通量向量；
- $Phi$：所有单元上的压力增量未知量组成的向量；
- $G$：离散梯度算子，把压力增量变成速度修正方向；
- $D$：离散散度算子，把速度或面通量变成每个控制体的质量残差；
- $M$：速度响应矩阵。

速度修正可以写成：

$ bold(U)^(n+1) = bold(U)^* - M^(-1) G Phi $

连续性约束是：

$ D bold(U)^(n+1) = 0 $

代入速度修正：

$ D (bold(U)^* - M^(-1) G Phi) = 0 $

得到：

$ D M^(-1) G Phi = D bold(U)^* $

这就是投影法最重要的离散矩阵形式。

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
如果你会构造 SIMPLE 的压力修正矩阵，那么可以把投影法的压力增量矩阵理解成：

$ "压力增量矩阵" = "散度算子" times "速度响应" times "梯度算子" $

也就是：

$ L = D M^(-1) G $
]

在最简单的一阶投影法中，$M$ 只是时间项带来的质量矩阵，所以：

$ M approx (rho / Delta t) I $

于是：

$ M^(-1) approx (Delta t / rho) I $

这就是为什么速度修正里会出现 $Delta t / rho$。

而在 SIMPLE 中，$M^(-1)$ 通常被近似成动量矩阵对角项的倒数：

$ M^(-1) approx 1 / a_P $

所以二者的“骨架”是一样的，只是速度响应的近似来源不同。


== 7. “压力”在不可压缩流里到底是什么？

理解投影法时，一个常见障碍是：为什么压力可以由一个 Poisson 方程求出来？

在可压缩流里，压力常常和密度、温度通过状态方程联系起来。例如理想气体：

$ p = rho R T $

但在恒密度不可压缩流里，密度已经固定，不能用密度变化来决定压力。

此时压力的主要角色是：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*压力是维持 $nabla dot bold(u)=0$ 这个约束的约束力，或者说 Lagrange multiplier。*
]

这和约束力的例子很像。比如一个小球被限制在某个曲面上运动，约束反力本身不是由独立的运动方程给出的，而是由“小球必须留在曲面上”这个约束决定的。

不可压缩流里的压力也类似：

- 对流、扩散、外力会试图改变速度；
- 这些变化可能产生非零散度；
- 压力梯度立刻调整速度，使最终速度仍然无散度。

所以投影法中的 Poisson 方程可以理解为：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*求一个压力增量 $phi$，让它的梯度正好抵消预测速度中的散度误差。*
]


== 8. Helmholtz-Hodge 分解：投影法的理论背景

这一节稍微抽象一点，但它能解释为什么算法叫“投影法”。

Helmholtz-Hodge 分解说，在合适的边界条件下，一个速度场可以分解成两部分：

$ bold(w) = bold(w)_perp + nabla q $

其中：

$ nabla dot bold(w)_perp = 0 $

也就是说，任意速度场 $bold(w)$ 可以拆成：

- 一个无散度部分 $bold(w)_perp$；
- 一个梯度部分 $nabla q$。

投影法做的事情就是：

$ bold(u)^* = bold(u)^(n+1) + (Delta t / rho) nabla phi $

或者写成：

$ bold(u)^(n+1) = bold(u)^* - (Delta t / rho) nabla phi $

这里：

- $bold(u)^*$ 是普通预测速度场；
- $(Delta t/rho) nabla phi$ 是需要去掉的梯度部分；
- $bold(u)^(n+1)$ 是剩下的无散度部分。

所以从这个角度看：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*压力增量不是为了“修压力”而存在，而是为了构造一个梯度场。把这个梯度场从预测速度里减掉，剩下的速度就满足不可压缩约束。*
]

这就是投影法最核心的理论图像。


== 9. 投影也可以看成一个最小修正问题

还有一个很有帮助的理解方式。

预测速度 $bold(u)^*$ 已经满足了动量方程的大部分信息：对流、扩散、边界驱动、旧压力梯度等。我们不希望大幅改变它，只希望修正掉它的散度。

所以可以提出一个问题：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
在所有满足 $nabla dot bold(u)=0$ 的速度场中，找一个离 $bold(u)^*$ 最近的速度场。
]

数学上可以写成：

$ min_(bold(u)) (1/2) integral_Omega (rho / Delta t) abs(bold(u) - bold(u)^*)^2 d Omega $

约束是：

$ nabla dot bold(u) = 0 $

对这个约束优化问题引入 Lagrange multiplier $phi$，可以得到：

$ (rho / Delta t) (bold(u) - bold(u)^*) + nabla phi = 0 $

于是：

$ bold(u) = bold(u)^* - (Delta t / rho) nabla phi $

再代入连续性：

$ nabla^2 phi = (rho / Delta t) nabla dot bold(u)^* $

这说明投影法有一个非常漂亮的含义：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*投影法不是随便修速度。它是在满足不可压缩约束的前提下，对预测速度做尽可能小的修正。*
]

这个解释对数值实现很重要。因为你会看到，投影法不是把动量预测完全推翻，而是只减去一个梯度场。


== 10. 投影法中的边界条件怎么理解？

边界条件是投影法最容易混乱的地方。

从速度修正式开始：

$ bold(u)^(n+1) = bold(u)^* - (Delta t / rho) nabla phi $

取边界法向分量：

$ bold(u)^(n+1) dot bold(n)
= bold(u)^* dot bold(n) - (Delta t / rho) partial_n phi $

如果边界要求法向速度为 $g_n$，即：

$ bold(u)^(n+1) dot bold(n) = g_n $

那么压力增量的法向导数应该满足：

$ partial_n phi = (rho / Delta t) (bold(u)^* dot bold(n) - g_n) $

对于封闭 cavity 的固壁：

$ g_n = 0 $

如果预测速度在壁面已经处理成没有法向穿透：

$ bold(u)^* dot bold(n) = 0 $

那么：

$ partial_n phi = 0 $

这就是常见的压力增量零法向梯度边界条件。

=== 10.1 moving lid 会不会改变这个条件？

lid-driven cavity 顶盖有切向速度，但没有法向速度。

比如上边界：

$ u = U_"lid", quad v = 0 $

这里 $u$ 是切向，$v$ 是法向。

压力增量 Neumann 条件主要来自法向速度约束。所以 moving lid 不会把法向压力增量条件变成非零。只要预测步骤已经保证边界法向速度为零，就仍然可以使用：

$ partial_n phi = 0 $

=== 10.2 为什么修正后还要重新施加速度边界？

投影修正主要保证连续性，尤其是法向通量平衡。它不一定自动保证所有切向 no-slip 条件仍然精确满足。

所以实际有限体积代码中常见做法是：

1. 预测前施加速度边界；
2. 预测后必要时处理 ghost cells；
3. 压力增量修正后再次施加速度边界；
4. 用修正后的速度重建 face velocity；
5. 再计算连续性残差。

这不是多余操作，而是为了保证投影步骤和边界条件相容。


== 11. 压力增量方程为什么需要钉住一个参考值？

不可压缩流中，速度只受压力梯度影响。

如果把压力整体加一个常数：

$ p -> p + C $

那么：

$ nabla p -> nabla p $

速度完全不变。

压力增量也一样。如果：

$ phi -> phi + C $

那么：

$ nabla phi -> nabla phi $

速度修正不变。

所以纯 Neumann 压力增量 Poisson 方程有一个常数零空间。离散矩阵会奇异。

为了解出唯一的代数解，通常需要：

- 固定一个单元的 $phi$，例如 $phi_(1,1)=0$；或
- 强制所有 $phi$ 的平均值为零；或
- 使用能处理 nullspace 的线性求解器。

这并不是物理上指定了某个绝对压力，而只是为了消除代数矩阵的常数自由度。


== 12. 物理时间步与伪时间步

投影法最早常用于非稳态问题，所以公式里有 $Delta t$。

如果你真的在算瞬态流动，那么：

$ Delta t = "真实物理时间步" $

它决定从 $t^n$ 到 $t^(n+1)$ 的物理演化。

但如果你用投影法求稳态 cavity，很多代码会把它当作伪时间推进：

$ "steady solution" = "pseudo-time marching 的收敛极限" $

这时 $Delta t$ 主要是数值参数：

- 小 $Delta t$：时间项强，矩阵更稳，但收敛可能慢；
- 大 $Delta t$：推进更激进，可能更快，也可能不稳定；
- 如果最终收敛到同一个离散稳态解，$Delta t$ 不应该改变最终物理解，只改变收敛路径。

在这种稳态伪时间投影中，动量预测方程可写成：

$ (rho / Delta t) (bold(u)^* - bold(u)^k)
+ C(bold(u)^k, bold(u)^*)
- D(bold(u)^*)
= - nabla p^k $

其中 $k$ 是迭代步，不一定是物理时间层。

收敛时：

$ bold(u)^(k+1) approx bold(u)^k $

于是伪时间项趋于很小，剩下的就是稳态动量平衡和连续性方程。


== 13. 和 SIMPLE 的更细致对比

从代码实现看，SIMPLE 和投影法很容易写成同一套框架。它们差异主要不在矩阵装配能力，而在算法解释和速度响应系数。

#figure(
  table(
    columns: 3,
    align: left,
    inset: 6pt,
    stroke: none,
    table.hline(),
    [*问题*], [*SIMPLE*], [*投影法*],
    table.hline(stroke: 0.6pt),

    [预测速度怎么来？],
    [用当前压力解稳态或伪稳态动量方程。],
    [用旧压力推进一个时间步或伪时间步，得到中间速度。],

    [压力修正方程怎么来？],
    [把速度修正式代入连续性。],
    [把投影修正式代入 $nabla dot bold(u)^(n+1)=0$。],

    [速度响应系数是什么？],
    [通常来自动量方程对角项和欠松弛。],
    [最简单形式是 $Delta t/rho$；更一般是动量算子的逆或近似逆。],

    [压力更新是否欠松弛？],
    [通常需要 $alpha_p$ 欠松弛。],
    [瞬态投影一般直接增量更新；伪稳态实现也可加入松弛。],

    [理论解释是什么？],
    [压力-速度耦合迭代。],
    [Helmholtz-Hodge 投影到无散度速度空间。],

    [最终检查什么？],
    [连续性残差、动量残差、变量变化量。],
    [同样检查连续性残差、动量残差、变量变化量。],
    table.hline(),
  ),
  caption: [SIMPLE 与投影法的算法差异。],
)

可以说：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*SIMPLE 是从稳态压力-速度耦合迭代的角度理解压力修正；投影法是从瞬态/伪瞬态速度空间投影的角度理解压力增量。*
]

二者都可以落到相似的有限体积压力修正矩阵上。


== 14. 投影法在当前代码中的变量对应

如果当前代码有如下变量，可以这样理解：

#figure(
  table(
    columns: 3,
    align: left,
    inset: 6pt,
    stroke: none,
    table.hline(),
    [*代码变量*], [*数学符号*], [*含义*],
    table.hline(stroke: 0.6pt),
    [`u`, `v`], [$u^k, v^k$ 或 $u^(k+1), v^(k+1)$], [当前修正后的速度场。],
    [`u_star`, `v_star`], [$u^*, v^*$], [动量预测得到的中间速度。],
    [`pressure`], [$p^k$ 或 $p^(k+1)$], [当前压力场。],
    [`pressure_correction`], [$phi$], [压力增量，也就是投影修正势函数。],
    [`d_u`, `d_v`], [$d = Delta t/rho$], [压力增量梯度造成的速度响应。],
    [`u_face`, `v_face`], [face velocity/flux], [用于对流、压力增量方程和连续性残差的面速度。],
    table.hline(),
  ),
  caption: [投影法变量与代码变量对应。],
)

投影法一次迭代可以理解成：

#figure(
  table(
    columns: 3,
    align: left,
    inset: 6pt,
    stroke: none,
    table.hline(),
    [*阶段*], [*数学含义*], [*代码意图*],
    table.hline(stroke: 0.6pt),
    [准备], [施加边界，准备 $d = Delta t/rho$], [设置 `d_u`, `d_v`，更新 ghost cells。],
    [预测 $u^*$], [解 x 动量预测方程], [得到 `u_star`。],
    [预测 $v^*$], [解 y 动量预测方程], [得到 `v_star`。],
    [构造压力增量方程], [由预测 face flux 的质量残差构造 $L Phi = b$], [填充 `pressure_correction` 线性系统。],
    [求解 $phi$], [得到要减去的梯度场], [线性求解器求 `pressure_correction`。],
    [修正速度], [$u = u^* - d partial_x phi$, $v = v^* - d partial_y phi$], [更新 `u`, `v`。],
    [修正压力], [$p = p + phi$], [更新 `pressure` 并去掉参考压力。],
    [检查收敛], [检查动量和连续性], [重建 face velocity，计算 residual。],
    table.hline(),
  ),
  caption: [投影法一次迭代流程。],
)


== 15. 动量预测矩阵如何理解

投影法的动量预测和你熟悉的动量方程装配基本一样，只是多了一个时间项或伪时间项。

对 $u$ 方程，有限体积离散后可以写成：

$ a_P u_P^* - a_E u_E^* - a_W u_W^* - a_N u_N^* - a_S u_S^* = b_P^u $

如果使用伪时间项，则中心系数多出：

$ a_P^"time" = rho V_P / Delta t $

右端项多出：

$ b_P^"time" = (rho V_P / Delta t) u_P^k $

因此：

$ a_P^"proj" = a_P^"base" + rho V_P / Delta t $

$ b_P^"proj" = b_P^"base" + (rho V_P / Delta t) u_P^k + "pressure source from" p^k $

这里的含义是：

- 对流、扩散项仍然像普通有限体积动量方程那样装配；
- 旧速度 $u^k$ 通过伪时间源项进入右端；
- 当前压力 $p^k$ 的梯度作为已知源项进入右端；
- 新压力增量 $phi$ 不在预测方程里，而在后面的投影步骤里处理。

这与 SIMPLE 的预测动量方程非常接近。


== 16. 压力增量矩阵如何理解

压力增量矩阵的核心原则只有一个：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*用同一套面通量修正公式来构造压力增量方程，并用同一套面通量定义来检查连续性残差。*
]

为什么这句话重要？

因为有限体积连续性是通过面通量表达的。如果你用一种 face velocity 构造压力增量矩阵，却用另一种 face velocity 检查 residual，就可能出现：

- 压力增量线性系统看起来解得很好；
- 但最终质量残差仍然不小；
- 或者 collocated grid 上出现压力-速度棋盘格解耦。

因此 collocated 网格上通常需要 Rhie-Chow 类型的面速度重构。它的目标不是改变控制方程，而是让压力和速度在同位网格上保持耦合，避免奇偶振荡。

在投影法中，压力增量方程可以记成：

$ L Phi = b $

其中：

$ L = D M^(-1) G $

右端项：

$ b = D bold(U)^* $

也就是预测速度的散度或预测面通量的质量不平衡。


== 17. 常见误解

=== 17.1 “投影法是不是不用压力修正？”

不是。

投影法当然有压力修正，只是通常叫 pressure increment 或 pressure correction，记作 $phi$。

它和 SIMPLE 的 $p'$ 很像，只是投影法更强调：$phi$ 的梯度用于把速度投影到无散度空间。

=== 17.2 “投影法是不是只能做瞬态，不能做稳态？”

不是。

投影法可以做瞬态，也可以作为伪时间迭代方法求稳态。做稳态时，时间步是数值推进参数，不一定代表真实物理时间。

=== 17.3 “压力 Poisson 方程是不是压力的物理方程？”

更准确地说，它是由动量方程和连续性约束导出的约束方程。

不可压缩流中压力不是由状态方程独立决定，而是为了维持无散度速度场而出现。

=== 17.4 “只要 cell-centered 速度修正了，连续性一定满足吗？”

不一定。

有限体积连续性真正检查的是 face flux balance。因此实现中更关键的是：

- 用预测 face flux 构造压力增量方程；
- 用压力增量修正 face flux；
- 用修正后的 face flux 检查连续性。

cell-centered velocity 对下一轮动量方程重要，但质量守恒主要体现在 face flux 上。

=== 17.5 “压力增量边界条件随便设成零可以吗？”

不能随便设。

零法向梯度 $partial_n phi = 0$ 通常成立的前提是：

- 边界是不可穿透壁面；
- 预测速度已经满足边界法向速度；
- 你希望修正后仍然保持这个法向速度。

如果入口、出口、开放边界或指定通量边界存在，压力增量边界条件要重新推导。


== 18. 你可以按这个顺序学习投影法

如果你已经掌握 SIMPLE，我建议按下面顺序建立投影法知识结构。

=== 第一步：先把它看成 SIMPLE 的近亲

记住：

$ bold(u) = bold(u)^* - d nabla p' $

对应：

$ bold(u)^(n+1) = bold(u)^* - (Delta t/rho) nabla phi $

先不要被“投影”这个名字吓到。实现上就是：预测速度、构造压力增量方程、修正速度和压力。

=== 第二步：理解 $d$ 的来源不同

SIMPLE 中：

$ d approx 1/a_P $

投影法最简单形式中：

$ d = Delta t/rho $

更一般地：

$ d " comes from " M^(-1) $

这一步能把矩阵装配和时间推进联系起来。

=== 第三步：理解压力增量方程就是连续性

投影法的 Poisson 方程不是额外假设，而是：

$ nabla dot (bold(u)^* - d nabla phi) = 0 $

也就是：

$ D M^(-1) G Phi = D bold(U)^* $

只要你能从 face flux balance 构造 SIMPLE 压力修正矩阵，就能构造投影法压力增量矩阵。

=== 第四步：理解“投影”的数学意义

最后再看 Helmholtz-Hodge 分解：

$ bold(u)^* = "divergence-free part" + "gradient part" $

投影法就是减去 gradient part，留下 divergence-free part。

这一步会让你明白：

- 为什么压力增量是一个标量场；
- 为什么要解 Poisson 方程；
- 为什么压力只差一个常数；
- 为什么边界条件非常重要。


== 19. 最小伪代码

下面用偏实现的方式总结投影法。

#block(fill: luma(248), inset: 8pt, radius: 4pt)[
```
Given u^k, v^k, p^k

1. Apply velocity and pressure boundary conditions.

2. Set projection response:
      d_u = d_v = projection_dt / rho

3. Solve momentum predictor:
      A_u u_star = b_u - G_x p^k + time_source(u^k)
      A_v v_star = b_v - G_y p^k + time_source(v^k)

4. Reconstruct predicted face velocities / fluxes from u_star, v_star.

5. Assemble pressure increment equation from predicted mass imbalance:
      L phi = D U_star
   where:
      L = D M^{-1} G

6. Pin pressure increment reference:
      phi_ref = 0

7. Solve for phi.

8. Correct velocity:
      u^{k+1} = u_star - d_u * grad_x(phi)
      v^{k+1} = v_star - d_v * grad_y(phi)

9. Correct pressure:
      p^{k+1} = p^k + phi
      p^{k+1} = p^{k+1} - p_reference

10. Reapply boundary conditions.

11. Reconstruct corrected face velocities and compute residuals.

12. Stop if continuity residual, momentum residual, and update norms are small.
```
]


== 20. 一句话总结

对已经理解 SIMPLE 的人，投影法可以这样记：

#block(fill: luma(245), inset: 8pt, radius: 4pt)[
*SIMPLE 和投影法都通过压力修正来消除预测速度的质量不平衡。SIMPLE 更像稳态压力-速度耦合迭代；投影法更像先按动量方程走一步，再用压力增量把速度场正交投影到无散度空间。有限体积实现中，二者最终都落到“用压力修正改变面通量，使每个控制体满足质量守恒”的矩阵构造上。*
]
