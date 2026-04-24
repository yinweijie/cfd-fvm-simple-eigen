# 架构图阅读指南

本文用于帮助阅读 `docs/` 下的三张 Mermaid 架构图：

- 模块关系图：[src_code_structure.mmd](src_code_structure.mmd)，渲染图：[src_code_structure.svg](src_code_structure.svg)
- SIMPLE 流程图：[simple_solver_flow.mmd](simple_solver_flow.mmd)，渲染图：[simple_solver_flow.svg](simple_solver_flow.svg)
- 数据归属图：[data_ownership.mmd](data_ownership.mmd)，渲染图：[data_ownership.svg](data_ownership.svg)

这三张图不是同一件事的重复表达。它们分别回答三个问题：

- `src_code_structure`：源码模块之间谁调用谁、谁依赖哪些核心数据。
- `simple_solver_flow`：一次完整求解过程中，尤其每轮 SIMPLE 迭代，按什么顺序发生。
- `data_ownership`：关键数据由谁创建、谁拥有、谁读、谁写，以及临时对象如何流动。

推荐阅读顺序是：先看模块关系图建立全局地图，再看 SIMPLE 流程图理解运行行为，最后看数据归属图确认生命周期和状态修改边界。

## 节点文字规则

为了避免把节点里的多行文字误读成多个模块，三张图采用两种不同的标签规则。

模块关系图 `src_code_structure` 使用“模块分组 + 职责节点”的写法：

- 黄色分组标题里的 `Module:` 或 `Header:` 表示源码文件或接口文件。
- 分组内部的 `R1`、`S1`、`D1`、`L1`、`O1` 是职责节点。
- 连线标签会重复对应职责编号和职责文字，例如 `R1 parse CLI into CavityCase`。
- 因此，一条线上的文字应该能直接回到同名职责节点；读者不需要再从一个长 `Role:` 描述里猜是哪一项职责。

数据归属图 `data_ownership` 仍采用“字段名 + 内容”的写法：

- `Module:` 表示实现文件，例如 `src/main.cpp` 或 `src/discretization.cpp`。
- `Header:` 表示主要接口头文件，例如 `include/cfd/field.hpp`。
- `Type:` 表示该文件中最重要的类型，例如 `SimpleSolver` 或 `FlowFields`。
- `Object:` / `Member:` / `Stack object:` 表示运行时对象或成员变量。
- `Role:` 表示这个模块或 helper 的职责。
- `Data:` / `Contains:` / `Names:` 表示该节点承载的数据内容。
- `Output:` / `Files:` 表示文件系统输出。

也就是说，模块关系图优先表达“哪个职责造成哪条关系”，数据归属图优先表达“这个节点是谁、里面有什么”。

## 总体行为

程序入口是 `src/main.cpp`。它解析命令行参数，得到一个 `CavityCase config`，然后构造 `SimpleSolver`。`SimpleSolver` 内部拥有三类核心状态：

- `config_`：求解参数、物理参数、松弛因子和输出路径。
- `grid_`：结构化网格尺寸、单元大小、二维索引到线性系统行号的映射。
- `fields_`：压力、速度、预测速度、压力修正、速度修正因子和面速度。

求解时，`SimpleSolver::run()` 反复执行 SIMPLE 迭代。每轮迭代先刷新边界条件和面速度，然后依次求解 `u` 动量方程、`v` 动量方程和压力修正方程，再用压力修正更新 `u / v / pressure`，最后计算连续性残差和收敛指标。

求解完成后，`main.cpp` 调用 `write_results()`，把最终场变量、中心线数据、残差历史和摘要写到 `results/<case>/`。

## 图 1：模块关系图

文件：

- [src_code_structure.mmd](src_code_structure.mmd)
- [src_code_structure.svg](src_code_structure.svg)

这张图是最高层的源码地图。它把项目分成四块：

- `Entry Point`：`src/main.cpp`，负责 CLI、配置和顶层编排。
- `Solver Runtime`：`SimpleSolver`、离散化函数和线性求解封装。
- `Core Data Model`：`CavityCase`、`StructuredGrid`、`FlowFields` 和 Eigen 数据结构。
- `Result Export`：`write_results()` 和最终输出文件。

这张图里的箭头不要全部理解成 `#include` 或编译依赖。它表达的是代码行为层面的关系。

模块关系图最重要的读法是：先找模块分组，再看分组内的职责节点，最后按职责编号追踪连线。

- `R*`：`src/main.cpp` 的职责。
- `S*`：`src/simple_solver.cpp` / `SimpleSolver` 的职责。
- `D*`：`src/discretization.cpp` 的职责。
- `L*`：`include/cfd/linear_system.hpp` 的职责。
- `O*`：`src/output.cpp` 的职责。

线条含义如下：

| 线条 | Mermaid 写法 | 含义 |
| --- | --- | --- |
| 实线箭头 | `A --> B` | 职责 `A` 主动调用、构造、更新或产出 `B`。 |
| 点线箭头 | `A -.-> B` | 职责 `A` 只读使用 `B`。 |
| 粗线箭头 | `A ==> B` | 职责 `A` 拥有 `B`，或把 `B` 作为产物导出。 |

几个关键读法：

- `R1` 表示 `src/main.cpp` 负责把 CLI 解析成 `CavityCase`。
- `R2` 和 `R3` 分开表达构造 solver 与运行 solver，避免把两个行为塞进一条线里。
- `S1` 表示 `SimpleSolver` 拥有 `config_ / grid_ / fields_`，所以它连到 `Case / Grid / Fields` 的线是粗线。
- `S2` 表示 `SimpleSolver` 调用 SIMPLE 各阶段 helper，所以它连到 `D1` 到 `D5`。
- `D1` 到 `D5` 分别表示边界、面速度、动量装配、压力修正装配、连续性残差；这些职责的连线标签与职责节点文字保持一致。
- `L1` 表示线性系统求解，它读取 `LinearSystem`，并使用 Eigen 迭代求解器。
- `O1` 表示输出模块读取最终求解状态，`O2` 表示真正写出结果文件。

## 图 2：SIMPLE 流程图

文件：

- [simple_solver_flow.mmd](simple_solver_flow.mmd)
- [simple_solver_flow.svg](simple_solver_flow.svg)

这张图使用 `sequenceDiagram`，重点是时间顺序。它最适合回答“运行时先做什么、后做什么”。

线条含义如下：

| 线条 | Mermaid 写法 | 含义 |
| --- | --- | --- |
| 调用箭头 | `A->>B: message` | `A` 调用 `B` 或请求 `B` 执行一步操作。 |
| 返回箭头 | `A-->>B: message` | `A` 向 `B` 返回结果。 |
| 循环块 | `loop ... end` | 这一段会重复执行，直到满足收敛条件或达到最大迭代数。 |
| 自调用 | `Solver->>Solver` | 对象内部完成的逻辑，例如收集指标并检查收敛。 |

这张图可以按下面的行为读：

1. `main.cpp` 构造 `SimpleSolver`。
2. 构造阶段先初始化边界条件和面速度，使初始场处于一致状态。
3. `main.cpp` 调用 `run()`。
4. 每轮 SIMPLE 迭代开始时，求解器刷新边界和 corrected face velocities。
5. 装配并求解 `u` 动量方程，得到 `u_star`，并用动量方程对角线更新 `d_u`。
6. 用 `u_star` 和当前 `v` 刷新 predictor face fluxes。
7. 装配并求解 `v` 动量方程，得到 `v_star`，并更新 `d_v`。
8. 装配并求解压力修正方程，得到 `pressure_correction`。
9. 用压力修正更新 `u / v / pressure`。
10. 再次应用边界条件并刷新 corrected face velocities。
11. 计算连续性残差，收集本轮指标并检查收敛。
12. 求解结束后返回 `SolveSummary`，再由 `write_results()` 写出结果文件。

注意：这张图强调顺序，不强调所有权。比如 `Fields` 出现在很多步骤里，是因为它不断被读写，不代表每个调用方都拥有它。

## 图 3：数据归属图

文件：

- [data_ownership.mmd](data_ownership.mmd)
- [data_ownership.svg](data_ownership.svg)

这张图用来回答数据生命周期问题：某个对象在哪里创建，谁持有它，谁只是读取它，谁会修改它。

线条含义如下：

| 线条 | Mermaid 写法 | 含义 |
| --- | --- | --- |
| 实线箭头 | `A --> B` | `A` 创建、传递或返回 `B`。 |
| 粗线箭头 | `A ==> B` | `A` 拥有、包含或导出 `B`。 |
| 点线箭头 | `A -.-> B` | `A` 读取 `B`，但不修改所有权。 |
| 虚线/点划箭头 | `A -.- B` | `A` 修改 `B` 的可变状态。 |

这张图中几个重要边界：

- `main()` 栈上先有 `CavityCase config`，它来自 CLI 解析。
- 构造 `SimpleSolver` 时，`config` 被复制进 `config_`。
- `SimpleSolver` 拥有 `config_`、`grid_` 和 `fields_`。
- `FlowFields` 包含三类矩阵状态：cell fields、correction factors、face fields。
- `discretization.cpp` 函数本身是无状态 helper。它读取 `config_ / grid_ / fields_`，也会修改 `fields_`。
- `MomentumAssembly`、`PressureCorrectionAssembly`、`LinearSolveResult`、`IterationMetrics` 都是每轮迭代中的临时对象。
- `LinearSolveResult` 的 solution vector 会被 scatter 回 `FlowFields`。
- `MomentumAssembly` 的 relaxed diagonal 会更新 `d_u / d_v`。
- `SolveSummary` 收集每轮 `IterationMetrics`，最后返回给 `main()`。
- `write_results()` 读取最终状态和 `SolveSummary`，把结果导出到文件系统。

这张图里最应该关注的是“谁会改 `FlowFields`”。主要有三类修改路径：

- 边界和面速度刷新：`discretization.cpp` 更新 ghost cells、`u_face`、`v_face`。
- 线性求解结果回填：`SimpleSolver` 把 `u_star`、`v_star`、`pressure_correction` 写回 `fields_`。
- 校正阶段：`SimpleSolver` 更新 corrected `u`、`v` 和 `pressure`。

## 三张图如何互补

如果你想定位一个函数属于哪个模块，看 `src_code_structure`。

如果你想知道 SIMPLE 的执行先后顺序，看 `simple_solver_flow`。

如果你想知道一个对象的生命周期和修改边界，看 `data_ownership`。

同一个箭头在不同图里可能服务于不同目的。例如：

- 模块图里的 `Solver --> Disc` 是“求解器调用离散化模块”。
- 流程图里的 `Solver->>Disc: assemble_u_momentum(...)` 是“这个调用发生在本轮迭代的具体位置”。
- 数据归属图里的 `Discretization -.- Fields` 是“离散化函数会修改 `FlowFields` 的可变状态”。

因此，不要把三张图里的所有箭头混合成一种单一语义。读图时先确认当前图回答的是结构、顺序，还是所有权。

## 常见误读

这些图不表示完整的 C++ include 图。实际头文件包含关系可以从源码和构建系统中查看。

这些图不表示线程关系。当前求解流程是单进程、顺序执行的数值算法流程。

这些图不表示每个函数的完整调用树。为了可读性，图中只保留了影响理解架构和数据流的主要关系。

这些图不表示数值公式的全部细节。SIMPLE、Rhie-Chow、压力修正方程和边界项推导见 [simple_solver_theory.md](simple_solver_theory.md)。
