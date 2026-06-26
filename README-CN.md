# Mini Stability Theory（迷你稳定性理论）

**从零开始、零依赖的 C 语言实现**，涵盖动力系统稳定性理论的完整体系——从 Lyapunov 直接法（1892）到现代非线性稳定性框架，包括输入-状态稳定性（Sontag 1989）、Lur'e 系统绝对稳定性、有限时间收敛以及互联系统稳定性分析。每个子模块对应 MIT、Stanford、UC Berkeley、Cambridge 等顶尖大学的控制理论课程。

## 子模块总览

| 子模块 | 主题 | 参考课程 |
|--------|------|----------|
| [mini-lyapunov-direct-method](mini-lyapunov-direct-method/) | Lyapunov 稳定性定理、正定函数、二次型/逆定理/非线性 Lyapunov 函数、Lyapunov 方程、全局渐近稳定性、径向无界性 | MIT 6.241J, Stanford AA 203, Berkeley EE 222 |
| [mini-lasalle-invariance](mini-lasalle-invariance/) | LaSalle 不变性原理、Barbashin-Krasovskii 定理、ω 极限集、不变集、吸引域估计、收敛性分析、离散时间 LaSalle | MIT 6.241J, Berkeley EE 222, Caltech CDS 140 |
| [mini-input-to-state-stability](mini-input-to-state-stability/) | ISS 框架（Sontag 1989）、ISS Lyapunov 函数、K/KL/K∞ 类函数、非线性小增益定理、互联系统 ISS | MIT 6.241J, Berkeley EE 222, Stanford AA 273 |
| [mini-input-output-stability](mini-input-output-stability/) | BIBO 稳定性、L2 增益、小增益定理、无源性（耗散性）、锥形扇区、Zames-Falb 乘子、圆判据、Nyquist 分析、结构奇异值 | MIT 6.241J, Caltech CDS 140, Cambridge 4F3 |
| [mini-absolute-stability](mini-absolute-stability/) | Lur'e 系统理论、Aizerman 猜想与反例、Kalman-Yakubovich-Popov（KYP）引理、扇区非线性、基于 Lyapunov 的绝对稳定性检验 | MIT 6.241J, Berkeley EE 222, Cambridge 4F3 |
| [mini-circle-popov-criterion](mini-circle-popov-criterion/) | 圆判据（时变非线性）、Popov 判据（时不变非线性）、回路变换、KYP 引理基础、Nyquist/Popov 图、增益与相位裕度估计 | MIT 6.241J, Cambridge 4F3, ETH 227-0216 |
| [mini-finite-time-stability](mini-finite-time-stability/) | 有限时间收敛（Bhat & Bernstein 2000）、固定时间稳定性（Polyakov 2012）、齐次系统、收敛时间估计、滑模面设计、轨迹分析 | MIT 6.832, Stanford AA 203, UC San Diego MAE 281A |
| [mini-stability-of-interconnected](mini-stability-of-interconnected/) | 复合系统稳定性、分散化 Lyapunov 函数、M 矩阵检验、基于无源性的互联（Hill-Moylan）、互联系统小增益、向量 Lyapunov 函数 | MIT 6.241J, Berkeley EE 222, ETH 227-0216 |

## 设计理念

- **从零开始、零依赖的纯 C 实现** — 每个模块仅使用 C 标准库（`math.h`、`stdlib.h`、`stdbool.h`），不依赖任何外部数值库、无 LAPACK/BLAS
- **模块自包含** — 每个 `mini-*` 目录都是独立项目，拥有自己的 `include/`、`src/`、`test/` 和构建基础设施
- **理论到代码的完整映射** — 头文件内嵌公式推导、关键定理（Lyapunov 1892、LaSalle 1960、Sontag 1989、Popov 1961、Zames 1966）及教材引用（Khalil、Vidyasagar、Isidori）
- **课程对齐的结构** — 每个子模块的 API 按照标准研究生课程大纲组织（MIT 6.241J 动力系统与控制、Berkeley EE 222 非线性系统、Stanford AA 273 非线性控制），可直接作为讲义配套代码使用

## 构建方式

每个子模块相互独立。使用任意 C11+ 编译器即可构建：

```bash
cd mini-lyapunov-direct-method
make
./build/test_runner
```

或手动编译：

```bash
gcc -std=c11 -O2 -Iinclude -o build/test_ldm src/*.c test/*.c -lm
./build/test_ldm
```

需要 **GCC**、**Clang** 或 **MSVC**，支持 C11 标准，链接数学库 `-lm`。

## 项目结构

```
6. mini-stability-theory/
├── mini-absolute-stability/            # Lur'e 系统、Aizerman 猜想、KYP 引理、扇区非线性
├── mini-circle-popov-criterion/        # 圆判据、Popov 判据、回路变换、Nyquist 分析
├── mini-finite-time-stability/         # 有限时间收敛、固定时间稳定、齐次系统
├── mini-input-output-stability/        # BIBO 稳定性、小增益定理、无源性、L2 增益、锥形扇区
├── mini-input-to-state-stability/      # ISS 框架、ISS Lyapunov 函数、非线性小增益
├── mini-lasalle-invariance/            # LaSalle 不变性原理、Barbashin-Krasovskii 定理、不变集、吸引域
├── mini-lyapunov-direct-method/        # Lyapunov 稳定性定理、二次型/逆定理 Lyapunov 函数
├── mini-stability-of-interconnected/   # 复合系统稳定性、分散化控制、无源性/小增益互联分析
├── .gitignore                          # 构建产物、IDE 文件、操作系统文件
├── README.md                           # 本文件（英文）
└── README-CN.md                        # 本文件（中文）
```

## 许可证

MIT
