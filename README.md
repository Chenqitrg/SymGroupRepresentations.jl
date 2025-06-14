# SymGroupRepresentations.jl

<!-- [![Build Status](https://github.com/Yue-Zhengyuan/SymGroupRepresentations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Yue-Zhengyuan/SymGroupRepresentations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Yue-Zhengyuan/SymGroupRepresentations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Yue-Zhengyuan/SymGroupRepresentations.jl) -->

Compute Clebsch-Gordan coefficients for small symmetric groups. Compatibility / interoperability with [TensorKit.jl](https://github.com/Jutho/TensorKit.jl).

## Conventions

The irreps of a symmetric group $S_n$ are labelled by partitions of $n$, and arranged in "reversed dictionary order". For example, when $n = 4$, there are 5 irreps:

<center>

|   Irrep   |   Partition of 4   | Dimension | Remark                  |
| :-------: | :----------------: | :-------: | :---------------------- |
|   $4_1$   |        $4$         |     1     | Trivial representation  |
| $3_1 1_1$ |      $3 + 1$       |     3     | Standard representation |
|   $2_2$   |      $2 \times 2$       |     2     |                         |
| $2_1 1_2$ | $2 + (1 \times 2)$ |     3     |                         |
|   $1_4$   |   $(1 \times 4)$   |     1     | Sign representation     |

</center>
