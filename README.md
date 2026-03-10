# Ising2D - 2次元イジングモデル モンテカルロシミュレーション

2次元正方格子上のイジングモデルのモンテカルロシミュレーションパッケージ

## Overview

本パッケージは、2次元正方格子上のイジングモデルの相転移・臨界現象を数値的に解析するためのシミュレーションフレームワークです。**メトロポリス法**を用いたモンテカルロシミュレーションを行い、各種物理量を計算します。

### Key Features

- **FDD (Functional Declarative Design)** アーキテクチャによる Pure/Impure 分離
- **メトロポリス法** によるモンテカルロシミュレーション
- **物理量の計算**: 磁化、エネルギー、比熱、磁化率
- **有限サイズスケーリング** による臨界指数の推定
- **相転移の解析**: 臨界温度 T_c = 2J/(k_B ln(1+√2)) 付近の振る舞い

## Installation

### Requirements

- Julia 1.10 以上
- 依存パッケージ: StaticArrays, JLD2, LinearAlgebra, Statistics

### Setup

```bash
cd /path/to/ising-2D_cp2026
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Model

### ハミルトニアン

```
H = -J Σ_{<i,j>} s_i s_j - h Σ_i s_i
```

- `s_i ∈ {+1, -1}`: 各格子点のスピン
- `J`: 交換相互作用定数 (J > 0: 強磁性)
- `h`: 外部磁場
- `<i,j>`: 最近接格子点対

### メトロポリス法

1. ランダムにスピンを1つ選択
2. スピン反転によるエネルギー変化 ΔE を計算
3. ΔE ≤ 0 なら反転を受理、ΔE > 0 なら確率 exp(-ΔE/k_BT) で受理
4. 1-3を繰り返す

### 臨界温度 (Onsager厳密解)

```
T_c = 2J / (k_B ln(1 + √2)) ≈ 2.269 J/k_B
```

## License

MIT
