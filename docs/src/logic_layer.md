# Logic Layer Documentation

`src/Ising/` ディレクトリは、2次元イジングモデルの物理量計算・格子操作・統計処理を担当する純粋関数群を含む。このレイヤーは **Pure** であり、副作用を持たない。

## ファイル一覧

- `Lattice.jl`: 格子状態の生成
- `Physics.jl`: 物理量（エネルギー・磁化・ΔE）の計算
- `Statistics.jl`: 熱力学量の導出

---

## 1. Lattice.jl

スピン格子の初期状態を生成する関数群。スピンは `+1`（上向き）と `-1`（下向き）の2値。

### 関数

| 関数 | シグネチャ | 説明 |
|------|-----------|------|
| `random_lattice` | `(L::Int) -> Matrix{Int}` | L×L のランダムスピン配置を生成 |
| `uniform_lattice` | `(L::Int; spin::Int=1) -> Matrix{Int}` | L×L の一様スピン配置を生成 |

```julia
lattice = random_lattice(32)       # 32×32 ランダム
lattice = uniform_lattice(16)      # 16×16 全スピン +1
lattice = uniform_lattice(16; spin=-1)  # 全スピン -1
```

---

## 2. Physics.jl

2次元イジングモデルの物理量を計算する純粋関数群。周期境界条件 (PBC) を使用。

### 関数

#### `magnetization(lattice) -> Float64`
1サイトあたりの磁化 $m = \frac{1}{N}\sum_i s_i$ を計算。

- $m \in [-1, 1]$
- $|m| \approx 1$: 秩序相（低温）
- $m \approx 0$: 無秩序相（高温）

#### `energy(lattice; J=1.0) -> Float64`
ハミルトニアン $H = -J \sum_{\langle i,j \rangle} s_i s_j$ を計算。

- 各ボンドは1回だけカウント（右隣・下隣のみ）
- 周期境界条件: `mod1` で実装

#### `delta_energy(lattice, i, j; J=1.0) -> Float64`
サイト $(i,j)$ のスピン反転時のエネルギー変化。

$$\Delta E = 2J \cdot s_{i,j} \cdot \sum_{\text{nn}} s_{\text{nn}}$$

- 4近傍スピンの和を計算
- メトロポリスステップの中核

### 実装の特徴

- `mod1` による周期境界条件
- `energy()` は `Int` で累積し最後に `-J * E` で `Float64` に変換
- `delta_energy()` はアロケーションゼロ（ベンチマーク: ~4.4 ns）

---

## 3. Statistics.jl

MCサンプリングの蓄積から熱力学量を導出する純粋関数。

### 関数

#### `compute_thermodynamics(avg, N, T) -> ThermodynamicQuantities`

`ThermalAverages` の累積値から以下を計算:

| 物理量 | 計算式 |
|--------|--------|
| ⟨e⟩ | `sum_e / n_samples` |
| ⟨\|m\|⟩ | `sum_abs_m / n_samples` |
| 比熱 C | $\frac{N}{T^2}(\langle e^2 \rangle - \langle e \rangle^2)$ |
| 帯磁率 χ | $\frac{N}{T}(\langle m^2 \rangle - \langle |m| \rangle^2)$ |
| ビンダーキュムラント $U_L$ | $1 - \frac{\langle m^4 \rangle}{3\langle m^2 \rangle^2}$ |

- 引数バリデーション: `n_samples > 0`, `N > 0`, `T > 0`
