# データ構造定義書

**プロジェクト名**: 2次元イジングモデル モンテカルロシミュレーション (Ising2D)
**バージョン**: 3.0
**作成日**: 2025-12-05
**最終更新日**: 2026-03-10
**設計方針**: FDD (Functional Declarative Design)

---

## 目次

1. [FDDに基づく型設計方針](#1-fddに基づく型設計方針)
2. [Domain Layer型定義](#2-domain-layer型定義)
3. [eDSL設計](#3-edsl設計)
4. [データフォーマット (JLD2)](#4-データフォーマット-jld2)
5. [データフロー](#5-データフロー)
6. [バリデーション設計](#6-バリデーション設計)

---

## 1. FDDに基づく型設計方針

### 1.1 型設計の基本原則

FDDでは、型を「本質的複雑さ」を表現する手段として活用し、「偶有的複雑さ」を最小化する。

```
Type Design Principles (FDD)
├── Immutability First: パラメータ・結果は不変データ構造
├── Make Illegal States Unrepresentable: Smart Constructorで不正値を排除
├── Algebraic Data Types: 直和型・直積型で表現
└── 格子は Matrix{Int} をそのまま使用（実用性優先）
```

### 1.2 Pure/Impure分離における型の役割

| 層 | 型の特徴 | 例 |
|----|----------|-----|
| **Domain Layer (Pure)** | 不変、検証済み、ドメイン意味を持つ | `IsingParams`, `Observables` |
| **Logic Layer (Pure)** | 純粋関数の入出力型 | `ThermalAverages`, `ThermodynamicQuantities` |
| **Interpreter Layer (Impure)** | 実行状態、乱数、IO | `SimulationState`, JLD2ファイル |

### 1.3 格子の表現

格子は `Matrix{Int}`（L×L、各要素 ±1）をそのまま使用する。
ラッパー型を導入しない理由:

- Julia の多重ディスパッチにより `Matrix{Int}` のまま専用関数を定義可能
- 格子はシミュレーション中に in-place 更新される（`!` 関数）ため、不変ラッパーは不適
- `random_lattice`, `uniform_lattice` が生成関数として格子の不変条件を保証

---

## 2. Domain Layer型定義

### 2.1 パラメータ構造体（Smart Constructors付き）

```julia
# src/Ising/Types.jl [Pure]

"""
イジングモデルのパラメータ（不変）

# Invariants
- L ≥ 2（格子サイズ）
- J ≠ 0（交換相互作用）
"""
struct IsingParams
    L::Int          # 格子サイズ（L×L）
    J::Float64      # 交換相互作用定数（正: 強磁性）
    h::Float64      # 外部磁場

    function IsingParams(L::Int, J::Float64=1.0, h::Float64=0.0)
        L ≥ 2 || throw(DomainError(L, "L must be ≥ 2"))
        J ≠ 0 || throw(DomainError(J, "J must be non-zero"))
        new(L, J, h)
    end
end
```

```julia
"""
シミュレーション設定（不変）

# Invariants
- n_sweeps > 0
- n_equilibration ≥ 0
- sample_interval ≥ 1
"""
struct SimulationConfig
    n_sweeps::Int           # 測定用MCスイープ数
    n_equilibration::Int    # 熱平衡化スイープ数（測定に含めない）
    sample_interval::Int    # サンプリング間隔（自己相関低減用）
    seed::Union{Int, Nothing}  # 乱数シード（Nothingならランダム）

    function SimulationConfig(;
        n_sweeps::Int=10000,
        n_equilibration::Int=1000,
        sample_interval::Int=1,
        seed::Union{Int, Nothing}=nothing
    )
        n_sweeps > 0 || throw(DomainError(n_sweeps, "n_sweeps must be > 0"))
        n_equilibration ≥ 0 || throw(DomainError(n_equilibration, "n_equilibration must be ≥ 0"))
        sample_interval ≥ 1 || throw(DomainError(sample_interval, "sample_interval must be ≥ 1"))
        new(n_sweeps, n_equilibration, sample_interval, seed)
    end
end
```

### 2.2 物理量の型

```julia
# src/Ising/Types.jl [Pure]

"""
1つの格子配置から計算される瞬時物理量

純粋関数 energy(), magnetization() の戻り値をまとめる。
"""
struct Observables
    energy_per_site::Float64        # e = E / N
    magnetization_per_site::Float64 # m = M / N
end
```

```julia
"""
熱平均の蓄積データ

MCサンプリング中に逐次更新される。
ビンダーキュムラント計算のため ⟨m⁴⟩ も保持する。
"""
struct ThermalAverages
    sum_e::Float64          # Σ e_i
    sum_e2::Float64         # Σ e_i²
    sum_abs_m::Float64      # Σ |m_i|
    sum_m2::Float64         # Σ m_i²
    sum_m4::Float64         # Σ m_i⁴
    n_samples::Int          # サンプル数

    ThermalAverages() = new(0.0, 0.0, 0.0, 0.0, 0.0, 0)
end
```

```julia
"""
熱力学量（ThermalAveragesから導出される不変値）

# フィールド
- specific_heat: C = (N/T²)(⟨e²⟩ - ⟨e⟩²)
- susceptibility: χ = (N/T)(⟨m²⟩ - ⟨|m|⟩²)
- binder_cumulant: U_L = 1 - ⟨m⁴⟩ / (3⟨m²⟩²)
"""
struct ThermodynamicQuantities
    mean_energy::Float64            # ⟨e⟩
    mean_abs_magnetization::Float64 # ⟨|m|⟩
    specific_heat::Float64          # C
    susceptibility::Float64         # χ
    binder_cumulant::Float64        # U_L
end
```

### 2.3 シミュレーション結果型（不変）

```julia
# src/Ising/Types.jl [Pure]

"""
単一温度でのシミュレーション結果（不変）
"""
struct SingleTemperatureResult
    params::IsingParams
    temperature::Float64                # T (= 1/β)
    averages::ThermalAverages
    thermodynamics::ThermodynamicQuantities
    acceptance_rate::Float64            # 受容率
end
```

```julia
"""
温度掃引の結果（不変）

各温度での結果を保持し、相転移解析に使用する。
"""
struct TemperatureSweepResult
    params::IsingParams                         # 共通パラメータ (L, J, h)
    config::SimulationConfig
    temperatures::Vector{Float64}               # 温度配列
    results::Vector{SingleTemperatureResult}     # 各温度の結果
end
```

### 2.4 性能最適化型

```julia
# src/Ising/Types.jl [Pure]

"""
ボルツマン因子のルックアップテーブル

h = 0 の場合、ΔE は {-8J, -4J, 0, 4J, 8J} の5通り。
exp(-β × ΔE) を事前計算して高速化する。
"""
struct BoltzmannTable
    factors::Dict{Int, Float64}   # ΔE/J(整数) → exp(-β × ΔE)
    beta::Float64
    J::Float64
end
```

---

## 3. eDSL設計

### 3.1 計算コマンドの定義

計算処理を「コマンド」として値で表現し、実行を「インタープリター」に委譲する。

```julia
# src/Ising/DSL.jl [Pure]

"""
シミュレーションコマンドの基底型
"""
abstract type SimulationCommand end

"""
単一温度でのMCシミュレーション実行
"""
struct RunSimulation <: SimulationCommand
    params::IsingParams
    temperature::Float64
    config::SimulationConfig
end

"""
温度掃引シミュレーション実行
"""
struct RunTemperatureSweep <: SimulationCommand
    params::IsingParams
    temperatures::Vector{Float64}
    config::SimulationConfig
end

"""
データ保存コマンド
"""
struct SaveResult <: SimulationCommand
    data::Any
    filename::String
    format::Symbol      # :jld2
end

"""
データ読み込みコマンド
"""
struct LoadResult <: SimulationCommand
    filename::String
    expected_type::Type
end
```

### 3.2 結果型（Either風）

```julia
# src/Ising/DSL.jl [Pure]

"""
計算結果を表す型（Either風）
"""
abstract type ComputationResult{T} end

struct Success{T} <: ComputationResult{T}
    value::T
end

struct Failure{T} <: ComputationResult{T}
    error::String
    context::Dict{String, Any}

    Failure{T}(error::String) where T = new{T}(error, Dict{String, Any}())
    Failure{T}(error::String, ctx::Dict{String, Any}) where T = new{T}(error, ctx)
end

# ユーティリティ関数
is_success(::Success) = true
is_success(::Failure) = false

unwrap(r::Success) = r.value
unwrap(r::Failure) = throw(ErrorException(r.error))

unwrap_or(r::Success, _) = r.value
unwrap_or(::Failure, default) = default

map_result(f, r::Success{T}) where T = Success(f(r.value))
map_result(_, r::Failure{T}) where T = r
```

---

## 4. データフォーマット (JLD2)

### 4.1 ファイル命名規則

```
ising2d_L{格子サイズ}_{内容}_{日付}.jld2

内容:
- T{温度}: 単一温度の結果
- sweep_T{T_min}-{T_max}: 温度掃引の結果

例:
- ising2d_L32_T2.27_20260310.jld2
- ising2d_L32_sweep_T1.50-3.50_20260310.jld2
```

### 4.2 単一温度結果ファイル

```julia
Dict(
    "type" => "SingleTemperatureResult",
    "version" => "3.0",
    "L" => Int,
    "J" => Float64,
    "h" => Float64,
    "T" => Float64,
    "n_sweeps" => Int,
    "n_equilibration" => Int,
    "seed" => Union{Int, Nothing},
    "thermodynamics" => Dict(
        "mean_energy" => Float64,
        "mean_abs_magnetization" => Float64,
        "specific_heat" => Float64,
        "susceptibility" => Float64,
        "binder_cumulant" => Float64,
    ),
    "acceptance_rate" => Float64,
    "metadata" => Dict(
        "created_at" => DateTime,
        "julia_version" => String,
    ),
)
```

### 4.3 温度掃引結果ファイル

```julia
Dict(
    "type" => "TemperatureSweepResult",
    "version" => "3.0",
    "L" => Int,
    "J" => Float64,
    "h" => Float64,
    "temperatures" => Vector{Float64},
    "mean_energy" => Vector{Float64},
    "mean_abs_magnetization" => Vector{Float64},
    "specific_heat" => Vector{Float64},
    "susceptibility" => Vector{Float64},
    "binder_cumulant" => Vector{Float64},
    "acceptance_rates" => Vector{Float64},
    "config" => Dict(
        "n_sweeps" => Int,
        "n_equilibration" => Int,
        "sample_interval" => Int,
        "seed" => Union{Int, Nothing},
    ),
    "metadata" => Dict(
        "created_at" => DateTime,
        "julia_version" => String,
    ),
)
```

---

## 5. データフロー

### 5.1 計算パイプライン

```
IsingParams + SimulationConfig + temperature
        │
        ▼
┌─────────────────────────────┐
│  格子生成 [Pure]             │  random_lattice(L) → Matrix{Int}
└─────────────────────────────┘
        │
        ▼
┌─────────────────────────────┐
│  熱平衡化 [Impure]           │  n_equilibration 回の sweep!
│  （測定なし）                │  格子を in-place 更新
└─────────────────────────────┘
        │
        ▼
┌─────────────────────────────┐
│  測定スイープ [Impure]       │  n_sweeps 回のループ:
│                              │    sweep!(lattice, beta)
│  各ステップで:               │    e = energy(lattice) / N
│  - Observables を計算 [Pure] │    m = magnetization(lattice)
│  - ThermalAverages に蓄積    │    ThermalAverages に加算
└─────────────────────────────┘
        │
        ▼
┌─────────────────────────────┐
│  熱力学量導出 [Pure]         │  ThermalAverages → ThermodynamicQuantities
│                              │  C, χ, U_L を計算
└─────────────────────────────┘
        │
        ▼
  SingleTemperatureResult（不変）
```

### 5.2 温度掃引パイプライン

```
temperatures = [T_min, ..., T_max]
        │
        ▼
  各温度 T に対して上記パイプラインを実行
        │
        ▼
  TemperatureSweepResult（不変）
        │
        ├──▶ JLD2保存 [Impure]
        └──▶ 可視化 [Impure]
             - e(T), |m|(T), C(T), χ(T), U_L(T) プロット
```

### 5.3 メモリ使用量

```
格子サイズ別メモリ推定

L = 32:  lattice = 32×32 × 8 bytes = 8 KB
L = 64:  lattice = 64×64 × 8 bytes = 32 KB
L = 128: lattice = 128×128 × 8 bytes = 128 KB

ThermalAverages: 48 bytes（6つのFloat64/Int）
SingleTemperatureResult: < 1 KB
TemperatureSweepResult (100温度点): < 100 KB

→ メモリはボトルネックにならない
```

---

## 6. バリデーション設計

### 6.1 型レベルバリデーション（Smart Constructors）

```julia
# IsingParams の不変条件
L ≥ 2         # 最小格子サイズ
J ≠ 0         # 相互作用が存在

# SimulationConfig の不変条件
n_sweeps > 0          # 測定ステップが存在
n_equilibration ≥ 0   # 非負
sample_interval ≥ 1   # 最低1ステップごと
```

### 6.2 物理量の検証

```julia
# src/Ising/Validation.jl [Pure]

"""
物理量の妥当性検証
"""
function validate_observables(obs::Observables, params::IsingParams)
    errors = String[]

    # エネルギーの範囲: e ∈ [-2J, 2J] (h=0, per site)
    e_min = -2.0 * abs(params.J)
    e_max = 2.0 * abs(params.J)
    (e_min ≤ obs.energy_per_site ≤ e_max) ||
        push!(errors, "energy_per_site out of range: $(obs.energy_per_site)")

    # 磁化の範囲: m ∈ [-1, 1]
    (-1.0 ≤ obs.magnetization_per_site ≤ 1.0) ||
        push!(errors, "magnetization_per_site out of range: $(obs.magnetization_per_site)")

    isempty(errors) ? Success(obs) : Failure{Observables}(join(errors, "; "))
end
```

### 6.3 既知の極限での検証

| 条件 | 期待値 | 検証方法 |
|------|--------|----------|
| 全スピン上向き (h=0) | E = -2JN, M = N | `energy(uniform_lattice(L))` |
| T → ∞ (β → 0) | ⟨\|m\|⟩ → 0 | 高温シミュレーション |
| T → 0 (β → ∞) | ⟨\|m\|⟩ → 1 | 低温シミュレーション |
| T = T_c ≈ 2.269 | U_L の交差 | 複数 L での温度掃引 |
| チェッカーボード | E = +2JN | `energy(checkerboard)` |

---

## 付録

### A. 型の依存関係

```
IsingParams ──────────┐
                      ├──▶ SingleTemperatureResult
SimulationConfig ─────┤
                      │
Observables ──────────┤
  (energy, magnetization から構成)
                      │
ThermalAverages ──────┤
  (Observables を逐次蓄積)
                      │
ThermodynamicQuantities
  (ThermalAverages から Pure に導出)
                      │
                      ▼
              TemperatureSweepResult
              (SingleTemperatureResult の配列)
```

### B. 関連ドキュメント

- [要件定義書](./requirement.md)
- [技術スタック仕様書](./tech_stack.md)

---

## 改訂履歴

| バージョン | 日付 | 変更内容 |
|-----------|------|---------|
| 1.0 | 2025-12-05 | 初版作成（HONetSync用） |
| 2.0 | 2025-12-05 | FDD対応版 |
| 3.0 | 2026-03-10 | 2次元イジングモデル用に全面書き換え |
