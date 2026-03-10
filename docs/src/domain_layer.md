# Domain Layer Documentation

`src/Types.jl` は、FDD (Functional Declarative Design) におけるドメインモデル（型定義）を担当する。このレイヤーは純粋 (Pure) であり、副作用を持たない。

## ファイル

- `src/Types.jl`: 全型定義

---

## 型一覧

### パラメータ型

#### `IsingParams`
イジングモデルの物理パラメータ。

```julia
struct IsingParams
    L::Int        # 格子サイズ (L ≥ 2)
    J::Float64    # 結合定数 (≠ 0)
end
```

- 内部コンストラクタでバリデーション: `L >= 2`, `J != 0.0`

#### `SimulationConfig`
モンテカルロシミュレーションの設定。

```julia
struct SimulationConfig
    n_sweeps::Int                   # 測定スイープ数 (> 0)
    n_equilibration::Int            # 平衡化スイープ数 (≥ 0)
    sample_interval::Int            # サンプリング間隔 (> 0)
    seed::Union{Int, Nothing}       # 乱数シード (optional)
end
```

### 物理量の型

#### `Observables`
1サンプルの観測量。

```julia
struct Observables
    energy_per_site::Float64
    magnetization_per_site::Float64
end
```

#### `ThermalAverages`
MCサンプリングの累積統計量。**mutable** — 測定ループ内で逐次更新される。

```julia
mutable struct ThermalAverages
    sum_e::Float64      # Σ e_i
    sum_e2::Float64     # Σ e_i²
    sum_abs_m::Float64  # Σ |m_i|
    sum_m2::Float64     # Σ m_i²
    sum_m4::Float64     # Σ m_i⁴
    n_samples::Int      # サンプル数
end
```

- `ThermalAverages()` でゼロ初期化

#### `ThermodynamicQuantities`
熱力学量の最終結果（immutable）。

```julia
struct ThermodynamicQuantities
    mean_energy::Float64              # ⟨e⟩
    mean_abs_magnetization::Float64   # ⟨|m|⟩
    specific_heat::Float64            # C = N/T² (⟨e²⟩ - ⟨e⟩²)
    susceptibility::Float64           # χ = N/T (⟨m²⟩ - ⟨|m|⟩²)
    binder_cumulant::Float64          # U_L = 1 - ⟨m⁴⟩/(3⟨m²⟩²)
end
```

### 結果型

#### `SingleTemperatureResult`
単一温度シミュレーションの結果。

```julia
struct SingleTemperatureResult
    params::IsingParams
    temperature::Float64
    averages::ThermalAverages
    thermodynamics::ThermodynamicQuantities
    acceptance_rate::Float64
end
```

#### `TemperatureSweepResult`
温度掃引シミュレーションの結果。

```julia
struct TemperatureSweepResult
    params::IsingParams
    config::SimulationConfig
    temperatures::Vector{Float64}
    results::Vector{SingleTemperatureResult}
end
```

### 性能最適化型

#### `BoltzmannTable`
ΔE → exp(-βΔE) の事前計算テーブル（未使用、将来の最適化用）。

```julia
struct BoltzmannTable
    factors::Dict{Int, Float64}
    beta::Float64
    J::Float64
end
```

---

## 設計原則

- **不正な状態を型で排除**: 内部コンストラクタでバリデーション (`L >= 2`, `n_sweeps > 0` 等)
- **Immutableが基本**: `ThermalAverages` のみ mutable（測定ループの効率のため）
- **具体的な型**: `Union` は `seed` フィールドの `Union{Int, Nothing}` のみに限定
