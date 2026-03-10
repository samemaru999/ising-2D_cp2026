# データ構造定義書

**プロジェクト名**: Ising2D
**最終更新日**: 2026-03-10

---

## 格子

格子状態は `Matrix{Int}`（L×L）で表現する。各要素は {-1, +1} の2値。

---

## パラメータ型

### IsingParams

| フィールド | 型 | 集合 | 説明 |
|-----------|-----|------|------|
| `L` | `Int` | {n ∈ Z \| n ≥ 2} | 格子サイズ |
| `J` | `Float64` | R \ {0} | 交換相互作用定数 |
| `h` | `Float64` | R | 外部磁場 |

### SimulationConfig

| フィールド | 型 | 集合 | 説明 |
|-----------|-----|------|------|
| `n_sweeps` | `Int` | Z⁺ | 測定用MCスイープ数 |
| `n_equilibration` | `Int` | Z≥0 | 熱平衡化スイープ数 |
| `sample_interval` | `Int` | Z⁺ | サンプリング間隔 |
| `seed` | `Union{Int, Nothing}` | Z ∪ {Nothing} | 乱数シード |

---

## 物理量の型

### Observables

1つの格子配置から計算される瞬時物理量。

| フィールド | 型 | 集合 | 説明 |
|-----------|-----|------|------|
| `energy_per_site` | `Float64` | [-2\|J\|, 2\|J\|] | e = E / N |
| `magnetization_per_site` | `Float64` | [-1, 1] | m = M / N |

### ThermalAverages

MCサンプリングの蓄積データ。

| フィールド | 型 | 集合 | 説明 |
|-----------|-----|------|------|
| `sum_e` | `Float64` | R | Σ eᵢ |
| `sum_e2` | `Float64` | R≥0 | Σ eᵢ² |
| `sum_abs_m` | `Float64` | R≥0 | Σ \|mᵢ\| |
| `sum_m2` | `Float64` | R≥0 | Σ mᵢ² |
| `sum_m4` | `Float64` | R≥0 | Σ mᵢ⁴ |
| `n_samples` | `Int` | Z≥0 | サンプル数 |

### ThermodynamicQuantities

ThermalAverages から導出される熱力学量。

| フィールド | 型 | 集合 | 説明 |
|-----------|-----|------|------|
| `mean_energy` | `Float64` | [-2\|J\|, 2\|J\|] | ⟨e⟩ |
| `mean_abs_magnetization` | `Float64` | [0, 1] | ⟨\|m\|⟩ |
| `specific_heat` | `Float64` | R≥0 | C = (N/T²)(⟨e²⟩ - ⟨e⟩²) |
| `susceptibility` | `Float64` | R≥0 | χ = (N/T)(⟨m²⟩ - ⟨\|m\|⟩²) |
| `binder_cumulant` | `Float64` | [0, 2/3] | U_L = 1 - ⟨m⁴⟩/(3⟨m²⟩²) |

---

## 結果型

### SingleTemperatureResult

| フィールド | 型 | 集合 | 説明 |
|-----------|-----|------|------|
| `params` | `IsingParams` | — | モデルパラメータ |
| `temperature` | `Float64` | R⁺ | T = 1/β |
| `averages` | `ThermalAverages` | — | 蓄積データ |
| `thermodynamics` | `ThermodynamicQuantities` | — | 導出熱力学量 |
| `acceptance_rate` | `Float64` | [0, 1] | メトロポリス受容率 |

### TemperatureSweepResult

| フィールド | 型 | 集合 | 説明 |
|-----------|-----|------|------|
| `params` | `IsingParams` | — | 共通パラメータ |
| `config` | `SimulationConfig` | — | シミュレーション設定 |
| `temperatures` | `Vector{Float64}` | (R⁺)ⁿ | 温度配列 |
| `results` | `Vector{SingleTemperatureResult}` | — | 各温度の結果 |

---

## 性能最適化型

### BoltzmannTable

h = 0 の場合、ΔE ∈ {-8J, -4J, 0, 4J, 8J} の5通り。

| フィールド | 型 | 集合 | 説明 |
|-----------|-----|------|------|
| `factors` | `Dict{Int, Float64}` | Z → R⁺ | ΔE/J → exp(-βΔE) |
| `beta` | `Float64` | R⁺ | 逆温度 |
| `J` | `Float64` | R \ {0} | 交換相互作用定数 |

---

## eDSL型

### SimulationCommand（抽象型）

| 具体型 | フィールド | 説明 |
|--------|-----------|------|
| `RunSimulation` | params::IsingParams, temperature::Float64 ∈ R⁺, config::SimulationConfig | 単一温度MC実行 |
| `RunTemperatureSweep` | params::IsingParams, temperatures::Vector{Float64} ⊂ (R⁺)ⁿ, config::SimulationConfig | 温度掃引実行 |
| `SaveResult` | data::Any, filename::String, format::Symbol ∈ {:jld2} | データ保存 |
| `LoadResult` | filename::String, expected_type::Type | データ読み込み |

### ComputationResult{T}（抽象型）

| 具体型 | フィールド | 説明 |
|--------|-----------|------|
| `Success{T}` | value::T | 成功値 |
| `Failure{T}` | error::String, context::Dict{String, Any} | エラー情報 |

---

## 型の依存関係

```
IsingParams ──────────┐
                      ├──▶ SingleTemperatureResult
SimulationConfig ─────┤
                      │
Observables           │
  ↓ (蓄積)            │
ThermalAverages ──────┤
  ↓ (導出)            │
ThermodynamicQuantities
                      │
                      ▼
              TemperatureSweepResult
```

---

## 改訂履歴

| バージョン | 日付 | 変更内容 |
|-----------|------|---------|
| 1.0 | 2025-12-05 | 初版作成 |
| 2.0 | 2025-12-05 | FDD対応版 |
| 3.0 | 2026-03-10 | 2次元イジングモデル用に全面書き換え |
| 3.1 | 2026-03-10 | 型定義（集合）のみに絞り込み |
