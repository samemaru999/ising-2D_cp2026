# Ising2D プロジェクト解析と API リファレンス

## 1. プロジェクト概要

**Ising2D** は、2次元正方格子イジングモデルのモンテカルロシミュレーションを行うJuliaパッケージ。メトロポリス法による単一温度シミュレーションと温度掃引を実装し、相転移（$T_c \approx 2.269$）を観測する。

**Functional Declarative Design (FDD)** アーキテクチャを採用し、純粋な物理計算（Ising/）と副作用を伴う処理（Runtime/）を分離している。

## 2. アーキテクチャ

```
┌─────────────────────────────────────────────────────┐
│  scripts/                                           │
│  quick_start.jl, temperature_sweep.jl,              │
│  run_simulation.jl, benchmark.jl, plot_results.jl   │
├─────────────────────────────────────────────────────┤
│                                                     │
│  ┌───────────────────────────────────────────────┐  │
│  │  Runtime Layer [Impure]                       │  │
│  │  MonteCarlo.jl  - MC更新、温度掃引            │  │
│  │  IO.jl          - JLD2保存・読み込み           │  │
│  │  Visualization.jl - GLMakieプロット            │  │
│  └───────────────────────────────────────────────┘  │
│                      ↑ 呼び出し                      │
│  ┌───────────────────────────────────────────────┐  │
│  │  Ising Layer [Pure]                           │  │
│  │  Lattice.jl    - 格子生成                      │  │
│  │  Physics.jl    - エネルギー、磁化、ΔE          │  │
│  │  Statistics.jl - 熱力学量の導出                 │  │
│  └───────────────────────────────────────────────┘  │
│                      ↑ 型参照                        │
│  ┌───────────────────────────────────────────────┐  │
│  │  Domain Layer [Pure]                          │  │
│  │  Types.jl - 全型定義                           │  │
│  └───────────────────────────────────────────────┘  │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### レイヤー構成

| レイヤー | ファイル | Pure/Impure | 役割 |
|---------|---------|------------|------|
| Domain | `Types.jl` | Pure | 型定義、バリデーション |
| Ising | `Lattice.jl`, `Physics.jl`, `Statistics.jl` | Pure | 格子生成、物理量計算、統計処理 |
| Runtime | `MonteCarlo.jl`, `IO.jl`, `Visualization.jl` | Impure | MC更新（乱数）、ファイルIO、プロット |

## 3. 主要データ構造

### パラメータ

```julia
params = IsingParams(32, 1.0)       # L=32, J=1.0
config = SimulationConfig(5000, 2000, 5)  # sweeps, equil, interval
```

### 結果

```julia
# 単一温度
SingleTemperatureResult(params, T, averages, thermodynamics, acceptance_rate)

# 温度掃引
TemperatureSweepResult(params, config, temperatures, results)
```

## 4. API リファレンス

### Pure関数 (Ising/)

```julia
# 格子生成
lattice = random_lattice(32)
lattice = uniform_lattice(32; spin=-1)

# 物理量
m = magnetization(lattice)              # ⟨s⟩ ∈ [-1, 1]
E = energy(lattice; J=1.0)              # H = -J Σ s_i s_j
dE = delta_energy(lattice, i, j; J=1.0) # スピン反転時の ΔE

# 統計
thermo = compute_thermodynamics(avg, N, T)  # ThermalAverages → ThermodynamicQuantities
```

### Impure関数 (Runtime/)

```julia
# MC更新
accepted = metropolis_step!(lattice, beta; J=1.0)  # 1ステップ
n_accept = sweep!(lattice, beta; J=1.0)            # 1スイープ (L²回)

# 高レベルシミュレーション
avg, total = run_single_temperature(lattice, beta, N, config; J=1.0)
result = run_temperature_sweep(params, config, temperatures)

# IO
save_result("file.jld2", result)
data = load_result("file.jld2")

# 可視化
fig = plot_thermodynamics(sweep_result)
fig = plot_spin_lattice(lattice)
fig = plot_binder_cumulant([sweep1, sweep2])
fig = plot_timeseries(energies, magnetizations)
```

## 5. 使用方法

### 単一温度シミュレーション

```julia
using Ising2D

lattice = random_lattice(32)
beta = 1.0 / 2.269
config = SimulationConfig(5000, 1000, 5)
avg, accepted = run_single_temperature(lattice, beta, 1024, config)
thermo = compute_thermodynamics(avg, 1024, 2.269)
```

### 温度掃引

```julia
using Ising2D

params = IsingParams(32, 1.0)
config = SimulationConfig(5000, 2000, 5)
temps = collect(1.5:0.05:3.5)
result = run_temperature_sweep(params, config, temps)
fig = plot_thermodynamics(result)
```

### CLIスクリプト

```bash
# デフォルト実行
julia --project scripts/run_simulation.jl

# カスタムパラメータ
julia --project scripts/run_simulation.jl --L 64 --Tmin 2.0 --Tmax 2.6 --dT 0.02

# ベンチマーク
julia --project scripts/benchmark.jl
```

## 6. テスト

```bash
# 全テスト
julia --project -e 'using Pkg; Pkg.test()'

# 個別テスト
julia --project -e 'using Test, Ising2D; @testset "Statistics" begin include("test/Unit/test_statistics.jl") end'
```

テスト構成:
- Unit tests (80%): `test/Unit/test_*.jl`
- Integration / E2E (20%): `test/Ising/`, `test/E2E/`

## 7. パフォーマンス

BenchmarkTools による計測結果 (L=32, Apple Silicon):

| 関数 | 時間 | アロケーション |
|------|------|-------------|
| `delta_energy` | ~4.4 ns | 0 |
| `metropolis_step!` | ~17 ns | 0 |
| `sweep!` | ~16 μs | 0 |
| 1スピン更新 | ~15 ns | 0 |
