# IO & Visualization Layer Documentation

`src/Runtime/IO.jl` と `src/Runtime/Visualization.jl` は、データ永続化と可視化を担当する **Impure** レイヤー。

## ファイル

- `src/Runtime/IO.jl`: JLD2形式での保存・読み込み
- `src/Runtime/Visualization.jl`: GLMakieによるプロット

---

## 1. IO.jl

`JLD2.jl` パッケージを使用してシミュレーション結果を保存・読み込みする。

### 保存関数

#### `save_result(filename, result::TemperatureSweepResult)`
温度掃引結果を保存。保存されるフィールド:
- `params`, `config`, `temperatures`, `results`

#### `save_result(filename, result::SingleTemperatureResult)`
単一温度結果を保存。保存されるフィールド:
- `params`, `temperature`, `averages`, `thermodynamics`, `acceptance_rate`

### 読み込み関数

#### `load_result(filename) -> Dict{String, Any}`
JLD2ファイルからデータをDictとして読み込む。

### ファイル命名規則

| パターン | 例 |
|---------|-----|
| 単一温度 | `L32_T2.269.jld2` |
| 温度掃引 | `sweep_L32.jld2` |

出力先: `data/output/{timestamp}/`

---

## 2. Visualization.jl

`GLMakie.jl` を使用したプロット関数群。全関数は `Figure` オブジェクトを返す。

### 関数一覧

#### `plot_spin_lattice(lattice) -> Figure`
スピン配置のヒートマップ。青 (-1) / 赤 (+1)。

#### `plot_thermodynamics(result::TemperatureSweepResult) -> Figure`
物理量の温度依存性を4パネルで表示:
- ⟨e⟩(T), ⟨|m|⟩(T), C(T), χ(T)

#### `plot_binder_cumulant(results::Vector{TemperatureSweepResult}) -> Figure`
複数格子サイズの $U_L(T)$ を重ね描き。交点が臨界温度 $T_c$ を示す。

#### `plot_timeseries(energies, magnetizations) -> Figure`
エネルギー・磁化のモンテカルロ時系列（平衡化の確認用）。

---

## スクリプト

### `scripts/quick_start.jl`
単一温度 (T=2.269) でのMCシミュレーション。

### `scripts/temperature_sweep.jl`
T=1.5〜3.5 を 0.05 刻みで掃引（パラメータ固定版）。

### `scripts/run_simulation.jl`
コマンドライン引数対応の温度掃引スクリプト。

```bash
julia --project scripts/run_simulation.jl --L 64 --Tmin 2.0 --Tmax 2.6 --dT 0.02
```

引数: `--L`, `--J`, `--Tmin`, `--Tmax`, `--dT`, `--n_sweeps`, `--n_eq`, `--sample_interval`, `--no-plot`

### `scripts/plot_results.jl`
保存済みJLD2ファイルを読み込んでプロット。

```bash
julia --project scripts/plot_results.jl data/output/20260310/sweep_L32.jld2
```

### `scripts/benchmark.jl`
ホットパス (`delta_energy`, `metropolis_step!`, `sweep!`) のBenchmarkTools計測。
