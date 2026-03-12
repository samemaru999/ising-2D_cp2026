# ディレクトリ構造ドキュメント

**プロジェクト名**: 2次元イジングモデル モンテカルロシミュレーション (Ising2D)
**最終更新日**: 2026-03-11

---

## プロジェクト全体構造

```
ising-2D_cp2026/
├── .claude/
│   ├── CLAUDE.md              # Claude Code プロジェクト設定 (FDD ワークフロー定義)
│   ├── settings.json
│   └── settings.local.json
├── .vscode/
│   └── settings.json
├── data/                      # シミュレーション入出力データ
│   └── output/
│       ├── <timestamp>/       # タイムスタンプ別出力ディレクトリ (例: 20260310_223737/)
│       │   └── *.jld2         # 各温度・温度掃引の保存結果
│       └── figs/              # 可視化出力画像
│           ├── m2_vs_T.png    # ⟨M²⟩ vs T プロット
│           └── thermodynamics_*.png
├── docs/                      # ドキュメント
│   ├── directory_structure.md # 本ファイル
│   ├── implementation_task_checklist.md
│   ├── spec/                  # 仕様ドキュメント
│   │   ├── requirement.md
│   │   ├── system_architecture.md
│   │   ├── tech_stack.md
│   │   └── data_structure.md
│   └── src/                   # ソースコード層別ドキュメント
│       ├── domain_layer.md
│       ├── fdd_implementation_guide.md
│       ├── interpreters_layer.md
│       ├── io_layer.md
│       ├── logic_layer.md
│       └── project_analysis.md
├── scripts/                   # 実行スクリプト群
│   ├── benchmark.jl           # ホットパスのベンチマーク計測
│   ├── plot_m2_vs_T.jl        # ⟨M²⟩ vs T プロット (quick_mc 使用)
│   ├── plot_results.jl        # 保存済みJLD2ファイルの可視化
│   ├── quick_start.jl         # 単一温度シミュレーション (入門用)
│   ├── run_simulation.jl      # 温度掃引 CLI (引数対応フルバージョン)
│   └── temperature_sweep.jl   # 温度掃引 (固定パラメータ版)
├── src/                       # Julia ソースコード
│   ├── Ising2D.jl             # メインモジュール (include + export)
│   ├── Types.jl               # 全型定義
│   ├── Ising/                 # Pure 層 (副作用なし)
│   │   ├── Lattice.jl         # 格子生成関数
│   │   ├── Physics.jl         # 物理量計算関数
│   │   └── Statistics.jl      # 熱力学量導出関数
│   └── Runtime/               # Impure 層 (副作用あり)
│       ├── MonteCarlo.jl      # メトロポリス法・MCループ・quick_mc
│       ├── IO.jl              # JLD2保存・読み込み
│       └── Visualization.jl   # GLMakie プロット関数
├── test/                      # テストスイート
│   ├── runtests.jl            # テストエントリーポイント
│   ├── Unit/                  # ユニットテスト (全体の80%)
│   │   ├── test_types.jl
│   │   ├── test_observables.jl
│   │   ├── test_results.jl
│   │   ├── test_lattice.jl
│   │   ├── test_physics.jl
│   │   ├── test_montecarlo.jl
│   │   └── test_statistics.jl
│   ├── Ising/                 # 統合テスト
│   │   └── test_integration.jl
│   └── E2E/                   # エンドツーエンドテスト
│       └── test_simulation.jl
├── .gitignore
├── Manifest.toml              # 依存関係のロックファイル
├── mise.toml                  # mise (ツールバージョン管理)
├── Project.toml               # Julia プロジェクト定義・依存宣言
└── README.md
```

---

## src/ ディレクトリの詳細

Pure/Impure 分離 (FDD アーキテクチャ) に基づいて設計されている。

### src/Types.jl - 全型定義

プロジェクト全体で使用するデータ型を一元管理する。

| 型名 | 種別 | 説明 |
|------|------|------|
| `IsingParams` | struct | 格子サイズ `L`、交換相互作用 `J` |
| `SimulationConfig` | struct | MC スイープ数、平衡化回数、サンプリング間隔、乱数シード |
| `Observables` | struct | 1 サンプルあたりのエネルギー・磁化 |
| `ThermalAverages` | mutable struct | MCサンプルの累積和（sum_e, sum_e2, sum_abs_m, sum_m2, sum_m4, n_samples） |
| `ThermodynamicQuantities` | struct | 熱平均から導出した巨視的物理量（⟨e⟩, ⟨\|m\|⟩, C, χ, U_L） |
| `SingleTemperatureResult` | struct | 単一温度シミュレーションの全結果 |
| `TemperatureSweepResult` | struct | 温度掃引全結果のコレクション |
| `BoltzmannTable` | struct | ボルツマン因子のルックアップテーブル（性能最適化用） |

### src/Ising/ - Pure 層

副作用のない純粋関数群。単体テストが容易で、決定論的な計算を保証する。

| ファイル | 主要関数 | 説明 |
|---------|---------|------|
| `Lattice.jl` | `random_lattice(L)` | L×L ランダムスピン配置を生成 |
| `Lattice.jl` | `uniform_lattice(L; spin=1)` | L×L 一様スピン配置を生成 |
| `Physics.jl` | `magnetization(lattice)` | 1サイトあたりの磁化 M = (1/N) Σ s_i |
| `Physics.jl` | `energy(lattice; J)` | 全エネルギー H = -J Σ s_i s_j (周期境界条件) |
| `Physics.jl` | `delta_energy(lattice, i, j; J)` | スピン反転時のエネルギー変化 ΔE |
| `Statistics.jl` | `compute_thermodynamics(avg, N, T)` | ThermalAverages から C, χ, U_L を導出 |

### src/Runtime/ - Impure 層

乱数・IO・外部ライブラリへの副作用を持つ関数群。

| ファイル | 主要関数 | 説明 |
|---------|---------|------|
| `MonteCarlo.jl` | `metropolis_step!(lattice, beta; J)` | 1サイトのメトロポリス更新 |
| `MonteCarlo.jl` | `sweep!(lattice, beta; J)` | 1モンテカルロスイープ (L² 回の更新試行) |
| `MonteCarlo.jl` | `run_single_temperature(...)` | 平衡化 + 測定ループ |
| `MonteCarlo.jl` | `run_temperature_sweep(params, config, temperatures)` | 温度配列全体をループする掃引 |
| `MonteCarlo.jl` | `quick_mc(T, L; ...)` | 簡易 MC インターフェース (→ Dict{String, Float64}) |
| `IO.jl` | `save_result(filename, result)` | JLD2 形式で結果を保存 |
| `IO.jl` | `load_result(filename)` | JLD2 ファイルを Dict として読み込み |
| `Visualization.jl` | `plot_spin_lattice(lattice)` | スピン配置のヒートマップ |
| `Visualization.jl` | `plot_thermodynamics(result)` | ⟨e⟩, ⟨\|m\|⟩, C, χ の 4 パネルプロット |
| `Visualization.jl` | `plot_binder_cumulant(results)` | 複数 L のビンダーキュムラント重ね描き |
| `Visualization.jl` | `plot_timeseries(energies, magnetizations)` | MC 時系列プロット |

---

## test/ ディレクトリの詳細

TDD 採用。ユニットテストが全体の約 80% を占める。

| ファイル | 対応ソース | テスト内容 |
|---------|-----------|-----------|
| `Unit/test_types.jl` | `Types.jl` | IsingParams, SimulationConfig のバリデーション |
| `Unit/test_observables.jl` | `Types.jl` | Observables, ThermalAverages, ThermodynamicQuantities の構築 |
| `Unit/test_results.jl` | `Types.jl` | BoltzmannTable, SingleTemperatureResult, TemperatureSweepResult の構築 |
| `Unit/test_lattice.jl` | `Ising/Lattice.jl` | 格子生成・サイズ・スピン値の検証 |
| `Unit/test_physics.jl` | `Ising/Physics.jl` | 磁化・エネルギー・ΔE の物理的正確性 |
| `Unit/test_montecarlo.jl` | `Runtime/MonteCarlo.jl` | metropolis_step!, sweep!, quick_mc のふるまい |
| `Unit/test_statistics.jl` | `Ising/Statistics.jl` | compute_thermodynamics の数値検証 |
| `Ising/test_integration.jl` | 複数ソース | MC ワークフロー全体の統合テスト |
| `E2E/test_simulation.jl` | 複数ソース | 高温/低温の磁化物理特性の E2E 検証 |

---

## scripts/ ディレクトリの詳細

| スクリプト | 目的 | 主な出力 |
|-----------|------|---------|
| `quick_start.jl` | 単一温度 MC の入門的実行 | `data/output/<timestamp>/L<L>_T<T>.jld2` |
| `temperature_sweep.jl` | 固定パラメータでの温度掃引 | `data/output/<timestamp>/sweep_L<L>.jld2` |
| `run_simulation.jl` | CLI 引数対応の温度掃引 | `data/output/<timestamp>/sweep_L<L>.jld2` |
| `plot_results.jl` | 保存済み JLD2 の可視化 | `data/output/figs/*.png` |
| `plot_m2_vs_T.jl` | `quick_mc` で ⟨M²⟩ vs T を作図 | `data/output/figs/m2_vs_T.png` |
| `benchmark.jl` | ホットパスの BenchmarkTools 計測 | 標準出力 (delta_energy, metropolis_step!, sweep!) |

---

## data/ ディレクトリの詳細

```
data/
└── output/
    ├── <yyyymmdd_HHMMSS>/     # 実行タイムスタンプ別ディレクトリ
    │   ├── L<L>_T<T>.jld2    # SingleTemperatureResult (quick_start.jl)
    │   └── sweep_L<L>.jld2   # TemperatureSweepResult (temperature_sweep.jl 等)
    └── figs/
        ├── m2_vs_T.png
        ├── thermodynamics_L<L>_<timestamp>.png
        └── single_L<L>_T<T>_<timestamp>.png
```

- JLD2 ファイルには `IsingParams`, `SimulationConfig` を含む全パラメータが保存される
- `data/output/figs/*.png` は git で追跡される（`.gitignore` にルール設定済み）

---

## アーキテクチャの原則

### Pure / Impure 分離

```
src/
├── Types.jl        # Pure (型定義)
├── Ising/          # Pure (副作用なし、決定論的)
│   ├── Lattice.jl
│   ├── Physics.jl
│   └── Statistics.jl
└── Runtime/        # Impure (乱数・IO・外部ライブラリ)
    ├── MonteCarlo.jl
    ├── IO.jl
    └── Visualization.jl
```

### 依存関係の方向

```
Impure (Runtime/) → Pure (Ising/, Types.jl)
```

Pure 層は Impure 層に依存しない。これにより Pure 関数の単体テストが副作用なしで実行できる。

---

## 関連ドキュメント

- [仕様: 要件定義](./spec/requirement.md)
- [仕様: システムアーキテクチャ](./spec/system_architecture.md)
- [仕様: 技術スタック](./spec/tech_stack.md)
- [仕様: データ構造](./spec/data_structure.md)
- [実装チェックリスト](./implementation_task_checklist.md)
- [FDD 実装ガイド](./src/fdd_implementation_guide.md)
