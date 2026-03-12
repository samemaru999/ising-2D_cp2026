# 実装タスクチェックリスト

**プロジェクト名**: 2次元イジングモデル モンテカルロシミュレーション (Ising2D)
**バージョン**: 0.0.0
**作成日**: 2026-03-11
**最終更新日**: 2026-03-11
**設計手法**: Functional Declarative Design (FDD)

---

## 目次

1. [フェーズ概要](#フェーズ概要)
2. [FDD Phase 1: 要件分析](#fdd-phase-1-要件分析)
3. [FDD Phase 2: アーキテクチャ設計](#fdd-phase-2-アーキテクチャ設計)
4. [FDD Phase 3: ドメインモデル設計](#fdd-phase-3-ドメインモデル設計)
5. [FDD Phase 4: インターフェース設計](#fdd-phase-4-インターフェース設計)
6. [FDD Phase 5: 実装](#fdd-phase-5-実装)
7. [FDD Phase 6: テスト](#fdd-phase-6-テスト)
8. [機能要件別タスク](#機能要件別タスク)
9. [非機能要件タスク](#非機能要件タスク)
10. [優先度と依存関係](#優先度と依存関係)

---

## フェーズ概要

### FDDフェーズと進捗

| FDD Phase | 名称 | 状態 |
|-----------|------|------|
| Phase 1 | 要件分析 | 完了 |
| Phase 2 | アーキテクチャ設計 | 完了 |
| Phase 3 | ドメインモデル設計 | 完了 |
| Phase 4 | インターフェース設計 | 完了 |
| Phase 5 | 実装 | 完了（コア機能） |
| Phase 6 | テスト | 完了（コア機能） |

### 機能要件との対応

| 機能要件ID | 機能名 | 実装状態 | 優先度 |
|-----------|--------|---------|--------|
| F-LAT | 格子生成 | 完了 | 最高 |
| F-PHY | 物理量計算 | 完了 | 最高 |
| F-MC | メトロポリス法 MC | 完了 | 最高 |
| F-STAT | 熱力学量導出 | 完了 | 高 |
| F-IO | JLD2 保存・読み込み | 完了 | 高 |
| F-VIS | 可視化 | 完了 | 中 |
| F-BENCH | ベンチマーク | 完了 | 低 |

---

## FDD Phase 1: 要件分析

### 1.1 ドメイン理解

- ✅ **TASK-1.1.1**: イジングモデルのドメイン概念の整理
  - ハミルトニアン H = -J Σ_{<i,j>} s_i s_j の定義
  - 周期境界条件の採用
  - スピン変数: s_i ∈ {+1, -1}

- ✅ **TASK-1.1.2**: 観測物理量の整理
  - 磁化 M = (1/N) Σ s_i
  - エネルギー E = -J Σ s_i s_j
  - 比熱 C = (N/T²)(⟨E²⟩ - ⟨E⟩²)
  - 帯磁率 χ = (N/T)(⟨M²⟩ - ⟨|M|⟩²)
  - ビンダーキュムラント U_L = 1 - ⟨M⁴⟩ / (3⟨M²⟩²)
  - 正確な臨界温度 Tc = 2/ln(1 + √2) ≈ 2.269

### 1.2 ユーザーシナリオ作成

- ✅ **TASK-1.2.1**: 基本シナリオの確認
  - 単一温度シミュレーション (`quick_start.jl`)
  - 温度掃引シミュレーション (`temperature_sweep.jl`, `run_simulation.jl`)
  - 簡易インターフェース (`quick_mc`)

---

## FDD Phase 2: アーキテクチャ設計

### 2.1 レイヤー構造設計

- ✅ **TASK-2.1.1**: Pure / Impure 分離アーキテクチャの確定
  - Pure 層: `src/Types.jl`, `src/Ising/`
  - Impure 層: `src/Runtime/`
  - 依存方向: Impure → Pure（逆方向なし）

- ✅ **TASK-2.1.2**: ディレクトリ構造の確定
  - `src/Ising/` - 副作用なし純粋関数（格子・物理・統計）
  - `src/Runtime/` - 副作用あり（MC・IO・可視化）

---

## FDD Phase 3: ドメインモデル設計

### 3.1 型定義 (src/Types.jl - Pure)

- ✅ **TASK-3.1.1**: パラメータ型の実装
  - `IsingParams` (L, J) - バリデーション付きコンストラクタ
  - `SimulationConfig` (n_sweeps, n_equilibration, sample_interval, seed)

- ✅ **TASK-3.1.2**: 物理量・蓄積型の実装
  - `Observables` - 1 サンプルの (energy_per_site, magnetization_per_site)
  - `ThermalAverages` (mutable) - sum_e, sum_e2, sum_abs_m, sum_m2, sum_m4, n_samples

- ✅ **TASK-3.1.3**: 結果型の実装
  - `ThermodynamicQuantities` - 熱力学量 (mean_energy, mean_abs_magnetization, specific_heat, susceptibility, binder_cumulant)
  - `SingleTemperatureResult` - 単一温度結果
  - `TemperatureSweepResult` - 温度掃引結果

- ✅ **TASK-3.1.4**: 性能最適化型の実装
  - `BoltzmannTable` - ボルツマン因子ルックアップテーブル (将来の最適化用)

---

## FDD Phase 4: インターフェース設計

### 4.1 モジュールインターフェース

- ✅ **TASK-4.1.1**: エクスポートインターフェースの定義 (`src/Ising2D.jl`)
  - 型: IsingParams, SimulationConfig, Observables, ThermalAverages, ThermodynamicQuantities, SingleTemperatureResult, TemperatureSweepResult, BoltzmannTable
  - Pure 関数: random_lattice, uniform_lattice, magnetization, energy, delta_energy, compute_thermodynamics
  - Impure 関数: metropolis_step!, sweep!, run_single_temperature, run_temperature_sweep, quick_mc, save_result, load_result
  - 可視化: plot_spin_lattice, plot_thermodynamics, plot_binder_cumulant, plot_timeseries

---

## FDD Phase 5: 実装

### 5.1 基盤実装

- ✅ **TASK-5.1.1**: Project.toml の作成
  - 依存: JLD2, GLMakie, BenchmarkTools, Statistics 等
  - `using Pkg; Pkg.instantiate()` で依存関係がインストールされる

- ✅ **TASK-5.1.2**: ディレクトリ構造の作成
  - `src/{Ising,Runtime}/` が存在する

- ✅ **TASK-5.1.3**: メインモジュール `src/Ising2D.jl` の作成
  - `using Ising2D` でモジュールが読み込める

### 5.2 Pure 層実装 (src/Ising/)

#### 5.2.1 格子生成 (Lattice.jl)

- ✅ **TASK-5.2.1**: `random_lattice(L::Int) -> Matrix{Int}`
  - L×L のランダムスピン配置（+1 / -1）を生成

- ✅ **TASK-5.2.2**: `uniform_lattice(L::Int; spin::Int=1) -> Matrix{Int}`
  - L×L の一様スピン配置を生成
  - spin ∉ {+1, -1} で AssertionError

#### 5.2.2 物理量計算 (Physics.jl)

- ✅ **TASK-5.2.3**: `magnetization(lattice::Matrix{Int}) -> Float64`
  - M = (1/N) Σ s_i を計算

- ✅ **TASK-5.2.4**: `energy(lattice::Matrix{Int}; J::Float64=1.0) -> Float64`
  - H = -J Σ_{<i,j>} s_i s_j を周期境界条件で計算
  - 各ボンドを 1 回のみカウント（右隣・下隣）

- ✅ **TASK-5.2.5**: `delta_energy(lattice::Matrix{Int}, i::Int, j::Int; J::Float64=1.0) -> Float64`
  - ΔE = 2J × s(i,j) × (4近傍スピンの和) を計算
  - メトロポリス更新の高速化に使用

#### 5.2.3 統計処理 (Statistics.jl)

- ✅ **TASK-5.2.6**: `compute_thermodynamics(avg::ThermalAverages, N::Int, T::Float64) -> ThermodynamicQuantities`
  - 比熱: C = (N/T²)(⟨E²⟩ - ⟨E⟩²)
  - 帯磁率: χ = (N/T)(⟨M²⟩ - ⟨|M|⟩²)
  - ビンダーキュムラント: U_L = 1 - ⟨M⁴⟩ / (3⟨M²⟩²)
  - n_samples=0, N=0, T≤0 で ArgumentError

### 5.3 Impure 層実装 (src/Runtime/)

#### 5.3.1 モンテカルロ法 (MonteCarlo.jl)

- ✅ **TASK-5.3.1**: `metropolis_step!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0) -> Bool`
  - ランダムにサイト (i,j) を選択
  - ΔE ≤ 0 または exp(-βΔE) > rand() で反転を受容
  - 受容した場合 true を返す

- ✅ **TASK-5.3.2**: `sweep!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0) -> Int`
  - L×L 回の metropolis_step! を実行（1 モンテカルロスイープ）
  - 受容回数を返す

- ✅ **TASK-5.3.3**: `run_single_temperature(lattice, beta, N, config; J) -> (ThermalAverages, Int)`
  - n_equilibration 回の平衡化スイープ
  - n_sweeps 回の測定スイープ（sample_interval ごとにサンプリング）
  - (ThermalAverages, total_accepted) を返す

- ✅ **TASK-5.3.4**: `run_temperature_sweep(params, config, temperatures) -> TemperatureSweepResult`
  - 温度配列を高温から順にループ
  - 前の温度の格子状態を次の温度に引き継ぐ（高温→低温の順に初期化）
  - 最終的に温度昇順で結果を返す

- ✅ **TASK-5.3.5**: `quick_mc(T::Float64, L::Int; ...) -> Dict{String, Float64}`
  - T > 0, L ≥ 2 のバリデーション付き
  - デフォルト: n_sweeps=2000, n_equilibration=500, sample_interval=5
  - `Dict("E" => ..., "M" => ..., "M2" => ...)` を返す
  - **追加日**: 2026-03-11

#### 5.3.2 データ I/O (IO.jl)

- ✅ **TASK-5.3.6**: `save_result(filename, result::TemperatureSweepResult) -> Nothing`
  - JLD2 形式で温度掃引結果を保存（params, config, temperatures, results を含む）

- ✅ **TASK-5.3.7**: `save_result(filename, result::SingleTemperatureResult) -> Nothing`
  - JLD2 形式で単一温度結果を保存（多重ディスパッチ）

- ✅ **TASK-5.3.8**: `load_result(filename::String) -> Dict{String, Any}`
  - JLD2 ファイルを Dict として読み込む

#### 5.3.3 可視化 (Visualization.jl)

- ✅ **TASK-5.3.9**: `plot_spin_lattice(lattice::Matrix{Int}) -> Figure`
  - スピン配置のヒートマップ（青: -1, 赤: +1）

- ✅ **TASK-5.3.10**: `plot_thermodynamics(result::TemperatureSweepResult) -> Figure`
  - ⟨e⟩, ⟨|m|⟩, C, χ の 4 パネルプロット (800×600)

- ✅ **TASK-5.3.11**: `plot_binder_cumulant(results::Vector{TemperatureSweepResult}) -> Figure`
  - 複数格子サイズ L のビンダーキュムラント U_L を重ね描き
  - 有限サイズスケーリング解析に使用

- ✅ **TASK-5.3.12**: `plot_timeseries(energies, magnetizations) -> Figure`
  - エネルギー・磁化の MC 時系列プロット

### 5.4 スクリプト群

- ✅ **TASK-5.4.1**: `scripts/quick_start.jl` - 単一温度 MC の入門スクリプト
- ✅ **TASK-5.4.2**: `scripts/temperature_sweep.jl` - 固定パラメータ温度掃引
- ✅ **TASK-5.4.3**: `scripts/run_simulation.jl` - CLI 引数対応温度掃引
- ✅ **TASK-5.4.4**: `scripts/plot_results.jl` - 保存済み JLD2 の可視化
- ✅ **TASK-5.4.5**: `scripts/benchmark.jl` - ホットパスの BenchmarkTools 計測
- ✅ **TASK-5.4.6**: `scripts/plot_m2_vs_T.jl` - quick_mc を用いた ⟨M²⟩ vs T プロット
  - **追加日**: 2026-03-11
  - Tc ≈ 2.269 の破線を自動描画
  - `data/output/figs/m2_vs_T.png` に保存

---

## FDD Phase 6: テスト

### 6.1 ユニットテスト (test/Unit/ - 全体の 80%)

- ✅ **TASK-6.1.1**: `test/Unit/test_types.jl` - 型バリデーション
  - IsingParams: L ≥ 2, J ≠ 0 の検証
  - SimulationConfig: n_sweeps > 0, n_equilibration ≥ 0, sample_interval > 0 の検証
  - 不正値で ArgumentError が発生することを確認

- ✅ **TASK-6.1.2**: `test/Unit/test_observables.jl` - 観測量型の構築テスト
  - Observables の構築と値アクセス
  - ThermalAverages のゼロコンストラクタと mutable 蓄積
  - ThermodynamicQuantities の構築

- ✅ **TASK-6.1.3**: `test/Unit/test_results.jl` - 結果型の構築テスト
  - BoltzmannTable の構築
  - SingleTemperatureResult の構築
  - TemperatureSweepResult の構築

- ✅ **TASK-6.1.4**: `test/Unit/test_lattice.jl` - 格子生成テスト
  - random_lattice: サイズ (L, L)、値 ∈ {+1, -1}
  - uniform_lattice: 全スピン一様、spin ∉ {+1,-1} で AssertionError

- ✅ **TASK-6.1.5**: `test/Unit/test_physics.jl` - 物理量計算テスト
  - magnetization: 全上向き → 1.0、全下向き → -1.0、チェッカーボード → 0.0
  - energy: 全上向き → -2JN、チェッカーボード → +2JN
  - delta_energy: energy との整合性確認、周期境界条件の動作確認

- ✅ **TASK-6.1.6**: `test/Unit/test_montecarlo.jl` - MC 更新テスト
  - metropolis_step!: β=0 (T=∞) で常に受容、スピン値 ±1 の保持
  - sweep!: 受容数 ∈ [0, N]、T=∞ で全受容
  - quick_mc: 返り値の型と必須キー、物理的な値の範囲、低温/高温の ⟨M²⟩ 特性
  - quick_mc: 不正引数で ArgumentError
  - **quick_mc テスト追加日**: 2026-03-11

- ✅ **TASK-6.1.7**: `test/Unit/test_statistics.jl` - 熱力学量導出テスト
  - 既知の ThermalAverages から C, χ, U_L を手計算と比較
  - n_samples=0, N=0, T<0 で ArgumentError

### 6.2 統合テスト (test/Ising/)

- ✅ **TASK-6.2.1**: `test/Ising/test_integration.jl` - MC ワークフロー統合テスト
  - 平衡化 → 測定 → 熱力学量計算の一連のフロー
  - T < Tc での秩序相（|m| > 0.3, e < 0）の確認

### 6.3 E2E テスト (test/E2E/)

- ✅ **TASK-6.3.1**: `test/E2E/test_simulation.jl` - 高温/低温の物理特性 E2E 検証
  - T = 1.0 (T << Tc): ⟨|m|⟩ > 0.8 (秩序相)
  - T = 4.0 (T >> Tc): ⟨|m|⟩ < 0.3 (無秩序相)

---

## 機能要件別タスク

### F-LAT: 格子生成

| タスクID | 関数 | Pure/Impure |
|---------|------|------------|
| TASK-5.2.1 | `random_lattice` | Pure |
| TASK-5.2.2 | `uniform_lattice` | Pure |

### F-PHY: 物理量計算

| タスクID | 関数 | Pure/Impure |
|---------|------|------------|
| TASK-5.2.3 | `magnetization` | Pure |
| TASK-5.2.4 | `energy` | Pure |
| TASK-5.2.5 | `delta_energy` | Pure |

### F-MC: モンテカルロ法

| タスクID | 関数 | Pure/Impure |
|---------|------|------------|
| TASK-5.3.1 | `metropolis_step!` | Impure |
| TASK-5.3.2 | `sweep!` | Impure |
| TASK-5.3.3 | `run_single_temperature` | Impure |
| TASK-5.3.4 | `run_temperature_sweep` | Impure |
| TASK-5.3.5 | `quick_mc` | Impure |

### F-STAT: 熱力学量導出

| タスクID | 関数 | Pure/Impure |
|---------|------|------------|
| TASK-5.2.6 | `compute_thermodynamics` | Pure |

### F-IO: データ保存・読み込み

| タスクID | 関数 | Pure/Impure |
|---------|------|------------|
| TASK-5.3.6 | `save_result` (sweep) | Impure |
| TASK-5.3.7 | `save_result` (single) | Impure |
| TASK-5.3.8 | `load_result` | Impure |

### F-VIS: 可視化

| タスクID | 関数 | Pure/Impure |
|---------|------|------------|
| TASK-5.3.9 | `plot_spin_lattice` | Impure |
| TASK-5.3.10 | `plot_thermodynamics` | Impure |
| TASK-5.3.11 | `plot_binder_cumulant` | Impure |
| TASK-5.3.12 | `plot_timeseries` | Impure |

---

## 非機能要件タスク

### NFR-PERF: 性能要件

| タスクID | 内容 | 状態 |
|---------|------|------|
| TASK-5.3.5 | quick_mc: デフォルト 2000 スイープで数秒以内 | 完了 |
| TASK-5.4.5 | benchmark.jl で delta_energy / sweep! の計測 | 完了 |

### NFR-REPR: 再現性要件

| タスクID | 内容 | 状態 |
|---------|------|------|
| TASK-5.3.6 | JLD2 に params, config を含めて保存 | 完了 |
| TASK-3.1.1 | SimulationConfig に seed フィールドを定義 | 完了（実際の seed 適用は未実装） |

### NFR-MAINT: 保守性要件

| タスクID | 内容 | 状態 |
|---------|------|------|
| TASK-3.1.* | 全型にバリデーション付きコンストラクタ | 完了 |
| TASK-6.1.* | ユニットテストが全関数をカバー | 完了 |

---

## 今後の実装課題

以下の項目は現在未実装であり、今後の開発対象となる。

### 優先度: 高

- [ ] **FUTURE-1**: `SimulationConfig.seed` を `Random.seed!()` に実際に適用する
  - 現在はフィールドとして保存されているが、シミュレーション実行時に使われていない
  - 完全な再現性を保証するために必要

- [ ] **FUTURE-2**: ビンダーキュムラントを用いた有限サイズスケーリング解析スクリプト
  - 複数 L (例: 16, 32, 64) の U_L 交点から Tc を精密に推定
  - `plot_binder_cumulant` は実装済み、スクリプトが未作成

### 優先度: 中

- [ ] **FUTURE-3**: `BoltzmannTable` を用いた metropolis_step! の高速化
  - ΔE は {-8J, -4J, 0, 4J, 8J} の 5 値のみ（2次元イジングモデル）
  - exp(-βΔE) をルックアップテーブルで事前計算すると速い

- [ ] **FUTURE-4**: Wolff クラスターアルゴリズムの実装
  - 臨界温度付近でのスローダウン（臨界緩和）を緩和する
  - 現在のメトロポリス法を補完する高速アルゴリズム

- [ ] **FUTURE-5**: マルチスレッド対応
  - 独立した温度点を並列実行する (`Threads.@threads` または `Distributed`)

### 優先度: 低

- [ ] **FUTURE-6**: 自己相関時間の計算と統計誤差の推定
  - ブロッキング法または Autocorrelation Time による誤差評価

- [ ] **FUTURE-7**: 設定ファイル (TOML) 対応
  - シミュレーションパラメータを設定ファイルで管理

---

## 優先度と依存関係

### 実装済みコンポーネントの依存関係図

```
Types.jl (Pure)
    |
    +-- Ising/Lattice.jl (Pure)
    |       |
    +-- Ising/Physics.jl (Pure)
    |       |
    +-- Ising/Statistics.jl (Pure)
    |       |
    +-- Runtime/MonteCarlo.jl (Impure) <-- Lattice + Physics + Statistics
    |       |
    +-- Runtime/IO.jl (Impure) <-- Types
    |
    +-- Runtime/Visualization.jl (Impure) <-- Types
```

### 現在の完成状態サマリ

```
コア実装      ████████████████████ 100%
テスト        ████████████████████ 100%
スクリプト    ████████████████████ 100%
ベンチマーク  ████████████████████ 100%
高度な最適化  ░░░░░░░░░░░░░░░░░░░░  0%  (FUTURE-3,4,5)
誤差解析      ░░░░░░░░░░░░░░░░░░░░  0%  (FUTURE-6)
```

---

## 付録

### A. Pure / Impure 分類基準

| 分類 | 基準 | 例 |
|------|------|-----|
| **Pure** | 副作用なし、同じ入力 → 同じ出力 | `magnetization`, `compute_thermodynamics` |
| **Impure** | 副作用あり（乱数, IO, 状態変更） | `metropolis_step!`, `save_result` |

### B. 変更履歴

| 日付 | 変更内容 |
|------|---------|
| 2026-03-11 | `quick_mc` を `src/Runtime/MonteCarlo.jl` に追加 |
| 2026-03-11 | `quick_mc` を `src/Ising2D.jl` に export 追加 |
| 2026-03-11 | `scripts/plot_m2_vs_T.jl` を新規作成 |
| 2026-03-11 | `test/Unit/test_montecarlo.jl` に `quick_mc` テストを追加 |
| 2026-03-11 | `data/output/figs/m2_vs_T.png` を生成・追跡 |

### C. 関連ドキュメント

- [要件定義書](./spec/requirement.md)
- [システムアーキテクチャ設計書](./spec/system_architecture.md)
- [技術スタック仕様書](./spec/tech_stack.md)
- [データ構造定義書](./spec/data_structure.md)
- [ディレクトリ構造](./directory_structure.md)
- [FDD 実装ガイド](./src/fdd_implementation_guide.md)

---

*本文書は Ising2D プロジェクトの実際の実装状態を反映しています。*
*最終更新: 2026-03-11*
