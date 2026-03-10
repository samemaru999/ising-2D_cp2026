> **このドキュメントの役目**: 関数（集合間の写像）の定義書。型間の関係を関数シグネチャとして記述し、各関数がどのレイヤーに属しPure/Impureのどちらかを明示する。

# システムアーキテクチャ設計書

**プロジェクト名**: Ising2D
**最終更新日**: 2026-03-10

---

## 1. レイヤー構造

```
Application Layer    [Impure]  CLI/Notebook, 設定, 可視化
Logic Layer          [Pure]    温度掃引制御, 統計導出
Domain Layer         [Pure]    格子生成, 物理量計算, ΔE計算
Interpreters Layer   [Impure]  MC実行エンジン (乱数)
Persistence Layer    [Impure]  JLD2ファイル入出力
```

---

## 2. 関数定義（型間の写像）

### 2.1 格子生成（Domain, Pure）

| 関数 | 型シグネチャ | 説明 |
|------|-------------|------|
| `random_lattice` | `Int → Matrix{Int}` | ランダムスピン配置 |
| `uniform_lattice` | `Int × Int → Matrix{Int}` | 一様スピン配置 |

### 2.2 物理量（Domain, Pure）

| 関数 | 型シグネチャ | 説明 |
|------|-------------|------|
| `magnetization` | `Matrix{Int} → Float64` | m = (1/N) Σ sᵢ |
| `energy` | `Matrix{Int} × Float64 → Float64` | H = -J Σ sᵢsⱼ |
| `delta_energy` | `Matrix{Int} × Int × Int × Float64 → Float64` | ΔE = 2J s(i,j) Σnn |

### 2.3 MC更新（Interpreters, Impure）

| 関数 | 型シグネチャ | 説明 |
|------|-------------|------|
| `metropolis_step!` | `Matrix{Int} × Float64 × Float64 → Bool` | 1サイト反転試行 |
| `sweep!` | `Matrix{Int} × Float64 × Float64 → Int` | L²回の試行（1MCS） |

### 2.4 蓄積・導出（Logic, Pure）— 未実装

| 関数 | 型シグネチャ | 説明 |
|------|-------------|------|
| `accumulate` | `ThermalAverages × Observables → ThermalAverages` | 測定値の蓄積 |
| `compute_thermodynamics` | `ThermalAverages × Float64 × Int → ThermodynamicQuantities` | C, χ, U_L の導出 |

### 2.5 シミュレーション実行（Interpreters, Impure）— 未実装

| 関数 | 型シグネチャ | 説明 |
|------|-------------|------|
| `interpret` | `RunSimulation → ComputationResult{SingleTemperatureResult}` | 単一温度MC実行 |
| `interpret` | `RunTemperatureSweep → ComputationResult{TemperatureSweepResult}` | 温度掃引実行 |

### 2.6 性能最適化（Domain, Pure）— 未実装

| 関数 | 型シグネチャ | 説明 |
|------|-------------|------|
| `build_boltzmann_table` | `Float64 × Float64 → BoltzmannTable` | exp(-βΔE) の事前計算 |
| `metropolis_step!` (高速版) | `Matrix{Int} × BoltzmannTable → Bool` | テーブル参照による判定 |

### 2.7 IO（Persistence, Impure）— 未実装

| 関数 | 型シグネチャ | 説明 |
|------|-------------|------|
| `interpret` | `SaveResult → ComputationResult{Nothing}` | JLD2保存 |
| `interpret` | `LoadResult → ComputationResult{Any}` | JLD2読み込み |

---

## 3. Pure/Impure境界

```
Pure（入力のみに依存、副作用なし）
  格子生成, 物理量計算, ΔE計算, 統計導出, バリデーション
      │
      │ 境界: 乱数生成、ファイルIO
      ▼
Impure（副作用あり）
  metropolis_step!, sweep!, interpret, JLD2保存/読み込み
```

---

## 改訂履歴

| バージョン | 日付 | 変更内容 |
|-----------|------|---------|
| 1.0 | 2026-03-10 | 初版作成 |
| 2.0 | 2026-03-10 | 関数（型間の写像）を中心に改訂、実装状況を反映 |
