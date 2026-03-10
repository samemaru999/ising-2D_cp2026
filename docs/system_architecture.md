# システムアーキテクチャ設計書

**プロジェクト名**: 2次元イジングモデル モンテカルロシミュレーション (Ising2D)
**バージョン**: 1.0
**作成日**: 2026-03-10
**設計手法**: Functional Declarative Design (FDD)

---

## 1. システム概要

### 1.1 目的

本システムは、2次元正方格子上のイジングモデルに対してメトロポリス法によるモンテカルロシミュレーションを実行し、相転移現象の解析・可視化を行うJuliaベースのシミュレーションツールである。

### 1.2 対象システムの数理モデル

```
ハミルトニアン:
H = -J Σ_{<i,j>} s_i s_j - h Σ_i s_i

スピン変数: s_i ∈ {+1, -1}
格子: L x L 正方格子、周期的境界条件
臨界温度: T_c = 2J / (k_B ln(1 + √2)) ≈ 2.269 J/k_B
```

### 1.3 主要な処理フロー

```
Phase 1: 格子構成
  格子生成 → 初期スピン配置 → パラメータ設定

Phase 2: モンテカルロシミュレーション
  熱平衡化 → サンプリング → 物理量計測

Phase 3: 解析
  統計処理 → 有限サイズスケーリング → 臨界指数推定

Phase 4: 可視化・出力
  スピン配置表示 → 物理量プロット → データ保存
```

---

## 2. FDD設計理念

### 2.1 本質的複雑さ vs 偶有的複雑さ

| 種類 | 定義 | 本プロジェクトでの例 |
|------|------|---------------------|
| **本質的複雑さ** | ドメイン固有の避けられない複雑さ | ハミルトニアン計算、メトロポリス法、統計力学的平均 |
| **偶有的複雑さ** | 実装選択による不必要な複雑さ | 過度な抽象化、不要な状態管理 |

### 2.2 FDD 5つの設計原則

| 原則 | 本プロジェクトでの適用 |
|------|----------------------|
| **Pure/Impure分離** | ハミルトニアン計算、物理量計算はPure。乱数生成、ファイルIOはImpure |
| **eDSLによる抽象化** | シミュレーションコマンドをデータとして表現 |
| **Service Handleパターン** | サービス関数の集合による抽象化 |
| **トップダウン反復** | 要件→設計→実装のフェーズ反復 |
| **型による設計** | スピン状態、格子構成、シミュレーションパラメータの型定義 |

---

## 3. FDDアーキテクチャレイヤー

### 3.1 レイヤー構造図

```
+-------------------------------------------------------------+
|                    Application Layer                         |
|   設定管理, ロギング, CLI/REPL/Notebook                       |
|   [src/App/] - Impure                                        |
+-------------------------------------------------------------+
|                    Business Logic Layer                      |
|   シミュレーションシナリオ, 温度掃引制御                       |
|   [src/Logic/] - Pure                                        |
+-------------------------------------------------------------+
|                    Service Layer                             |
|   機能インターフェース, eDSL定義                              |
|   [src/Services/] - Mixed                                    |
+-------------------------------------------------------------+
|                    Domain Model Layer                        |
|   ドメイン型定義, ハミルトニアン, 物理量計算                   |
|   [src/Domain/] - 完全にPure                                 |
+-------------------------------------------------------------+
|                    Interoperability Layer                    |
|   eDSLインタープリター, モンテカルロ実行エンジン               |
|   [src/Interpreters/] - Impure                               |
+-------------------------------------------------------------+
|                    Persistence Layer                         |
|   JLD2ファイル入出力, データシリアライゼーション                |
|   [src/Runtime/] - Impure                                    |
+-------------------------------------------------------------+
```

---

## 4. Pure/Impure分離

### 4.1 分離構造

```
Pure Layer (宣言的、テスト容易、参照透過)
+-------------------------------------------+
|  Domain/Types.jl       # スピン状態、格子型 |
|  Domain/Hamiltonian.jl # エネルギー計算     |
|  Domain/Observables.jl # 物理量計算         |
|  Domain/Validation.jl  # パラメータ検証     |
|  Logic/Scenarios.jl    # シナリオ定義       |
+-------------------------------------------+
               |
               v
Impure Layer (副作用あり、最小限に保つ)
+-------------------------------------------+
|  Interpreters/MCInterpreter.jl  # MC実行   |
|  Runtime/IOHandler.jl           # ファイルIO|
|  App/Main.jl                    # エントリ  |
+-------------------------------------------+
```

### 4.2 Pure関数の例

```julia
# ハミルトニアン計算 (Pure - 入力のみに依存)
function compute_energy(spins::Matrix{Int8}, J::Float64, h::Float64)::Float64
    L = size(spins, 1)
    E = 0.0
    for i in 1:L, j in 1:L
        s = spins[i, j]
        # 右と下の隣接スピンのみカウント（重複回避）
        s_right = spins[i, mod1(j+1, L)]
        s_down = spins[mod1(i+1, L), j]
        E -= J * s * (s_right + s_down)
        E -= h * s
    end
    return E
end

# エネルギー変化の計算 (Pure)
function compute_delta_E(spins::Matrix{Int8}, i::Int, j::Int, J::Float64, h::Float64)::Float64
    L = size(spins, 1)
    s = spins[i, j]
    nn_sum = spins[mod1(i-1,L), j] + spins[mod1(i+1,L), j] +
             spins[i, mod1(j-1,L)] + spins[i, mod1(j+1,L)]
    return 2.0 * s * (J * nn_sum + h)
end
```

---

## 5. サブシステム設計

### 5.1 サブシステム構成

| サブシステム | 責務 | Pure/Impure |
|-------------|------|------------|
| Lattice | 格子生成、近傍関係、スピン配置操作 | Pure |
| Hamiltonian | エネルギー計算、ΔE計算 | Pure |
| MonteCarlo | メトロポリス法の実行、熱平衡化 | Impure (乱数) |
| Observables | 磁化、比熱、磁化率、ビンダーキュムラント | Pure |
| Statistics | 統計平均、エラー推定、ブロック平均法 | Pure |
| FiniteSizeScaling | 臨界温度推定、臨界指数推定 | Pure |
| Visualization | スピン配置表示、物理量プロット | Impure (描画) |
| IO | JLD2保存・読み込み | Impure (ファイル) |

### 5.2 eDSL設計

```julia
# コマンドをデータとして表現 (Pure)
abstract type SimCommand end

struct InitializeLattice <: SimCommand
    L::Int
    init_config::Symbol  # :all_up, :all_down, :random
end

struct Thermalize <: SimCommand
    n_steps::Int
end

struct Sample <: SimCommand
    n_steps::Int
    sampling_interval::Int
end

struct ComputeObservables <: SimCommand end

struct TemperatureSweep <: SimCommand
    T_min::Float64
    T_max::Float64
    dT::Float64
end

const Scenario = Vector{SimCommand}
```

---

## 6. 外部パッケージ統合

| パッケージ | レイヤー | 用途 |
|-----------|---------|------|
| Random (stdlib) | Interpreters | 乱数生成（メトロポリス法） |
| Statistics (stdlib) | Logic | 統計量計算 |
| LinearAlgebra (stdlib) | Domain | 行列演算 |
| JLD2.jl | Runtime | データ永続化 |
| GLMakie.jl | App | 可視化 |
| ProgressMeter.jl | App | 進捗表示 |

---

## 7. 実行環境

| モード | 用途 | 起動方法 |
|--------|------|---------|
| REPL | インタラクティブ開発 | `julia --project=.` |
| Script | バッチ実行 | `julia --project=. scripts/run.jl` |
| Notebook | 探索的解析 | `jupyter notebook notebooks/` |
| Test | テスト実行 | `julia --project=. -e 'using Pkg; Pkg.test()'` |

---

## 付録

### A. 関連ドキュメント

- [要件定義書](./requirement.md)
- [技術スタック仕様書](./tech_stack.md)
- [データ構造定義書](./data_structure.md)
- [ディレクトリ構造](./directory_structure.md)
- [実装タスクチェックリスト](./implementation_task_checklist.md)

### B. 参考文献

1. L. Onsager, "Crystal Statistics. I.", Phys. Rev. 65, 117 (1944)
2. N. Metropolis et al., "Equation of State Calculations by Fast Computing Machines", J. Chem. Phys. 21, 1087 (1953)

### C. 用語集

| 用語 | 定義 |
|------|------|
| FDD | Functional Declarative Design - 関数型宣言的設計 |
| Pure関数 | 副作用なし、同じ入力に対して常に同じ出力を返す関数 |
| Impure関数 | 副作用を持つ関数（IO、乱数、状態変更など） |
| eDSL | embedded Domain Specific Language - 埋め込みドメイン固有言語 |
| メトロポリス法 | 詳細釣り合い条件を満たすMCMCアルゴリズム |
