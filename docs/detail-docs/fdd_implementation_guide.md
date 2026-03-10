# FDD (Functional Declarative Design) と本プロジェクト実装の関係

このドキュメントでは、FDD（Functional Declarative Design）の思想と、Ising2Dプロジェクトの実際の実装がどのようにつながっているかを解説します。

---

## 目次

1. [FDDとは何か](#1-fddとは何か)
2. [eDSL（埋め込みドメイン特化言語）の実装](#2-edsl埋め込みドメイン特化言語の実装)
3. [interpret関数の仕組み](#3-interpret関数の仕組み)
4. [Result型とunwrap](#4-result型とunwrap)
5. [Pure/Impure分離の実際](#5-pureimpure分離の実際)
6. [全体の処理フロー図](#6-全体の処理フロー図)
7. [実践的なコード例](#7-実践的なコード例)

---

## 1. FDDとは何か

### 1.1 基本理念

FDD (Functional Declarative Design) は、関数型プログラミングの原則に基づいた設計手法です。主な目標は：

1. **本質的複雑さと偶有的複雑さの分離**
   - 本質的複雑さ：ドメイン固有の避けられない複雑さ（例：位相縮約理論の数学）
   - 偶有的複雑さ：実装選択による不必要な複雑さ（例：過度な抽象化）

2. **Pure（純粋）とImpure（非純粋）の分離**
   - 副作用のないコードと、副作用を持つコードを明確に分離

3. **宣言的記述**
   - 「何を計算するか」と「どのように計算するか」を分離

### 1.2 本プロジェクトにおけるFDDアーキテクチャ

```
┌─────────────────────────────────────────────────────────────┐
│                     利用者コード                              │
│    run_optimization, compare_strategies など                 │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│   ┌─────────────────────────────────────────────────────┐   │
│   │           Interpreters Layer [Impure]               │   │
│   │  SimulationInterpreter.jl, PhaseInterpreter.jl      │   │
│   │  - eDSLコマンドを実際の計算に変換                     │   │
│   │  - ODEソルバー呼び出し                               │   │
│   │  - Result型の生成                                    │   │
│   └─────────────────────────────────────────────────────┘   │
│                           ↑                                 │
│                       interpret()                           │
│                           ↓                                 │
│   ┌─────────────────────────────────────────────────────┐   │
│   │              Domain Layer [Pure]                    │   │
│   │  Types.jl, DSL.jl, Validation.jl, Constants.jl      │   │
│   │  - 型定義（NetworkParams, SPS, PSF等）               │   │
│   │  - eDSLコマンド定義（ComputeSPS, OptimizeLinear...） │   │
│   │  - 検証ロジック                                      │   │
│   └─────────────────────────────────────────────────────┘   │
│                           ↑                                 │
│                    Pure関数呼び出し                          │
│                           ↓                                 │
│   ┌─────────────────────────────────────────────────────┐   │
│   │              Logic Layer [Pure]                     │   │
│   │  FHN.jl, Coupling.jl, PhaseReduction.jl,           │   │
│   │  Optimization.jl                                    │   │
│   │  - 数学的アルゴリズム                                │   │
│   │  - 純粋な計算関数                                    │   │
│   └─────────────────────────────────────────────────────┘   │
│                                                             │
├─────────────────────────────────────────────────────────────┤
│   ┌─────────────────────────────────────────────────────┐   │
│   │           Runtime/IO Layer [Impure]                 │   │
│   │  IOHandler.jl                                       │   │
│   │  - ファイル保存・読み込み                             │   │
│   └─────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

---

## 2. eDSL（埋め込みドメイン特化言語）の実装

### 2.1 eDSLとは

eDSL（embedded Domain Specific Language）は、ホスト言語（Julia）内で実装される特定ドメイン向けの小さな言語です。

**核心的なアイデア**：計算処理を「実行するもの」ではなく「データとして表現するもの」として扱う。

```julia
# 従来のアプローチ：直接実行
sps = compute_stable_periodic_solution(params, 101, 1000.0)  # 即座に計算実行

# eDSLアプローチ：まずコマンドを「データ」として作成
cmd = ComputeSPS(params, 101, 1000.0)  # 計算の「宣言」のみ、実行はまだ
result = interpret(cmd)                  # ここで初めて実行
```

### 2.2 実装箇所：`src/Domain/DSL.jl`

このファイルにすべてのeDSLコマンドが定義されています。

```julia
# ============================
# コマンドの基底型
# ============================
abstract type ComputationCommand end

# ============================
# 具体的なコマンド定義
# ============================

# 安定周期解計算コマンド
struct ComputeSPS <: ComputationCommand
    params::NetworkParams
    Ntheta::Int
    t_relax::Float64
end

# 位相感受関数計算コマンド
struct ComputePSF <: ComputationCommand
    sps::StablePeriodicSolution
    n_iterations::Int
    coupling_type::CouplingType
end

# 最適化コマンド
struct OptimizeLinearStability <: OptimizationCommand
    individual_pcfs::Vector{IndividualPCF}
    target_q::Float64
    N::Int
    Ntheta::Int
end

# ... その他のコマンド
```

### 2.3 コマンドが「データ」である利点

| 従来のアプローチ | eDSLアプローチ |
|------------------|----------------|
| 関数を呼ぶと即座に計算 | コマンドは単なるデータ構造 |
| 実行前の検証が困難 | 実行前にコマンドを検査可能 |
| テストにはモックが必要 | コマンド生成ロジックを純粋にテスト可能 |
| 計算の記録が困難 | コマンドをログや再現に使用可能 |

```julia
# 例：コマンドの検査（実行前）
cmd = ComputeSPS(params, 101, 1000.0)
println(cmd.Ntheta)     # 101 - 設定を確認可能
println(cmd.t_relax)    # 1000.0

# コマンドのリストとして計算パイプラインを表現
script = [
    ComputeSPS(params, 101, 1000.0),
    # PSF, IndividualPCF等は動的に生成
]
```

### 2.4 Deep Embedding vs Shallow Embedding

本プロジェクトは**Deep Embedding**を採用しています：

```julia
# Deep Embedding（本プロジェクト）
# コマンドはデータ構造として中間表現を形成
cmd = ComputeSPS(params)  # ASTノードを作成
result = interpret(cmd)    # 後でインタープリターが解釈

# Shallow Embedding（参考）
# 操作が即座に実行される
result = compute_sps(params)  # 直接実行
```

---

## 3. interpret関数の仕組み

### 3.1 インタープリターパターンの概念

インタープリターパターンは、**コマンド（データ）を受け取り、実際の処理を実行する**デザインパターンです。

```
┌─────────────┐       interpret()        ┌─────────────┐
│   コマンド    │  ─────────────────────>  │   結果       │
│   (データ)    │                         │ (Result型)   │
└─────────────┘                          └─────────────┘
     Pure                                      ↓
                                         Success/Failure
```

### 3.2 実装箇所：`src/Interpreters/`

```julia
# SimulationInterpreter.jl
function interpret(cmd::ComputeSPS)
    try
        # Impure: ODEソルバー呼び出し
        sps = compute_stable_periodic_solution(cmd.params, cmd.Ntheta, cmd.t_relax)
        return Success(sps)
    catch e
        return Failure{StablePeriodicSolution}("SPS computation failed: $e")
    end
end

# PhaseInterpreter.jl
function interpret(cmd::ComputeIndividualPCF)
    try
        # Pure関数を呼び出し
        pcf = compute_individual_pcf(cmd.sps, cmd.psf, cmd.i, cmd.j, cmd.k)
        return Success(pcf)
    catch e
        return Failure{IndividualPCF}("IndividualPCF computation failed: $e")
    end
end
```

### 3.3 多重ディスパッチによるインタープリターの実装

Juliaの多重ディスパッチを活用し、コマンドの型に応じて適切な`interpret`関数が自動選択されます：

```julia
# 同じ関数名 interpret だが、引数の型で処理が分かれる
interpret(cmd::ComputeSPS)              # → SPS計算
interpret(cmd::ComputePSF)              # → PSF計算
interpret(cmd::ComputeIndividualPCF)    # → IndividualPCF計算
interpret(cmd::OptimizeLinearStability) # → 線形安定性最適化
interpret(cmd::SaveData)                # → データ保存
```

これにより、新しいコマンドを追加する際は：
1. `DSL.jl`に新しいコマンド型を定義
2. 対応する`interpret`関数を追加

既存コードを変更せずに拡張できます（Open-Closed Principle）。

### 3.4 なぜインタープリターを使うのか

| 利点 | 説明 |
|------|------|
| **副作用の隔離** | Impureな処理（ODE計算、IO）をインタープリター層に集約 |
| **テスタビリティ** | ビジネスロジックをPure関数として独立テスト可能 |
| **拡張性** | 新コマンドの追加が容易 |
| **エラーハンドリング** | 統一されたResult型で一貫したエラー処理 |
| **ログ・監視** | インタープリターレベルで実行を監視可能 |

---

## 4. Result型とunwrap

### 4.1 Result型（Either型）の概念

Result型は、関数型言語で一般的なエラーハンドリングパターンです。計算結果が「成功」か「失敗」かを型で表現します。

```julia
# 従来の例外ベースのエラーハンドリング
try
    result = risky_computation()
    # 成功時の処理
catch e
    # エラー時の処理
end

# Result型ベースのエラーハンドリング
result = risky_computation()  # Result{T}を返す
if is_success(result)
    value = unwrap(result)    # 成功時：値を取り出す
else
    error_msg = result.error  # 失敗時：エラー情報を取得
end
```

### 4.2 実装箇所：`src/Domain/DSL.jl`

```julia
# 計算結果を表す型（Either風）
abstract type ComputationResult{T} end

# 成功結果
struct Success{T} <: ComputationResult{T}
    value::T
end

# 失敗結果
struct Failure{T} <: ComputationResult{T}
    error::String
    context::Dict{String, Any}
end

# ユーティリティ関数
is_success(::Success) = true
is_success(::Failure) = false

function unwrap(r::Success{T}) where T
    return r.value  # 値を取り出す
end

function unwrap(r::Failure{T}) where T
    throw(ErrorException("Unwrap failed: $(r.error)"))  # エラーをスロー
end

function unwrap_or(r::Success{T}, default::T) where T
    return r.value
end

function unwrap_or(r::Failure{T}, default::T) where T
    return default  # 失敗時はデフォルト値
end
```

### 4.3 unwrapの使い方

```julia
# 基本的な使用パターン
result = interpret(ComputeSPS(params))

# パターン1: 明示的なチェック
if is_success(result)
    sps = unwrap(result)
    println("周期: $(sps.T)")
else
    println("エラー: $(result.error)")
end

# パターン2: デフォルト値を使用
sps = unwrap_or(result, default_sps)

# パターン3: 連鎖的な処理（Railway-Oriented Programming）
result = interpret(ComputeSPS(params)) >>
         (sps -> interpret(ComputePSF(sps))) >>
         (psf -> interpret(ComputeIndividualPCF(sps, psf, 1, 1, 1)))
```

### 4.4 モナド的操作

Result型は**モナド**として振る舞い、エラー処理を合成可能にします：

```julia
# map_result: 成功値に関数を適用（失敗はそのまま伝播）
result2 = map_result(sps -> sps.T, result)
# Success(42.5) or Failure{Float64}("...")

# bind_result (>>演算子): Result返す関数を連鎖
result3 = result >> (sps -> interpret(ComputePSF(sps)))
# 最初のエラーで処理が停止し、そのエラーが伝播
```

### 4.5 なぜ例外ではなくResult型を使うのか

| 例外ベース | Result型ベース |
|------------|----------------|
| エラーが暗黙的 | エラーが型で明示 |
| 呼び出し側でtry-catchが必要 | 戻り値の型でエラーの可能性が分かる |
| 制御フローが予測困難 | 線形な制御フロー |
| テストが複雑 | テストが容易 |

---

## 5. Pure/Impure分離の実際

### 5.1 分離の原則

- **Pure（純粋）関数**：同じ入力に対して常に同じ出力、副作用なし
- **Impure（非純粋）関数**：外部状態に依存、または副作用を持つ

### 5.2 本プロジェクトでの分離

```
Pure Layer（副作用なし）
├── src/Domain/Types.jl      # データ型定義
├── src/Domain/DSL.jl        # コマンド定義
├── src/Domain/Validation.jl # 検証ロジック
├── src/Logic/FHN.jl         # FHNモデルの方程式
├── src/Logic/Coupling.jl    # 結合関数
├── src/Logic/PhaseReduction.jl # 位相縮約理論
└── src/Logic/Optimization.jl   # 最適化アルゴリズム

Impure Layer（副作用あり）
├── src/Interpreters/SimulationInterpreter.jl  # ODEソルバー呼び出し
├── src/Interpreters/PhaseInterpreter.jl       # インタープリター
└── src/Runtime/IO/IOHandler.jl                # ファイルIO
```

### 5.3 具体例：FHN.jl（Pure）

```julia
# src/Logic/FHN.jl - 完全にPure

# 純粋関数：同じ入力 → 同じ出力
function fhn_single_dynamics(u::Float64, v::Float64, I::Float64, fhn::FHNParams)
    du = fhn.delta * (fhn.a + v - fhn.b * u)
    dv = v - v^3 / 3 - u + I
    return (du, dv)
end

# ヤコビアン計算も純粋関数
function fhn_single_jacobian(u::Float64, v::Float64, fhn::FHNParams)
    J = zeros(2, 2)
    J[1, 1] = -fhn.delta * fhn.b
    J[1, 2] = fhn.delta
    J[2, 1] = -1.0
    J[2, 2] = 1.0 - v^2
    return J
end
```

### 5.4 具体例：SimulationInterpreter.jl（Impure）

```julia
# src/Interpreters/SimulationInterpreter.jl - Impure

function interpret(cmd::ComputeSPS)
    try
        # Impure: ODE問題の構築と解法
        prob = ODEProblem(dynamics!, X0, (0.0, t_relax))
        sol = solve(prob, Tsit5())  # ← 副作用（数値計算）

        # Pure関数を呼び出してデータ変換
        sps = create_sps(sol, params)

        return Success(sps)
    catch e
        return Failure{StablePeriodicSolution}("...")
    end
end
```

### 5.5 分離の利点

```julia
# Logic層の関数は単独でテスト可能
@testset "Pure Functions" begin
    fhn = FHNParams(0.08, 0.7, 0.8)
    du, dv = fhn_single_dynamics(0.0, 0.0, 0.5, fhn)
    @test du ≈ 0.056  # 0.08 * (0.7 + 0.0 - 0.8*0.0)
    @test dv ≈ 0.5    # 0.0 - 0.0 - 0.0 + 0.5
end

# モック不要、セットアップ不要、高速
```

---

## 6. 全体の処理フロー図

### 6.1 典型的な計算フロー

```
利用者コード
     │
     │  run_optimization(4; coupling_type=:linear, target_value=-0.1)
     ▼
┌────────────────────────────────────────────────────────────────────┐
│  High-Level API (Ising2D.jl)                                     │
│                                                                    │
│  1. params = NetworkParams(4)                                      │
│  2. results = run_full_pipeline(params; coupling_type=:linear, ..) │
└────────────────────────────────────────────────────────────────────┘
                     │
                     ▼
┌────────────────────────────────────────────────────────────────────┐
│  run_full_pipeline (PhaseInterpreter.jl)                          │
│                                                                    │
│  Step 1: sps_cmd = ComputeSPS(params)         ← コマンド作成      │
│          sps_result = interpret(sps_cmd)       ← 実行             │
│          sps = unwrap(sps_result)              ← 値取り出し       │
│                                                                    │
│  Step 2: psf_cmd = ComputePSF(sps)            ← コマンド作成      │
│          psf_result = interpret(psf_cmd)       ← 実行             │
│          psf = unwrap(psf_result)              ← 値取り出し       │
│                                                                    │
│  Step 3: individual_pcfs = compute_all_...     ← 64個のPCF計算    │
│                                                                    │
│  Step 4: opt_cmd = OptimizeLinearStability(...) ← コマンド作成    │
│          opt_result = interpret(opt_cmd)        ← 実行            │
│          opt = unwrap(opt_result)               ← 値取り出し      │
│                                                                    │
│  Step 5: pcf_cmd = ComputeTotalPCF(...)        ← コマンド作成     │
│          total_pcf = unwrap(interpret(pcf_cmd)) ← 実行 & 取り出し │
│                                                                    │
│  return (sps=sps, psf=psf, total_pcf=total_pcf, metrics=..., ...)  │
└────────────────────────────────────────────────────────────────────┘
                     │
                     ▼
┌────────────────────────────────────────────────────────────────────┐
│  interpret(cmd::ComputeSPS) [SimulationInterpreter.jl]            │
│                                                                    │
│  try                                                               │
│      # Impure: ODEソルバー呼び出し                                  │
│      prob = ODEProblem(dynamics!, X0, tspan)                       │
│      sol = solve(prob, Tsit5())                                    │
│                                                                    │
│      # Pure: データ構造作成                                         │
│      sps = StablePeriodicSolution(N, Ntheta, T, omega, Xs, params) │
│                                                                    │
│      return Success(sps)                                           │
│  catch e                                                           │
│      return Failure{StablePeriodicSolution}("...")                 │
│  end                                                               │
└────────────────────────────────────────────────────────────────────┘
```

### 6.2 データの流れ

```
NetworkParams (入力)
      │
      ▼ interpret(ComputeSPS)
StablePeriodicSolution (SPS)
      │
      ▼ interpret(ComputePSF)
PhaseSensitivityFunction (PSF)
      │
      ▼ interpret(ComputeIndividualPCF) × N³
Vector{IndividualPCF}
      │
      ▼ interpret(OptimizeLinearStability)
OptimizationResult (C_opt)
      │
      ▼ interpret(ComputeTotalPCF)
TotalPCF (Γ関数)
      │
      ▼ compute_stability_metrics
StabilityMetrics (Λ, R, is_stable)
```

---

## 7. 実践的なコード例

### 7.1 最も簡単な使用例

```julia
using Ising2D

# 高レベルAPIを使用（内部でeDSL + interpretが動作）
results = run_optimization(4; coupling_type=:linear, target_value=-0.1)

println("周期: $(results.sps.T)")
println("安定性指標 Λ: $(results.metrics.Lambda)")
```

### 7.2 ステップバイステップの使用例

```julia
using Ising2D

# 1. パラメータ設定（Domain層）
params = NetworkParams(4)

# 2. コマンドを作成して実行（eDSL + Interpreter）
sps_cmd = ComputeSPS(params, 101, 500.0)
sps_result = interpret(sps_cmd)

# 3. 結果を確認してunwrap
if is_success(sps_result)
    sps = unwrap(sps_result)
    println("周期 T = $(sps.T)")
else
    println("エラー: $(sps_result.error)")
end

# 4. PSF計算
psf_cmd = ComputePSF(sps)
psf = unwrap(interpret(psf_cmd))

# 5. IndividualPCF計算
pcfs = IndividualPCF[]
for i in 1:4, j in 1:4, k in 1:4
    cmd = ComputeIndividualPCF(sps, psf, i, j, k)
    pcf = unwrap(interpret(cmd))
    push!(pcfs, pcf)
end

# 6. 最適化
opt_cmd = OptimizeLinearStability(pcfs, -0.1, 4, 101)
opt = unwrap(interpret(opt_cmd))

println("最適化された結合テンソルのFrobenius Norm: $(opt.frobenius_norm)")
```

### 7.3 エラーハンドリングの例

```julia
# パターン1: 明示的なチェック
result = interpret(ComputeSPS(params))
if !is_success(result)
    @error "SPS計算失敗" error=result.error
    return nothing
end
sps = unwrap(result)

# パターン2: unwrap_orでデフォルト値
sps = unwrap_or(interpret(ComputeSPS(params)), cached_sps)

# パターン3: モナド的チェイン
final_result = interpret(ComputeSPS(params)) >>
    sps -> interpret(ComputePSF(sps)) >>
    psf -> Success((sps, psf))  # 両方を保持

if is_success(final_result)
    sps, psf = unwrap(final_result)
end
```

### 7.4 カスタムコマンドの追加方法

新しい計算処理を追加する場合：

```julia
# Step 1: DSL.jlに新しいコマンドを定義
struct MyCustomComputation <: ComputationCommand
    input_data::SomeType
    parameter::Float64
end

# Step 2: 対応するinterpret関数を追加
function interpret(cmd::MyCustomComputation)
    try
        # Pure関数を呼び出し
        result = my_pure_computation(cmd.input_data, cmd.parameter)
        return Success(result)
    catch e
        return Failure{ResultType}("MyCustomComputation failed: $e")
    end
end

# Step 3: Logic層にPure関数を実装
function my_pure_computation(data::SomeType, param::Float64)::ResultType
    # 純粋な計算処理
end
```

---

## まとめ

| FDD概念 | 本プロジェクトでの実装 |
|---------|------------------------|
| **eDSL** | `src/Domain/DSL.jl` の `ComputationCommand` とその派生型 |
| **インタープリター** | `src/Interpreters/` の `interpret` 関数群 |
| **Result型** | `Success{T}`, `Failure{T}` と `unwrap`, `is_success` 等 |
| **Pure/Impure分離** | Domain/Logic層 = Pure, Interpreters/Runtime層 = Impure |

この設計により：
- **テスタビリティ**：Logic層の関数は単独でテスト可能
- **拡張性**：新コマンドの追加が容易
- **保守性**：責務の分離により変更の影響範囲が限定
- **再現性**：コマンドとしての計算記録が可能

---

## 参考資料

- [Domain-specific language - Wikipedia](https://en.wikipedia.org/wiki/Domain-specific_language)
- [Interpreter Design Pattern - GeeksforGeeks](https://www.geeksforgeeks.org/system-design/interpreter-design-pattern/)
- [Pure vs Impure Functions - freeCodeCamp](https://www.freecodecamp.org/news/pure-function-vs-impure-function/)
- [CLAUDE.md](../CLAUDE.md) - 本プロジェクトのFDD方法論ガイドライン
