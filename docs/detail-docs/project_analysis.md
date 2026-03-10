# Ising2D プロジェクト解析とAPIリファレンス

このドキュメントでは、`Ising2D` パッケージの詳細な解析結果と、APIの使用方法について解説します。

## 1. プロジェクト概要

**Ising2D** は、高次結合（Higher-Order Coupling）を持つ集団振動ネットワークの同期現象を解析・最適化するための Julia パッケージです。
**Functional Declarative Design (FDD)** アーキテクチャを採用しており、純粋な数学的定義（Domain/Logic）と、副作用を伴う計算処理（Interpreters/Runtime）が厳密に分離されています。

## 2. アーキテクチャ詳細

本プロジェクトは以下の4層構造で構成されています。

### Layer 1: Domain Layer [Pure]
- **役割**: ドメインモデル、定数、eDSL（コマンド定義）、バリデーションロジックの定義。
- **特徴**: 副作用なし、他レイヤーへの依存なし。
- **主要ファイル**:
    - `src/Domain/Types.jl`: `NetworkParams`, `StablePeriodicSolution` などの型定義。
    - `src/Domain/DSL.jl`: `ComputeSPS`, `OptimizeLinearStability` などのコマンド定義。

### Layer 2: Logic Layer [Pure]
- **役割**: コアとなる数学的アルゴリズム、ダイナミクス、物理法則の実装。
- **特徴**: 純粋関数のみで構成。
- **主要ファイル**:
    - `src/Logic/FHN.jl`: FitzHugh-Nagumoモデルの微分方程式。
    - `src/Logic/PhaseReduction.jl`: 位相縮約、安定性指標の計算式。

### Layer 3: Interpreters Layer [Impure]
- **役割**: eDSLコマンドの解釈と実行。ODEソルバーの呼び出しなど、計算負荷の高い処理を担当。
- **特徴**: DomainとLogicを使用し、計算結果を生成する。
- **主要ファイル**:
    - `src/Interpreters/SimulationInterpreter.jl`: 微分方程式の数値積分。
    - `src/Interpreters/PhaseInterpreter.jl`: 位相感受関数や結合関数の計算パイプライン。

### Layer 4: Runtime Layer [Impure]
- **役割**: ファイルI/Oなどの外部とのやり取り。
- **特徴**: システムの状態変更や永続化を担当。
- **主要ファイル**:
    - `src/Runtime/IO/IOHandler.jl`: JLD2形式でのデータ保存・読み込み。

---

## 3. 主要データ構造 (API Reference)

ユーザーが直接扱う主要なデータ型 (`src/Domain/Types.jl`) です。

### `NetworkParams`
ネットワーク設定を保持する構造体。
```julia
struct NetworkParams
    N::Int                    # 振動子数
    fhn::FHNParams           # FHNパラメータ (δ, a, b)
    I::Vector{Float64}       # 外部入力
    K::Matrix{Float64}       # ネットワーク内結合行列
    topology::TopologyType   # トポロジー
end
```
**作成方法**:
```julia
# デフォルト設定で作成（N=4）
params = NetworkParams(4)
```

### `StablePeriodicSolution` (SPS)
安定周期解の計算結果。
```julia
struct StablePeriodicSolution
    N::Int
    Ntheta::Int             # 位相の離散化数
    T::Float64              # 周期
    omega::Float64          # 角周波数
    Xs::Matrix{Float64}     # リミットサイクル軌道
    params::NetworkParams
end
```

### `TotalPCF`
全体位相結合関数 $\Gamma(\phi, \psi)$。
```julia
struct TotalPCF
    N::Int
    Gamma::Matrix{Float64}  # Γ関数
    C::Array{Float64, 3}    # 結合テンソル C_ijk
    dGamma_dphi::Float64    # ∂Γ/∂φ (0,0)
    dGamma_dpsi::Float64    # ∂Γ/∂ψ (0,0)
    # ...
end
```

---

## 4. eDSL コマンド (API Reference)

計算処理は「コマンド」として表現され、`interpret` 関数に渡されて実行されます。

### 基本的な計算フロー
1. **`ComputeSPS`**: 安定周期解を計算
2. **`ComputePSF`**: 位相感受関数 (Adjoint) を計算
3. **`ComputeIndividualPCF`**: 個別の結合関数 $\gamma_{ijk}$ を計算
4. **`OptimizeLinearStability` / `OptimizeRotation`**: 最適な結合テンソル $C$ を計算

### コマンド例

#### `ComputeSPS`
```julia
cmd = ComputeSPS(params, 101, 1000.0)
# params: NetworkParams
# 101: 位相分割数 (Ntheta)
# 1000.0: 緩和時間 (t_relax)
```

#### `ComputeIndividualPCF`
特定の $i, j, k$ トリプレットに対する結合関数を計算します。
```julia
cmd = ComputeIndividualPCF(sps, psf, i, j, k)
```

#### `ObserveLinearStability`
線形安定性指標 $\Lambda$ を最大化（負の方向に大きく）する結合テンソルを求めます。
```julia
cmd = OptimizeLinearStability(individual_pcfs, target_q, N, Ntheta)
# target_q: 目標とする Γ1 + Γ2 の値
```

---

## 5. 使用方法ガイド

### 高レベル API の使用
最も簡単な使用方法は、`Ising2D.jl` で定義されている高レベル関数を使用することです。

#### 最適化の実行 (`run_optimization`)
全パイプライン（SPS計算 $\to$ 最適化）を一括で実行します。

```julia
using Ising2D

# N=4 のネットワークで、線形安定性最適化を実行
# target_value は目標とする安定性指標
results = run_optimization(4; 
    coupling_type=:linear, 
    target_value=-0.1, 
    save_results=true
)
```

#### 戦略比較 (`compare_strategies`)
一様結合、線形最適化、回転最適化の性能を比較します。

```julia
using Ising2D

results = compare_strategies(4)
```

### 低レベル API の使用 (ステップバイステップ)
詳細な制御が必要な場合は、コマンドを個別に作成して `interpret` します。

```julia
using Ising2D

# 1. パラメータ設定
params = NetworkParams(4)

# 2. 安定周期解の計算
result_sps = interpret(ComputeSPS(params))
if !is_success(result_sps)
    error("SPS calculation failed")
end
sps = unwrap(result_sps)

# 3. 位相感受関数の計算
result_psf = interpret(ComputePSF(sps))
psf = unwrap(result_psf)

# 4. 個別PCFの計算 (全組み合わせ)
pcfs = IndividualPCF[]
for i in 1:4, j in 1:4, k in 1:4
    cmd = ComputeIndividualPCF(sps, psf, i, j, k)
    push!(pcfs, unwrap(interpret(cmd)))
end

# 5. 線形安定性の最適化
cmd_opt = OptimizeLinearStability(pcfs, -0.2, 4, 101)
result_opt = interpret(cmd_opt)
opt_result = unwrap(result_opt)

println("Optimized Frobenius Norm: $(opt_result.frobenius_norm)")
```

## 6. IO とデータ保存

計算結果は `Runtime/IO/IOHandler.jl` を通じて JLD2 ファイルとして保存可能です。

- `Generic save`: `interpret(SaveData(data, "filename.jld2"))`
- `Specific helpers`: `save_sps`, `save_optimization_result` などがエクスポートされています。

## 7. エラーハンドリング

`interpret` 関数は `ComputationResult` 型（`Success` または `Failure`）を返します。
例外をスローするのではなく、結果型でのハンドリングが推奨されます。

```julia
result = interpret(cmd)
if is_success(result)
    val = unwrap(result)
else
    println("Error: $(result.error)")
end
```
