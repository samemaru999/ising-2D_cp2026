# Domain Layer Documentation

`src/Domain/` ディレクトリは、FDD (Functional Declarative Design) におけるドメインモデル、定数、型定義、およびバリデーションロジックを担当します。このレイヤーは純粋 (Pure) であり、副作用を持ちません。

## ファイル一覧

- `Constants.jl`: 物理定数、計算パラメータ、事前定義された結合行列
- `Types.jl`: ドメインモデル（型定義）
- `DSL.jl`: eDSL (embedded Domain Specific Language) コマンド定義
- `Validation.jl`: パラメータと計算結果のバリデーションロジック

---

## 1. Constants.jl

システム全体で使用される定数とデフォルトパラメータを定義します。

### 主要定数

| 定数名 | 値 | 説明 |
| :--- | :--- | :--- |
| `N_OSC_DEFAULT` | 4 | デフォルトの振動子数 |
| `N_NET` | 3 | ネットワーク数 (A, B, C) |
| `Ds` | 2 | 状態次元 (u, v) |
| `Nθ_DEFAULT` | 101 | 位相分割数 |
| `TREL` | 500.0 | 初期緩和時間 |
| `dt_DEFAULT` | 1e-3 | 時間刻み |
| `ε_DEFAULT` | 2e-4 | 結合強度 |

### FHNパラメータ（デフォルト）

- `δ_DEFAULT`: 0.08 (時間スケール分離)
- `a_DEFAULT`: 0.7
- `b_DEFAULT`: 0.8

### 事前定義結合行列

- `INTRA_K_4`: N=4 用のネットワーク内結合行列
- `N=10` 用の `INTRA_K_10` も利用可能

### 事前定義高次結合テンソル

- `HO_LOPT_C_4`: 線形安定性最適化済みの結合テンソル (N=4)
- `HO_ROPT_C_4`: 回転特性最適化済みの結合テンソル (N=4)

### ヘルパー関数

- `get_intra_K(N)`: 指定されたNに対応する結合行列を返します。
- `get_default_external_input(N)`: ペースメーカーを含む外部入力ベクトルを生成します。
- `create_uniform_coupling_tensor(N)`: 一様な値を持つ結合テンソルを生成します。

---

## 2. Types.jl

アプリケーションで使用される全てのデータ構造を定義します。全ての構造体は不変 (Immutable) です。

### パラメータ型

- **`FHNParams`**: FitzHugh-Nagumoモデルのパラメータ (δ, a, b)。
- **`NetworkParams`**: ネットワーク全体のパラメータ (N, FHNParams, I, K, Topology)。
- **`ThreeNetworkParams`**: 3つのネットワーク間の結合パラメータ (C, ε)。

### 計算結果型

- **`StablePeriodicSolution` (SPS)**: 安定周期解のリミットサイクル軌道などを保持。
- **`PhaseSensitivityFunction` (PSF)**: 位相感受関数 $Z(\theta)$ (コード内では `Q`) を保持。
- **`IndividualPCF`**: 個別の3体位相結合関数 $\gamma_{ijk}(\phi, \psi)$。
- **`TotalPCF`**: 全体位相結合関数 $\Gamma(\phi, \psi)$。
- **`OptimizationResult`**: 最適化計算の結果 (最適化されたテンソル $C_{opt}$ など)。
- **`PhaseTimeSeries`**: 位相ダイナミクスシミュレーションの時系列データ。

### ユーティリティ関数

- `get_oscillator_state(X, i)`: 状態ベクトルから特定の振動子の状態を抽出。
- `index_to_phase(k, Ntheta)`: インデックスを位相 $[0, 2\pi)$ に変換。
- `phase_to_index(theta, Ntheta)`: 位相をインデックスに変換。

---

## 3. DSL.jl

計算処理を表すコマンドを定義します。これらはすべて `ComputationCommand` のサブタイプです。

### コマンド一覧

| コマンド名 | 役割 |
| :--- | :--- |
| **`ComputeSPS`** | 安定周期解 (Stable Periodic Solution) を計算 |
| **`ComputePSF`** | 位相感受関数 (Adjoint method) を計算 |
| **`ComputeIndividualPCF`** | 3体位相結合関数 $\gamma_{ijk}$ を計算 |
| **`ComputeTotalPCF`** | 全体位相結合関数 $\Gamma$ を構成 |
| **`OptimizeLinearStability`** | 線形安定性を最大化する結合テンソルを探索 |
| **`OptimizeRotation`** | 回転特性を最大化する結合テンソルを探索 |
| **`SimulatePhaseDynamics`** | 縮約位相モデルのシミュレーションを実行 |
| **`SimulateFullDynamics`** | 元の微分方程式によるシミュレーションを実行 |
| **`SaveData` / `LoadData`** | データの保存・読み込み |

### 結果型 (`ComputationResult`)

計算結果は `Success` または `Failure` でラップされます。
`unwrap` 関数で値を取り出すか、`is_success` で成功判定を行います。

---

## 4. Validation.jl

データの整合性をチェックするための関数群です。

### パラメータ検証

- `validate_fhn_params(delta, a, b)`
- `validate_network_params(N, I, K)`
- `validate_coupling_tensor(C, N)`

### 結果検証

- **`validate_sps`**: 周期性 ($X(0) \approx X(T)$)、角周波数の整合性などをチェック。
- **`validate_psf`**: 正規化条件 ($Q \cdot F = \omega$)、周期性をチェック。
- **`validate_optimization_result`**: 目標値の達成度、不適切な値 (NaN/Inf) の有無をチェック。

これらの関数は、バリデーションエラーのリスト (`Vector{ValidationError}`) または `Nothing` (問題なし) を返します。
