# Interpreters Layer Documentation

`src/Interpreters/` ディレクトリは、eDSLコマンドを解釈し、実際の数値計算（副作用を伴う計算）を実行するインプリメンテーション層です。
FDDアーキテクチャにおいて、数学的な定義（Domain/Logic）と、それを実行するための計算機処理（MPI, ODEソルバー等）を分離する役割を持ちます。

## ファイル一覧

- `SimulationInterpreter.jl`: 微分方程式の数値積分（SPS等）
- `PhaseInterpreter.jl`: 位相計算パイプライン（PSF, PCF, 最適化）

---

## 1. SimulationInterpreter.jl

`DifferentialEquations.jl` を使用して、力学系のシミュレーションコマンドを実行します。

### SPS計算 (`ComputeSPS`)
安定周期解 (Stable Periodic Solution) を求めます。

1. **緩和**: 初期状態から `t_relax` 時間だけシミュレーションを実行し、過渡状態を脱します。
2. **周期推定**: `estimate_period` 関数でピーク検出を行い、周期 $T$ を概算します。
3. **軌道サンプリング**: 推定された周期 $T$ で再度シミュレーションを行い、$N_\theta$ 点の等時間間隔データを取得します。
4. **補正**: 終点を始点と一致させ、閉軌道として `StablePeriodicSolution` を返します。

### PSF計算 (`ComputePSF`)
随伴方程式 (Adjoint Equation) を解いて位相感受関数 $Z(\theta)$ (コード内 `Q`) を求めます。

- アルゴリズム: 後退差分法による反復解法
    - $\frac{dZ}{d\theta} = -J(X_s(\theta))^T Z(\theta)$
- **`iterate_adjoint_equation`**: 時間を逆行 ($\theta: 2\pi \to 0$) して1周期分積分します。
- **`normalize_psf`**: $Z(\theta) \cdot F(X_s(\theta)) = \omega$ を満たすように正規化します。

### フルダイナミクス (`SimulateFullDynamics`)
高次結合を含む3ネットワーク全体のシミュレーションを実行します。

- 状態ベクトルを連結 ($3 \times N \times D_s$) して ODE を解きます。

### 位相ダイナミクス (`SimulatePhaseDynamics`)
縮約された位相方程式をシミュレーションします。

- 定義された `total_pcf` ($\Gamma$) を用いて $\dot{\phi}$ を計算します。

---

## 2. PhaseInterpreter.jl

位相縮約理論に基づく計算パイプラインを実行します。

### 個別PCF計算 (`ComputeIndividualPCF`)
振動子のトリプレット $(i, j, k)$ ごとの結合関数 $\gamma_{ijk}$ を計算します。
Logic層の `compute_individual_pcf` を呼び出し、結果を返します。

- **`compute_all_individual_pcfs_parallel`**: 全ての組み合わせ ($N^3$個) を計算します。

### 全体PCF計算 (`ComputeTotalPCF`)
結合テンソル $C$ と個別PCF $\gamma_{ijk}$ から全体結合関数 $\Gamma$ を構成します。

### 最適化コマンド
- **`OptimizeLinearStability`**: 線形安定性を最適化するテンソルを探索。
- **`OptimizeRotation`**: 回転特性を最適化するテンソルを探索。
- それぞれ Logic層の最適化関数を呼び出し、結果の検証 (`validate_optimization_result`) と性能評価 (`evaluate_coupling_performance`) を行います。

### パイプライン実行 (`run_full_pipeline`)
一連の解析フローを自動実行する高レベル関数です。

1. **SPS計算**: リミットサイクル取得
2. **PSF計算**: 位相感受関数取得
3. **IndividualPCFs計算**: $N^3$ 個の要素関数計算
4. **最適化**: 指定された戦略 (`:linear`, `:rotation`, `:uniform`) で結合テンソル $C$ を決定
5. **TotalPCF & Metrics**: 最終的な結合関数と安定性を評価

### 実験用 (`run_comparison_experiment`)
一様結合、線形最適化、回転最適化の3条件を同一パラメータで実行し、性能比較を行います。
