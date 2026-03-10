# 実装タスクチェックリスト

**プロジェクト名**: 2次元イジングモデル モンテカルロシミュレーション (Ising2D)
**バージョン**: 2.0
**作成日**: 2025-12-05
**最終更新日**: 2025-12-05
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

### FDDフェーズと推定工数

| FDD Phase | 名称 | 内容 | 推定工数 |
|-----------|------|------|----------|
| Phase 1 | 要件分析 | ドメイン理解、要件モデル構築 | 2日 |
| Phase 2 | アーキテクチャ設計 | レイヤー構造、サブシステム設計 | 2日 |
| Phase 3 | ドメインモデル設計 | 型定義、eDSL設計 | 3日 |
| Phase 4 | インターフェース設計 | サービス境界、依存性注入 | 2日 |
| Phase 5 | 実装 | Pure/Impure関数の実装 | 12日 |
| Phase 6 | テスト | ユニット/統合テスト | 4日 |

**合計推定工数**: 約25日（約1ヶ月）

### 機能要件との対応

| 機能要件ID | 機能名 | FDD実装Phase | 優先度 |
|-----------|--------|-------------|--------|
| F-SIM | 元力学系シミュレーション | Phase 5 | 最高 |
| F-PSF | 位相感受関数導出 | Phase 5 | 高 |
| F-OPT | 最適結合テンソル導出 | Phase 5 | 高 |
| F-PCF | 位相結合関数導出 | Phase 5 | 高 |
| F-PHASE | 位相縮約後シミュレーション | Phase 5 | 中 |
| F-VIS | 可視化 | Phase 5 | 中 |

---

## FDD Phase 1: 要件分析

### 1.1 マインドマップ分析

- [ ] **TASK-1.1.1**: ドメイン概念の抽出
  - 複雑度: 低
  - 依存: なし
  - 成果物: `docs/requirements/mind_maps/domain_concepts.md`
  - 受け入れ条件:
    - 主要概念（FHN振動子、ネットワーク、位相縮約等）が列挙されている
    - 概念間の関連が明確化されている

- [ ] **TASK-1.1.2**: 数学モデルの整理
  - 複雑度: 中
  - 依存: TASK-1.1.1
  - 成果物: `docs/requirements/mind_maps/mathematical_model.md`
  - 受け入れ条件:
    - FHN方程式が明記されている
    - 結合関数（ペアワイズ、高次）が定義されている
    - 位相縮約の数式が整理されている

### 1.2 ユーザーシナリオ作成

- [ ] **TASK-1.2.1**: 基本シナリオの作成
  - 複雑度: 低
  - 依存: TASK-1.1.1
  - 成果物: `docs/requirements/scenarios/basic_simulation.md`
  - 受け入れ条件:
    - 安定周期解探索のシナリオが定義されている
    - 位相結合関数計算のシナリオが定義されている

- [ ] **TASK-1.2.2**: 最適化シナリオの作成
  - 複雑度: 低
  - 依存: TASK-1.2.1
  - 成果物: `docs/requirements/scenarios/optimization.md`
  - 受け入れ条件:
    - 線形安定性最適化のシナリオが定義されている
    - 回転特性最適化のシナリオが定義されている

### 1.3 Q&A形式の要件整理

- [ ] **TASK-1.3.1**: 未決定事項の整理
  - 複雑度: 低
  - 依存: TASK-1.1.2
  - 成果物: `docs/requirements/qa.md`
  - 受け入れ条件:
    - N=10での結合行列Kの具体値が確認されている
    - 外部入力I_iの設定が明確化されている
    - 並列計算方針（Distributed vs Threads）が決定されている

---

## FDD Phase 2: アーキテクチャ設計

### 2.1 レイヤー構造設計

- [ ] **TASK-2.1.1**: FDDレイヤー定義
  - 複雑度: 中
  - 依存: Phase 1完了
  - 成果物: `docs/architecture/layers.md`
  - 受け入れ条件:
    - 6層構造（App, Logic, Service, Domain, Interpreters, Runtime）が定義されている
    - 各レイヤーの責務が明確化されている
    - Pure/Impure分離が明記されている

- [ ] **TASK-2.1.2**: レイヤー間依存関係の定義
  - 複雑度: 中
  - 依存: TASK-2.1.1
  - 成果物: `docs/architecture/dependencies.md`
  - 受け入れ条件:
    - 依存方向（Pure -> Impure）が明記されている
    - 禁止される依存が列挙されている

### 2.2 サブシステム設計

- [ ] **TASK-2.2.1**: コアサブシステムの特定
  - 複雑度: 中
  - 依存: TASK-2.1.1
  - 成果物: `docs/architecture/subsystems.md`
  - 受け入れ条件:
    - Simulation, PSF, PCF, Optimization, Phase Dynamicsサブシステムが定義されている
    - 各サブシステムのインターフェースが概略定義されている

---

## FDD Phase 3: ドメインモデル設計

### 3.1 ドメイン型定義 (Domain Layer - Pure)

- [ ] **TASK-3.1.1**: 基本型の設計
  - 複雑度: 中
  - 依存: Phase 2完了
  - 成果物: `src/Domain/Types.jl` の設計仕様
  - 受け入れ条件:
    - FHNState, FHNParams 型が設計されている
    - OscillatorState = SVector{2, Float64} が定義されている
    - NetworkState 型が設計されている

- [ ] **TASK-3.1.2**: 結合構造型の設計
  - 複雑度: 中
  - 依存: TASK-3.1.1
  - 成果物: `src/Domain/Types.jl` の設計仕様
  - 受け入れ条件:
    - CouplingMatrix 型が設計されている
    - CouplingTensor 型が設計されている
    - NetworkTopology 型が設計されている

- [ ] **TASK-3.1.3**: 計算結果型の設計
  - 複雑度: 中
  - 依存: TASK-3.1.1
  - 成果物: `src/Domain/Types.jl` の設計仕様
  - 受け入れ条件:
    - StablePeriodicSolution 型が設計されている
    - PhaseSensitivityFunction 型が設計されている
    - PhaseCouplingFunction 型が設計されている
    - OptimizedCouplingTensor 型が設計されている

### 3.2 eDSL設計 (Domain Layer - Pure)

- [ ] **TASK-3.2.1**: シミュレーションコマンドの設計
  - 複雑度: 中
  - 依存: TASK-3.1.3
  - 成果物: `src/Domain/DSL.jl` の設計仕様
  - 受け入れ条件:
    - abstract type SimCommand が設計されている
    - FindPeriodicSolution, ComputePSF, ComputePCF 等のコマンドが設計されている
    - Scenario = Vector{SimCommand} が定義されている

- [ ] **TASK-3.2.2**: スマートコンストラクタの設計
  - 複雑度: 低
  - 依存: TASK-3.1.2
  - 成果物: `src/Domain/Constructors.jl` の設計仕様
  - 受け入れ条件:
    - 検証付きコンストラクタが設計されている
    - 不正状態の生成が防止される設計になっている

### 3.3 検証ロジック設計 (Domain Layer - Pure)

- [ ] **TASK-3.3.1**: パラメータ検証の設計
  - 複雑度: 低
  - 依存: TASK-3.1.1
  - 成果物: `src/Domain/Validation.jl` の設計仕様
  - 受け入れ条件:
    - validate_fhn_params が設計されている
    - validate_coupling_strength (0 < epsilon << 1) が設計されている
    - validate_topology が設計されている

---

## FDD Phase 4: インターフェース設計

### 4.1 サービスインターフェース設計 (Service Layer)

- [ ] **TASK-4.1.1**: SimulationService インターフェース
  - 複雑度: 中
  - 依存: Phase 3完了
  - 成果物: `src/Services/Interfaces.jl` の設計仕様
  - 受け入れ条件:
    - find_periodic_solution 関数シグネチャが定義されている
    - simulate_network 関数シグネチャが定義されている
    - assign_phase 関数シグネチャが定義されている

- [ ] **TASK-4.1.2**: PhaseReductionService インターフェース
  - 複雑度: 中
  - 依存: TASK-4.1.1
  - 成果物: `src/Services/Interfaces.jl` の設計仕様
  - 受け入れ条件:
    - compute_psf 関数シグネチャが定義されている
    - compute_individual_pcf 関数シグネチャが定義されている
    - compute_total_pcf 関数シグネチャが定義されている

- [ ] **TASK-4.1.3**: OptimizationService インターフェース
  - 複雑度: 中
  - 依存: TASK-4.1.2
  - 成果物: `src/Services/Interfaces.jl` の設計仕様
  - 受け入れ条件:
    - optimize_linear_stability 関数シグネチャが定義されている
    - optimize_rotation 関数シグネチャが定義されている

### 4.2 Service Handleパターン設計

- [ ] **TASK-4.2.1**: サービスハンドル構造体の設計
  - 複雑度: 中
  - 依存: TASK-4.1.3
  - 成果物: `src/Services/ServiceHandle.jl` の設計仕様
  - 受け入れ条件:
    - HONetSyncServices 構造体が設計されている
    - create_production_services 関数が設計されている
    - create_mock_services 関数が設計されている

---

## FDD Phase 5: 実装

### 5.1 基盤実装

#### 5.1.1 プロジェクト初期化

- [ ] **TASK-5.1.1**: Project.toml の作成
  - 複雑度: 低
  - 依存: なし
  - 成果物: `Project.toml`
  - 受け入れ条件: `Pkg.instantiate()` で依存関係がインストールされる

- [ ] **TASK-5.1.2**: ディレクトリ構造の作成
  - 複雑度: 低
  - 依存: なし
  - 成果物: FDD準拠のディレクトリ構造
  - 受け入れ条件: `src/{Domain,Services,Logic,Interpreters,Runtime,App}` が存在する

- [ ] **TASK-5.1.3**: メインモジュール HONetSync.jl の作成
  - 複雑度: 低
  - 依存: TASK-5.1.2
  - 成果物: `src/HONetSync.jl`
  - 受け入れ条件: `using HONetSync` でモジュールが読み込める

#### 5.1.2 定数・パラメータ定義

- [ ] **TASK-5.1.4**: constants.jl の実装
  - 複雑度: 低
  - 依存: TASK-5.1.2
  - 成果物: `src/utils/constants.jl`
  - 受け入れ条件:
    - N_osc, N_net, Ntheta, Ds 等の定数が定義されている
    - N=4 と N=10 の切り替えが容易

- [ ] **TASK-5.1.5**: 事前定義結合行列の実装
  - 複雑度: 低
  - 依存: TASK-5.1.4
  - 成果物: `src/utils/constants.jl`
  - 受け入れ条件:
    - intraK_4, intraK_10 が定義されている
    - ho_interC_uni_4, ho_interC_opt_4 等が定義されている

### 5.2 Domain Layer 実装 (Pure)

#### 5.2.1 型定義

- [ ] **TASK-5.2.1**: Types.jl の実装
  - 複雑度: 中
  - 依存: TASK-5.1.4, Phase 3設計
  - 成果物: `src/Domain/Types.jl`
  - Pure/Impure: **Pure**
  - 受け入れ条件:
    - FHNState, FHNParams 型が実装されている
    - NetworkTopology 型が実装されている
    - 全計算結果型が実装されている

#### 5.2.2 数学モデル (Pure関数)

- [ ] **TASK-5.2.2**: FHNモデルの実装
  - 複雑度: 中
  - 依存: TASK-5.2.1
  - 成果物: `src/Domain/Models.jl`
  - Pure/Impure: **Pure**
  - 対応要件: F-SIM-01
  - 受け入れ条件:
    - `fhn_dynamics(state, params)` がPure関数として実装されている
    - du/dt = delta*(a + v - b*u) が正しく計算される
    - dv/dt = v - v^3/3 - u + I が正しく計算される

- [ ] **TASK-5.2.3**: ネットワーク内結合関数の実装
  - 複雑度: 中
  - 依存: TASK-5.2.2
  - 成果物: `src/Domain/Models.jl`
  - Pure/Impure: **Pure**
  - 対応要件: F-SIM-02
  - 受け入れ条件:
    - `intra_coupling(x_i, x_j, K_ij)` がPure関数として実装されている
    - g_ij = (0, K_ij*(v_j - v_i)) が計算される

- [ ] **TASK-5.2.4**: 高次相互作用関数の実装
  - 複雑度: 中
  - 依存: TASK-5.2.2
  - 成果物: `src/Domain/Models.jl`
  - Pure/Impure: **Pure**
  - 対応要件: F-SIM-03
  - 受け入れ条件:
    - `higher_order_interaction(x_i, x_j, x_k)` がPure関数として実装されている
    - H_ijk = (0, sin(v_j + v_k - 2*v_i)) が計算される

- [ ] **TASK-5.2.5**: ヤコビ行列の解析的実装
  - 複雑度: 高
  - 依存: TASK-5.2.3
  - 成果物: `src/Domain/Models.jl`
  - Pure/Impure: **Pure**
  - 対応要件: F-PSF-01
  - 受け入れ条件:
    - `jacobian_fhn(state, params)` がPure関数として実装されている
    - `jacobian_network(states, params, K)` が実装されている
    - 自動微分を使用していない

#### 5.2.3 eDSL実装

- [ ] **TASK-5.2.6**: DSL.jl の実装
  - 複雑度: 中
  - 依存: TASK-5.2.1
  - 成果物: `src/Domain/DSL.jl`
  - Pure/Impure: **Pure**
  - 受け入れ条件:
    - SimCommand 抽象型と具象コマンドが実装されている
    - コマンドは不変のデータ構造として表現されている

#### 5.2.4 検証ロジック

- [ ] **TASK-5.2.7**: Validation.jl の実装
  - 複雑度: 低
  - 依存: TASK-5.2.1
  - 成果物: `src/Domain/Validation.jl`
  - Pure/Impure: **Pure**
  - 受け入れ条件:
    - `validate_params` がPure関数として実装されている
    - `validate_coupling_strength` が 0 < epsilon << 1 を検証する

### 5.3 Services Layer 実装 (Mixed)

- [ ] **TASK-5.3.1**: Interfaces.jl の実装
  - 複雑度: 中
  - 依存: Phase 4設計
  - 成果物: `src/Services/Interfaces.jl`
  - Pure/Impure: **Pure** (インターフェース定義のみ)
  - 受け入れ条件:
    - サービスハンドル構造体が実装されている
    - 関数シグネチャが型アノテーション付きで定義されている

### 5.4 Logic Layer 実装 (Pure)

- [ ] **TASK-5.4.1**: シナリオビルダーの実装
  - 複雑度: 中
  - 依存: TASK-5.2.6, TASK-5.3.1
  - 成果物: `src/Logic/Scenarios.jl`
  - Pure/Impure: **Pure**
  - 受け入れ条件:
    - `build_full_analysis_scenario(params)` がPure関数として実装されている
    - シナリオはコマンドのリストとして表現される

- [ ] **TASK-5.4.2**: 最適化ロジックの実装
  - 複雑度: 高
  - 依存: TASK-5.2.5
  - 成果物: `src/Logic/Optimization.jl`
  - Pure/Impure: **Pure**
  - 対応要件: F-OPT-02, F-OPT-03
  - 受け入れ条件:
    - `compute_optimal_linear_coupling(V1, V2, q)` がPure関数として実装されている
    - `compute_optimal_rotation_coupling(V1, V2, mu)` がPure関数として実装されている
    - 解析解が正しく計算される

### 5.5 Interpreters Layer 実装 (Impure)

#### 5.5.1 ODE インタープリター

- [ ] **TASK-5.5.1**: ODEインタープリターの実装
  - 複雑度: 高
  - 依存: TASK-5.2.4
  - 成果物: `src/Interpreters/ODEInterpreter.jl`
  - Pure/Impure: **Impure**
  - 受け入れ条件:
    - DifferentialEquations.jl を使用したODE求解が実装されている
    - 適応的ステップサイズ制御が有効

- [ ] **TASK-5.5.2**: 安定周期解探索の実装
  - 複雑度: 高
  - 依存: TASK-5.5.1
  - 成果物: `src/Interpreters/ODEInterpreter.jl`
  - Pure/Impure: **Impure**
  - 対応要件: F-SIM-04
  - 受け入れ条件:
    - 初期緩和処理が実装されている
    - 位相原点検出が実装されている
    - 周期T、角周波数omegaが正しく計算される
    - Xs が Ntheta点で離散化される

#### 5.5.2 位相縮約インタープリター

- [ ] **TASK-5.5.3**: 随伴方程式積分の実装
  - 複雑度: 高
  - 依存: TASK-5.2.5, TASK-5.5.2
  - 成果物: `src/Interpreters/PhaseInterpreter.jl`
  - Pure/Impure: **Impure**
  - 対応要件: F-PSF-01
  - 受け入れ条件:
    - dQ/dt = -J^T * Q の逆向き積分が実装されている
    - 正規化条件 Q*F = omega が適用される
    - REP回の反復で収束する

- [ ] **TASK-5.5.4**: 個別位相結合関数計算の実装
  - 複雑度: 高
  - 依存: TASK-5.5.3
  - 成果物: `src/Interpreters/PhaseInterpreter.jl`
  - Pure/Impure: **Impure**
  - 対応要件: F-PCF-01
  - 受け入れ条件:
    - gamma_ijk(phi, psi) の数値積分が実装されている
    - 並列計算（@distributed）が適用される

- [ ] **TASK-5.5.5**: 位相割り当ての実装
  - 複雑度: 中
  - 依存: TASK-5.5.2
  - 成果物: `src/Interpreters/PhaseInterpreter.jl`
  - Pure/Impure: **Pure** (計算自体は純粋)
  - 対応要件: F-SIM-05
  - 受け入れ条件:
    - 最近傍法が実装されている
    - ゼロクロス法が実装されている
    - 3ネットワーク間位相差が計算される

#### 5.5.3 シナリオインタープリター

- [ ] **TASK-5.5.6**: SimInterpreter.jl の実装
  - 複雑度: 高
  - 依存: TASK-5.5.4, TASK-5.2.6
  - 成果物: `src/Interpreters/SimInterpreter.jl`
  - Pure/Impure: **Impure**
  - 受け入れ条件:
    - `interpret(cmd::SimCommand, ctx)` が多重ディスパッチで実装されている
    - 各コマンド型に対するインタープリターが実装されている

### 5.6 Runtime Layer 実装 (Impure)

- [ ] **TASK-5.6.1**: IOHandler.jl の実装
  - 複雑度: 中
  - 依存: TASK-5.2.1
  - 成果物: `src/Runtime/IOHandler.jl`
  - Pure/Impure: **Impure**
  - 対応要件: NFR-REPR-02
  - 受け入れ条件:
    - JLD2形式での保存が実装されている
    - JLD2形式での読み込みが実装されている
    - パラメータも含めて保存される

- [ ] **TASK-5.6.2**: StateManager.jl の実装
  - 複雑度: 中
  - 依存: TASK-5.6.1
  - 成果物: `src/Runtime/StateManager.jl`
  - Pure/Impure: **Impure**
  - 受け入れ条件:
    - RuntimeState 構造体が実装されている
    - 状態更新関数が実装されている

### 5.7 App Layer 実装 (Impure)

- [ ] **TASK-5.7.1**: Main.jl の実装
  - 複雑度: 中
  - 依存: TASK-5.5.6, TASK-5.6.2
  - 成果物: `src/App/Main.jl`
  - Pure/Impure: **Impure**
  - 受け入れ条件:
    - アプリケーションエントリーポイントが実装されている
    - 設定読み込みが実装されている

- [ ] **TASK-5.7.2**: 可視化関数の実装
  - 複雑度: 中
  - 依存: TASK-5.2.1
  - 成果物: `src/visualization/plotting.jl`
  - Pure/Impure: **Impure**
  - 対応要件: F-VIS-01〜05
  - 受け入れ条件:
    - リミットサイクルプロットが実装されている
    - 位相感受関数プロットが実装されている
    - 位相結合関数の3D/ヒートマップが実装されている
    - 位相空間軌道プロットが実装されている

---

## FDD Phase 6: テスト

### 6.1 ユニットテスト (Pure関数)

- [ ] **TASK-6.1.1**: Domain型のテスト
  - 複雑度: 低
  - 依存: TASK-5.2.1
  - 成果物: `test/Unit/test_types.jl`
  - 受け入れ条件:
    - 型の構築が正しく動作する
    - 不正な値でエラーが発生する

- [ ] **TASK-6.1.2**: Pure数学関数のテスト
  - 複雑度: 中
  - 依存: TASK-5.2.4
  - 成果物: `test/Unit/test_models.jl`
  - 受け入れ条件:
    - fhn_dynamics のテストがある
    - higher_order_interaction のテストがある
    - jacobian_fhn のテストがある

- [ ] **TASK-6.1.3**: 検証ロジックのテスト
  - 複雑度: 低
  - 依存: TASK-5.2.7
  - 成果物: `test/Unit/test_validation.jl`
  - 受け入れ条件:
    - 有効なパラメータで成功する
    - 無効なパラメータでエラーを返す

- [ ] **TASK-6.1.4**: 最適化ロジックのテスト
  - 複雑度: 中
  - 依存: TASK-5.4.2
  - 成果物: `test/Unit/test_optimization.jl`
  - 受け入れ条件:
    - 最適化後の Lambda が目標値 q に一致する
    - 最適化後の R が目標値 mu に一致する

### 6.2 インテグレーションテスト (モック使用)

- [ ] **TASK-6.2.1**: サービスモックの実装
  - 複雑度: 中
  - 依存: TASK-5.3.1
  - 成果物: `test/Integration/mocks.jl`
  - 受け入れ条件:
    - MockSimulationService が実装されている
    - MockPSFService が実装されている

- [ ] **TASK-6.2.2**: シナリオ実行テスト
  - 複雑度: 中
  - 依存: TASK-6.2.1, TASK-5.5.6
  - 成果物: `test/Integration/test_scenarios.jl`
  - 受け入れ条件:
    - モックを使用したシナリオ実行が成功する
    - 各ステップの結果が検証される

### 6.3 N=4 テストケース

- [ ] **TASK-6.3.1**: N=4での安定周期解探索テスト
  - 複雑度: 中
  - 依存: TASK-5.5.2
  - 成果物: `test/E2E/test_n4_sps.jl`
  - 対応要件: NFR-MAINT-03
  - 受け入れ条件:
    - N=4で数秒以内に計算完了
    - 周期性が確認できる

- [ ] **TASK-6.3.2**: N=4での位相感受関数テスト
  - 複雑度: 中
  - 依存: TASK-5.5.3
  - 成果物: `test/E2E/test_n4_psf.jl`
  - 対応要件: NFR-MAINT-03
  - 受け入れ条件:
    - N=4で数十秒以内に計算完了
    - 正規化条件 Q*F = omega が満たされている

- [ ] **TASK-6.3.3**: N=4での位相結合関数テスト
  - 複雑度: 中
  - 依存: TASK-5.5.4
  - 成果物: `test/E2E/test_n4_pcf.jl`
  - 対応要件: NFR-MAINT-03
  - 受け入れ条件:
    - N=4で数分以内に計算完了
    - 既存データとの整合性確認

- [ ] **TASK-6.3.4**: 既存データとの比較検証
  - 複雑度: 中
  - 依存: TASK-6.3.3
  - 成果物: `test/E2E/test_comparison.jl`
  - 受け入れ条件:
    - ho_interC_uni_4, ho_interC_opt_4 との整合性確認
    - ho_loptC_4, ho_roptC_4 との比較

### 6.4 N=10 本番テスト

- [ ] **TASK-6.4.1**: N=10での全パイプライン実行
  - 複雑度: 高
  - 依存: Phase 5全て
  - 成果物: `test/E2E/test_n10_full.jl`
  - 対応要件: NFR-PERF-01
  - 受け入れ条件:
    - 安定周期解: 5分以内
    - 位相感受関数: 10分以内
    - 位相結合関数: 1時間以内

- [ ] **TASK-6.4.2**: メモリ使用量の確認
  - 複雑度: 中
  - 依存: TASK-6.4.1
  - 成果物: `test/E2E/test_memory.jl`
  - 対応要件: NFR-PERF-02
  - 受け入れ条件: 最大メモリ使用量16GB以内

### 6.5 数値精度検証

- [ ] **TASK-6.5.1**: 周期解の周期性テスト
  - 複雑度: 低
  - 依存: TASK-5.5.2
  - 成果物: `test/Property/test_periodicity.jl`
  - 受け入れ条件: ||Xs[:,1] - Xs[:,Ntheta]|| < tolerance

- [ ] **TASK-6.5.2**: 正規化条件のテスト
  - 複雑度: 低
  - 依存: TASK-5.5.3
  - 成果物: `test/Property/test_normalization.jl`
  - 受け入れ条件: |Q*F - omega| < tolerance

---

## 機能要件別タスク

### F-SIM: 元力学系シミュレーション

| タスクID | タスク名 | 対応要件 | Pure/Impure |
|---------|---------|---------|------------|
| TASK-5.2.2 | FHNモデルの実装 | F-SIM-01 | Pure |
| TASK-5.2.3 | ネットワーク内結合関数 | F-SIM-02 | Pure |
| TASK-5.2.4 | 高次相互作用関数 | F-SIM-03 | Pure |
| TASK-5.5.2 | 安定周期解探索 | F-SIM-04 | Impure |
| TASK-5.5.5 | 位相割り当て | F-SIM-05 | Pure |

### F-PSF: 位相感受関数導出

| タスクID | タスク名 | 対応要件 | Pure/Impure |
|---------|---------|---------|------------|
| TASK-5.2.5 | ヤコビ行列の解析的実装 | F-PSF-01 | Pure |
| TASK-5.5.3 | 随伴方程式積分 | F-PSF-01 | Impure |

### F-OPT: 最適結合テンソル導出

| タスクID | タスク名 | 対応要件 | Pure/Impure |
|---------|---------|---------|------------|
| TASK-5.4.2 | 最適化ロジック (線形安定性) | F-OPT-02 | Pure |
| TASK-5.4.2 | 最適化ロジック (回転特性) | F-OPT-03 | Pure |

### F-PCF: 位相結合関数導出

| タスクID | タスク名 | 対応要件 | Pure/Impure |
|---------|---------|---------|------------|
| TASK-5.5.4 | 個別位相結合関数計算 | F-PCF-01 | Impure |
| TASK-5.4.1 | 全体位相結合関数計算 | F-PCF-02 | Pure |

### F-PHASE: 位相縮約後シミュレーション

| タスクID | タスク名 | 対応要件 | Pure/Impure |
|---------|---------|---------|------------|
| TASK-5.5.1 | 位相差方程式の積分 | F-PHASE-01 | Impure |

### F-VIS: 可視化

| タスクID | タスク名 | 対応要件 | Pure/Impure |
|---------|---------|---------|------------|
| TASK-5.7.2 | 可視化関数 | F-VIS-01〜05 | Impure |

---

## 非機能要件タスク

### NFR-PERF: 性能要件

| タスクID | タスク名 | 対応要件 | 受け入れ条件 |
|---------|---------|---------|-------------|
| TASK-5.5.4 | 並列計算対応 | NFR-PERF-01 | 1時間以内（N=10） |
| TASK-5.2.1 | StaticArrays使用 | NFR-PERF-01 | 高速な配列演算 |

### NFR-EXT: 拡張性要件

| タスクID | タスク名 | 対応要件 | 受け入れ条件 |
|---------|---------|---------|-------------|
| TASK-5.1.4 | N可変対応 | NFR-EXT-01 | N=4, 10, 20対応 |
| TASK-5.2.2 | 多重ディスパッチ設計 | NFR-EXT-02 | モデル差し替え可能 |

### NFR-MAINT: 保守性要件

| タスクID | タスク名 | 対応要件 | 受け入れ条件 |
|---------|---------|---------|-------------|
| TASK-5.2.* | 型アノテーション | NFR-MAINT-01 | 全関数に型付き |
| TASK-6.3.* | N=4テストケース | NFR-MAINT-03 | 全機能動作確認 |

### NFR-REPR: 再現性要件

| タスクID | タスク名 | 対応要件 | 受け入れ条件 |
|---------|---------|---------|-------------|
| TASK-5.6.1 | JLD2保存 | NFR-REPR-02 | パラメータ含め保存 |

---

## 優先度と依存関係

### 依存関係図

```
FDD Phase 1 (要件分析)
    |
    v
FDD Phase 2 (アーキテクチャ設計)
    |
    v
FDD Phase 3 (ドメインモデル設計)
    |
    v
FDD Phase 4 (インターフェース設計)
    |
    v
FDD Phase 5 (実装)
    |
    +-- 5.1 基盤
    |     |
    |     v
    +-- 5.2 Domain Layer (Pure) ----+
    |     |                          |
    |     v                          |
    +-- 5.3 Services Layer           |
    |     |                          |
    |     v                          |
    +-- 5.4 Logic Layer (Pure) <-----+
    |     |
    |     v
    +-- 5.5 Interpreters Layer (Impure)
    |     |
    |     v
    +-- 5.6 Runtime Layer (Impure)
    |     |
    |     v
    +-- 5.7 App Layer (Impure)
          |
          v
FDD Phase 6 (テスト)
    |
    +-- 6.1 Unit (Pure関数)
    +-- 6.2 Integration (モック)
    +-- 6.3 N=4 E2E
    +-- 6.4 N=10 E2E
    +-- 6.5 Property
```

### クリティカルパス

最も時間がかかる経路:

```
TASK-5.2.1 (Types)
    -> TASK-5.2.4 (高次相互作用)
    -> TASK-5.2.5 (ヤコビ行列)
    -> TASK-5.5.2 (安定周期解探索)
    -> TASK-5.5.3 (位相感受関数)
    -> TASK-5.5.4 (位相結合関数) [最長]
    -> TASK-5.4.2 (最適化)
    -> TASK-6.4.1 (N=10テスト)
```

### 優先度マトリクス

| フェーズ/タスク群 | 重要度 | 緊急度 | 優先度 |
|-----------------|--------|--------|--------|
| Phase 1-2 (設計) | 高 | 高 | 最優先 |
| Phase 3 (ドメインモデル) | 高 | 高 | 最優先 |
| Phase 4 (インターフェース) | 中 | 高 | 高 |
| 5.1-5.2 (基盤・Domain) | 高 | 高 | 最優先 |
| 5.3-5.4 (Services・Logic) | 中 | 中 | 中 |
| 5.5 (Interpreters) | 高 | 中 | 高 |
| 5.6-5.7 (Runtime・App) | 低 | 低 | 後回し可 |
| Phase 6 (テスト) | 中 | 低 | 継続的 |

---

## 付録

### A. タスク複雑度の定義

| 複雑度 | 推定時間 | 説明 |
|--------|----------|------|
| 低 | 1-2時間 | 単純な実装、既存パターンの適用 |
| 中 | 半日-1日 | 設計判断が必要、テストが必要 |
| 高 | 1-3日 | 複雑なアルゴリズム、デバッグ時間含む |

### B. Pure/Impure分類基準

| 分類 | 基準 | 例 |
|------|------|-----|
| **Pure** | 副作用なし、同じ入力→同じ出力 | fhn_dynamics, jacobian_fhn |
| **Impure** | 副作用あり（IO, 状態変更, 乱数） | solve_ode, save_result |

### C. 関連ドキュメント

- [要件定義書](./requirement.md)
- [システムアーキテクチャ設計書](./system_architecture.md)
- [技術スタック仕様書](./tech_stack.md)
- [データ構造定義書](./data_structure.md)
- [ディレクトリ構造](./directory_structure.md)

---

*本文書はrequirement.mdおよびCLAUDE.mdに基づいてFDD Phase 1-6に沿って作成されました。*
