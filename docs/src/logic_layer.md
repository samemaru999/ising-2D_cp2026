# Logic Layer Documentation

`src/Logic/` ディレクトリは、FDD (Functional Declarative Design) におけるビジネスロジック、数学的アルゴリズム、力学、および最適化理論を担当します。このレイヤーは純粋 (Pure) であり、データベースやファイルシステムへの副作用を持ちません。

## ファイル一覧

- `FHN.jl`: FitzHugh-Nagumoモデルの微分方程式とヤコビアン
- `Coupling.jl`: 拡散結合および高次結合関数
- `PhaseReduction.jl`: 位相縮約理論、安定性指標、縮約位相方程式
- `Optimization.jl`: 結合テンソルの解析的最適化アルゴリズム

---

## 1. FHN.jl

FitzHugh-Nagumo (FHN) 振動子の力学系を定義します。

### 単一振動子力学
FHNモデル:
- $\dot{u} = \delta(a + v - bu)$
- $\dot{v} = v - v^3/3 - u + I$

- **`fhn_single_dynamics(u, v, I, fhn)`**: 時間微分 $(\dot{u}, \dot{v})$ を計算。
- **`fhn_single_jacobian(u, v, fhn)`**: ヤコビアン行列 $J$ を計算。

### ネットワーク力学
- **`network_intrinsic_dynamics(X, params)`**: 結合なしの個々の振動子の力学。
- **`network_full_dynamics(X, params)`**: ネットワーク内結合 $K_{ij}$ を含む完全な力学系。
- **`network_full_jacobian(X, params)`**: ネットワーク全体のヤコビアン行列（$2N \times 2N$）。

### 3ネットワーク結合系
- **`three_network_dynamics(X_A, X_B, X_C, three_params)`**: 3つのネットワーク間の高次結合を含む時間発展方程式を計算します。
    - 高次相互作用: $H_{ijk}(x_i, x_j, x_k) = (0, \sin(v_j + v_k - 2v_i))^T$

---

## 2. Coupling.jl

結合項の具体的な計算式を提供します。

### ネットワーク内結合 (Diffusive Coupling)
- **`diffusive_coupling(v_i, v_j, K_ij)`**: $K_{ij}(v_j - v_i)$ を計算。
- **`compute_intra_coupling(X, K, N)`**: 全振動子の結合項ベクトルを計算。

### ネットワーク間結合 (High-Order Interaction)
- **`higher_order_H(v_i, v_j, v_k)`**: $\sin(v_j + v_k - 2v_i)$ を計算。
- **`compute_interaction_tensor(X_self, X1, X2, N)`**: 状態ベクトルから相互作用テンソル $H_{ijk}$ を生成。
- **`compute_three_network_coupling`**: 3ネットワーク系における各ネットワークへの結合入力ベクトルを計算。

### 結合行列生成
- `create_full_coupling_matrix`: 全結合行列。
- `create_ring_coupling_matrix`: リング型結合行列。
- `create_uniform_coupling_tensor`: 一様結合テンソル。

---

## 3. PhaseReduction.jl

位相縮約法に基づく理論計算を行います。

### 位相結合関数 (PCF)
- **`compute_individual_pcf(sps, psf, i, j, k)`**:
    - 定義: $\gamma_{ijk}(\phi, \psi) = \frac{1}{T} \int_0^T Z(\theta) \cdot H(X_i(\theta), X_j(\theta-\phi), X_k(\theta-\psi)) d\theta$
    - リミットサイクル $X_s$ と位相感受関数 $Z$ (コード内では `psf.Q`) を用いて計算。
- **`compute_total_pcf(individual_pcfs, C, ...)`**:
    - 全体PCF $\Gamma(\phi, \psi) = \sum C_{ijk} \gamma_{ijk}(\phi, \psi)$ を計算。

### 安定性・回転特性
- **`compute_linear_stability(Γ₁, Γ₂)`**: $\Lambda = -\frac{3}{2}(\Gamma_1 + \Gamma_2)$
    - $\Gamma_1 = \partial_\phi \Gamma(0,0), \quad \Gamma_2 = \partial_\psi \Gamma(0,0)$
- **`compute_rotation_characteristic(Γ₁, Γ₂)`**: $R = \frac{\sqrt{3}}{2}|\Gamma_1 - \Gamma_2|$

### 縮約位相方程式
- **`reduced_phase_dynamics(φ₁, φ₂, pcf, ε)`**:
    - $\dot{\phi}_1 = \epsilon [\Gamma(\phi_1, \phi_2) - \Gamma(\phi_2, \phi_1)]$
    - $\dot{\phi}_2 = \epsilon [\Gamma(\phi_2, -\phi_1-\phi_2) - \Gamma(-\phi_1-\phi_2, \phi_2)]$

---

## 4. Optimization.jl

結合テンソル $C_{ijk}$ の解析的最適化を行います。

### 導関数の計算
- **`compute_all_gamma_derivatives`**: 各 $\gamma_{ijk}$ の原点での偏微分係数 $\gamma_{1,ijk}, \gamma_{2,ijk}$ を計算し、テンソルとして返します。

### 最適化アルゴリズム
- **`optimize_linear_stability`**: 線形安定性 $\Lambda$ を最大化（負の方向に大きく）。
    - 解: $C_{ijk} \propto (\gamma_{1,ijk} + \gamma_{2,ijk})$
- **`optimize_rotation`**: 回転特性 $R$ を最大化。
    - 解: $C_{ijk} \propto (\gamma_{1,ijk} - \gamma_{2,ijk})$
- **`optimize_combined`**: 両者の加重和を最適化。

### 補助関数
- **`evaluate_coupling_performance`**: 与えられた $C$ に対する $\Lambda, R$ を計算。
- **`compare_coupling_strategies`**: 一様結合、線形最適化、回転最適化の結果を比較。
