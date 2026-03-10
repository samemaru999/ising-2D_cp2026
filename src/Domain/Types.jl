"""
    Domain/Types.jl - ドメイン型定義 [Pure]

FDDの代数的データ型（ADT）によるドメインモデル
全て不変（immutable）構造体として定義
"""

using StaticArrays
using LinearAlgebra

# =============================================================================
# 基本型エイリアス
# =============================================================================

"""単一振動子の状態（2次元ベクトル: u, v）"""
const OscillatorState = SVector{2,Float64}

"""ネットワーク状態のベクトル表現（length = Ds * N）"""
const NetworkStateVector = Vector{Float64}

# =============================================================================
# 結合タイプ型（直和型 - 多重ディスパッチ用）
# =============================================================================

"""
    CouplingType

結合タイプの抽象型

FDD: 結合タイプを型で表現し、多重ディスパッチで計算方法を切り替える
これにより、拡散結合の場合は簡略化されたヤコビアン計算を使用可能
"""
abstract type CouplingType end

"""
    DiffusiveCoupling <: CouplingType

拡散結合（線形結合）

g_ij(x_i, x_j) = K_ij(v_j - v_i)

この結合タイプでは M_ij と N_ij が結合行列 L に統合され、
PSF.m スタイルの簡略化が適用される:
    J = diag(1 - u²) + L
"""
struct DiffusiveCoupling <: CouplingType
    variable::Symbol  # :u または :v（どの変数に結合が作用するか）

    function DiffusiveCoupling(variable::Symbol=:v)
        variable in (:u, :v) || throw(ArgumentError("variable must be :u or :v"))
        new(variable)
    end
end

"""デフォルト: v成分への拡散結合"""


"""
    GeneralCoupling <: CouplingType

一般的な非線形結合

M_ij, N_ij を明示的に計算する必要がある場合
"""
struct GeneralCoupling <: CouplingType end

# =============================================================================
# トポロジー型（直和型）
# =============================================================================

"""ネットワークトポロジーの抽象型"""
abstract type TopologyType end

"""全結合トポロジー"""
struct FullCoupling <: TopologyType end

"""ツリー型トポロジー"""
struct TreeCoupling <: TopologyType
    parent_indices::Vector{Int}

    function TreeCoupling(parent_indices::Vector{Int})
        # バリデーション: 親インデックスの妥当性
        for (i, p) in enumerate(parent_indices)
            if p != 0 && (p < 0 || p > length(parent_indices))
                throw(DomainError(p, "Invalid parent index at position $i"))
            end
            if p == i
                throw(DomainError(p, "Self-reference at position $i"))
            end
        end
        new(parent_indices)
    end
end

# =============================================================================
# パラメータ構造体（Smart Constructors付き）
# =============================================================================

"""
    FHNParams

FitzHugh-Nagumo振動子のパラメータ

# Fields
- `delta::Float64`: 時間スケール分離パラメータ (δ > 0)
- `a::Float64`: FHNパラメータ a
- `b::Float64`: FHNパラメータ b (0 < b < 1)
"""
struct FHNParams
    delta::Float64
    a::Float64
    b::Float64

    function FHNParams(delta::Float64, a::Float64, b::Float64)
        delta > 0 || throw(DomainError(delta, "delta must be positive"))
        0 < b < 1 || throw(DomainError(b, "b must be in (0, 1)"))
        new(delta, a, b)
    end
end

"""デフォルトFHNパラメータを作成"""
FHNParams() = FHNParams(δ_DEFAULT, a_DEFAULT, b_DEFAULT)

"""
    NetworkParams

ネットワークシミュレーションのパラメータ

# Fields
- `N::Int`: 振動子数
- `fhn::FHNParams`: FHN振動子パラメータ
- `I::Vector{Float64}`: 外部入力ベクトル (長さN)
- `K::Matrix{Float64}`: ネットワーク内結合行列 (N×N)
- `topology::TopologyType`: トポロジータイプ
"""
struct NetworkParams
    N::Int
    fhn::FHNParams
    I::Vector{Float64}
    K::Matrix{Float64}
    topology::TopologyType

    function NetworkParams(N::Int, fhn::FHNParams, I::Vector{Float64},
        K::Matrix{Float64}, topology::TopologyType)
        N > 0 || throw(DomainError(N, "N must be positive"))
        length(I) == N || throw(DimensionMismatch("I must have length N=$N"))
        size(K) == (N, N) || throw(DimensionMismatch("K must be $N x $N"))
        all(diag(K) .== 0) || @warn "Diagonal of K should be zero (self-coupling)"
        new(N, fhn, I, K, topology)
    end
end

"""デフォルトネットワークパラメータを作成"""
function NetworkParams(N::Int; topology::TopologyType=FullCoupling())
    fhn = FHNParams()
    I = get_default_external_input(N)
    K = get_intra_K(N)
    NetworkParams(N, fhn, I, K, topology)
end

"""
    ThreeNetworkParams

3ネットワーク間高次結合のパラメータ

# Fields
- `network::NetworkParams`: 各ネットワークの共通パラメータ
- `C::Array{Float64,3}`: 高次結合テンソル (N×N×N)
- `epsilon::Float64`: 結合強度 (0 < ε << 1)
"""
struct ThreeNetworkParams
    network::NetworkParams
    C::Array{Float64,3}
    epsilon::Float64

    function ThreeNetworkParams(network::NetworkParams, C::Array{Float64,3}, epsilon::Float64)
        N = network.N
        size(C) == (N, N, N) || throw(DimensionMismatch("C must be $N x $N x $N"))
        0 < epsilon < 0.1 || @warn "epsilon should satisfy 0 < ε << 1 for phase reduction validity"
        !any(isnan, C) || throw(DomainError(C, "C must not contain NaN"))
        !any(isinf, C) || throw(DomainError(C, "C must not contain Inf"))
        new(network, C, epsilon)
    end
end

# =============================================================================
# 計算結果構造体（Immutable）
# =============================================================================

"""
    StablePeriodicSolution

安定周期解（SPS）の計算結果

# Fields
- `N::Int`: 振動子数
- `Ntheta::Int`: 位相分割数
- `T::Float64`: 周期 [秒]
- `omega::Float64`: 角周波数 [rad/s] (ω = 2π/T)
- `Xs::Matrix{Float64}`: 安定周期解 (Ds*N × Ntheta)
- `params::NetworkParams`: 使用したパラメータ
"""
struct StablePeriodicSolution
    N::Int
    Ntheta::Int
    T::Float64
    omega::Float64
    Xs::Matrix{Float64}
    params::NetworkParams
    X0_initial::Vector{Float64}  # 緩和前の初期値

    function StablePeriodicSolution(N::Int, Ntheta::Int, T::Float64, omega::Float64,
        Xs::Matrix{Float64}, params::NetworkParams,
        X0_initial::Vector{Float64}=zeros(Ds * N))
        T > 0 || throw(DomainError(T, "Period T must be positive"))
        size(Xs) == (Ds * N, Ntheta) ||
            throw(DimensionMismatch("Xs must be $(Ds*N) x $Ntheta"))
        isapprox(omega, 2π / T; rtol=1e-10) ||
            throw(DomainError(omega, "omega must equal 2π/T"))
        length(X0_initial) == Ds * N ||
            throw(DimensionMismatch("X0_initial must have length $(Ds*N)"))
        new(N, Ntheta, T, omega, Xs, params, X0_initial)
    end
end

"""
    PhaseSensitivityFunction

位相感受関数（PSF）の計算結果

# Fields
- `N::Int`: 振動子数
- `Ntheta::Int`: 位相分割数
- `Q::Matrix{Float64}`: 位相感受関数 (Ds*N × Ntheta)
- `T::Float64`: 周期
- `omega::Float64`: 角周波数
"""
struct PhaseSensitivityFunction
    N::Int
    Ntheta::Int
    Q::Matrix{Float64}
    T::Float64
    omega::Float64

    function PhaseSensitivityFunction(N::Int, Ntheta::Int, Q::Matrix{Float64},
        T::Float64, omega::Float64)
        size(Q) == (Ds * N, Ntheta) ||
            throw(DimensionMismatch("Q must be $(Ds*N) x $Ntheta"))
        new(N, Ntheta, Q, T, omega)
    end
end

"""
    IndividualPCF

個別位相結合関数 γ_ijk(φ, ψ)

# Fields
- `gamma::Matrix{Float64}`: γ_ijk (Ntheta × Ntheta)
- `i::Int`: 基底振動子インデックス
- `j::Int`: 第1結合振動子インデックス
- `k::Int`: 第2結合振動子インデックス
"""
struct IndividualPCF
    gamma::Matrix{Float64}
    i::Int
    j::Int
    k::Int
end

"""
    TotalPCF

全体位相結合関数 Γ(φ, ψ)

# Fields
- `N::Int`: 振動子数
- `Ntheta::Int`: 位相分割数
- `Gamma::Matrix{Float64}`: Γ(φ, ψ) (Ntheta × Ntheta)
- `C::Array{Float64,3}`: 使用した結合テンソル
- `dGamma_dphi::Float64`: ∂Γ/∂φ at (0,0) = Γ₁
- `dGamma_dpsi::Float64`: ∂Γ/∂ψ at (0,0) = Γ₂
"""
struct TotalPCF
    N::Int
    Ntheta::Int
    Gamma::Matrix{Float64}
    C::Array{Float64,3}
    dGamma_dphi::Float64
    dGamma_dpsi::Float64

    function TotalPCF(N::Int, Ntheta::Int, Gamma::Matrix{Float64}, C::Array{Float64,3},
        dGamma_dphi::Float64, dGamma_dpsi::Float64)
        size(Gamma) == (Ntheta, Ntheta) ||
            throw(DimensionMismatch("Gamma must be $Ntheta x $Ntheta"))
        size(C) == (N, N, N) ||
            throw(DimensionMismatch("C must be $N x $N x $N"))
        new(N, Ntheta, Gamma, C, dGamma_dphi, dGamma_dpsi)
    end
end

"""
    OptimizationResult

最適化結果

# Fields
- `N::Int`: 振動子数
- `C_opt::Array{Float64,3}`: 最適化された結合テンソル
- `opt_type::Symbol`: 最適化タイプ (:linear または :rotation)
- `target::Float64`: 目標値 (q または μ)
- `achieved::Float64`: 達成値
- `frobenius_norm::Float64`: ||C||_F
"""
struct OptimizationResult
    N::Int
    C_opt::Array{Float64,3}
    opt_type::Symbol
    target::Float64
    achieved::Float64
    frobenius_norm::Float64

    function OptimizationResult(N::Int, C_opt::Array{Float64,3}, opt_type::Symbol,
        target::Float64, achieved::Float64, frobenius_norm::Float64)
        opt_type in (:linear, :rotation) ||
            throw(ArgumentError("opt_type must be :linear or :rotation"))
        size(C_opt) == (N, N, N) ||
            throw(DimensionMismatch("C_opt must be $N x $N x $N"))
        new(N, C_opt, opt_type, target, achieved, frobenius_norm)
    end
end

"""
    PhaseTimeSeries

位相時系列データ

# Fields
- `time::Vector{Float64}`: 時間軸
- `phases::Matrix{Float64}`: 位相 (3 × T_steps) for A, B, C
- `phase_diffs::Matrix{Float64}`: 位相差 (2 × T_steps) for φ₁, φ₂
"""
struct PhaseTimeSeries
    time::Vector{Float64}
    phases::Matrix{Float64}
    phase_diffs::Matrix{Float64}

    function PhaseTimeSeries(time::Vector{Float64}, phases::Matrix{Float64},
        phase_diffs::Matrix{Float64})
        T_steps = length(time)
        size(phases, 2) == T_steps ||
            throw(DimensionMismatch("phases must have $T_steps columns"))
        size(phases, 1) == 3 ||
            throw(DimensionMismatch("phases must have 3 rows (A, B, C)"))
        size(phase_diffs, 2) == T_steps ||
            throw(DimensionMismatch("phase_diffs must have $T_steps columns"))
        size(phase_diffs, 1) == 2 ||
            throw(DimensionMismatch("phase_diffs must have 2 rows (φ₁, φ₂)"))
        new(time, phases, phase_diffs)
    end
end

# =============================================================================
# ユーティリティ関数（Pure）
# =============================================================================

"""
    get_oscillator_state(X::Vector{Float64}, i::Int) -> OscillatorState

状態ベクトルからi番目の振動子状態を取得（Pure関数）

# Arguments
- `X::Vector{Float64}`: ネットワーク状態 (length = Ds*N)
- `i::Int`: 振動子インデックス (1-based)

# Returns
- `OscillatorState`: (u_i, v_i)
"""
function get_oscillator_state(X::Vector{Float64}, i::Int)
    return OscillatorState(X[2i-1], X[2i])
end

"""
    set_oscillator_state!(X::Vector{Float64}, i::Int, state::OscillatorState)

状態ベクトルのi番目の振動子状態を設定
"""
function set_oscillator_state!(X::Vector{Float64}, i::Int, state::OscillatorState)
    X[2i-1] = state[1]
    X[2i] = state[2]
end

"""
    index_to_phase(k::Int, Ntheta::Int) -> Float64

位相インデックスから位相値への変換（Pure関数）

θ = 2π(k - 1) / (Ntheta - 1)  where k ∈ 1:Ntheta
"""
function index_to_phase(k::Int, Ntheta::Int)
    return 2π * (k - 1) / (Ntheta - 1)
end

"""
    phase_to_index(theta::Float64, Ntheta::Int) -> Int

位相値から最近傍インデックスへの変換（Pure関数）
"""
function phase_to_index(theta::Float64, Ntheta::Int)
    theta_normalized = mod(theta, 2π)
    k = round(Int, (Ntheta - 1) * theta_normalized / (2π)) + 1
    return clamp(k, 1, Ntheta)
end

"""
    normalize_phase(theta::Float64) -> Float64

位相を [-π, π] の範囲に正規化（Pure関数）
"""
function normalize_phase(theta::Float64)
    theta = mod(theta, 2π)
    return theta > π ? theta - 2π : theta
end

"""
    normalize_phase_positive(theta::Float64) -> Float64

位相を [0, 2π) の範囲に正規化（Pure関数）
"""
function normalize_phase_positive(theta::Float64)
    return mod(theta, 2π)
end
