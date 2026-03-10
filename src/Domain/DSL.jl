"""
    Domain/DSL.jl - 埋め込みDSL定義 [Pure]

FDDのeDSL（embedded Domain Specific Language）
計算処理をコマンド（値）として表現し、実行をインタープリターに委譲
"""

# =============================================================================
# 計算コマンド基底型
# =============================================================================

"""
計算コマンドの基底型（eDSL）

すべての計算は純粋なデータ構造として表現される。
実際の実行はインタープリター層で行う。
"""
abstract type ComputationCommand end

# =============================================================================
# 安定周期解計算コマンド
# =============================================================================

"""
    ComputeSPS

安定周期解（Stable Periodic Solution）計算コマンド

# Fields
- `params::NetworkParams`: ネットワークパラメータ
- `Ntheta::Int`: 位相分割数
- `t_relax::Float64`: 緩和時間
"""
struct ComputeSPS <: ComputationCommand
    params::NetworkParams
    Ntheta::Int
    t_relax::Float64

    function ComputeSPS(params::NetworkParams, Ntheta::Int, t_relax::Float64)
        Ntheta > 1 || throw(DomainError(Ntheta, "Ntheta must be > 1"))
        t_relax > 0 || throw(DomainError(t_relax, "t_relax must be positive"))
        new(params, Ntheta, t_relax)
    end
end

"""デフォルト値でComputeSPSを作成"""
ComputeSPS(params::NetworkParams) = ComputeSPS(params, Nθ_DEFAULT, TREL)

# =============================================================================
# 位相感受関数計算コマンド
# =============================================================================

"""
    ComputePSF

位相感受関数（Phase Sensitivity Function）計算コマンド

FDD: 結合タイプを型で表現し、多重ディスパッチで計算方法を切り替える

# Fields
- `sps::StablePeriodicSolution`: 安定周期解
- `n_iterations::Int`: 随伴方程式の反復回数
- `coupling_type::CouplingType`: 結合タイプ（計算方法の選択に使用）

# Coupling Types
- `DiffusiveCoupling()`: 拡散結合用の簡略化計算（PSF.mスタイル）
- `GeneralCoupling()`: 一般的な非線形結合用の完全計算
"""
struct ComputePSF <: ComputationCommand
    sps::StablePeriodicSolution
    n_iterations::Int
    coupling_type::CouplingType

    function ComputePSF(sps::StablePeriodicSolution, n_iterations::Int,
                        coupling_type::CouplingType)
        n_iterations > 0 || throw(DomainError(n_iterations, "n_iterations must be positive"))
        new(sps, n_iterations, coupling_type)
    end
end

"""拡散結合（デフォルト）でComputePSFを作成"""
ComputePSF(sps::StablePeriodicSolution) = ComputePSF(sps, REP, DiffusiveCoupling())

"""反復回数指定でComputePSFを作成（拡散結合）"""
ComputePSF(sps::StablePeriodicSolution, n_iterations::Int) =
    ComputePSF(sps, n_iterations, DiffusiveCoupling())

"""結合タイプ指定でComputePSFを作成（デフォルト反復回数）"""
ComputePSF(sps::StablePeriodicSolution, coupling_type::CouplingType) =
    ComputePSF(sps, REP, coupling_type)

# =============================================================================
# 個別位相結合関数計算コマンド
# =============================================================================

"""
    ComputeIndividualPCF

個別位相結合関数 γ_ijk(φ, ψ) 計算コマンド

# Fields
- `sps::StablePeriodicSolution`: 安定周期解
- `psf::PhaseSensitivityFunction`: 位相感受関数
- `i::Int`: 基底振動子インデックス
- `j::Int`: 第1結合振動子インデックス
- `k::Int`: 第2結合振動子インデックス
"""
struct ComputeIndividualPCF <: ComputationCommand
    sps::StablePeriodicSolution
    psf::PhaseSensitivityFunction
    i::Int
    j::Int
    k::Int

    function ComputeIndividualPCF(sps::StablePeriodicSolution, psf::PhaseSensitivityFunction,
                                   i::Int, j::Int, k::Int)
        N = sps.N
        1 <= i <= N || throw(DomainError(i, "i must be in 1:$N"))
        1 <= j <= N || throw(DomainError(j, "j must be in 1:$N"))
        1 <= k <= N || throw(DomainError(k, "k must be in 1:$N"))
        sps.Ntheta == psf.Ntheta || throw(DimensionMismatch("Ntheta mismatch"))
        new(sps, psf, i, j, k)
    end
end

# =============================================================================
# 全体位相結合関数計算コマンド
# =============================================================================

"""
    ComputeTotalPCF

全体位相結合関数 Γ(φ, ψ) 計算コマンド

# Fields
- `individual_pcfs::Vector{IndividualPCF}`: 個別位相結合関数のリスト
- `C::Array{Float64,3}`: 結合テンソル
- `N::Int`: 振動子数
- `Ntheta::Int`: 位相分割数
"""
struct ComputeTotalPCF <: ComputationCommand
    individual_pcfs::Vector{IndividualPCF}
    C::Array{Float64,3}
    N::Int
    Ntheta::Int

    function ComputeTotalPCF(individual_pcfs::Vector{IndividualPCF}, C::Array{Float64,3},
                              N::Int, Ntheta::Int)
        size(C) == (N, N, N) || throw(DimensionMismatch("C must be $N x $N x $N"))
        !isempty(individual_pcfs) || throw(ArgumentError("individual_pcfs cannot be empty"))
        new(individual_pcfs, C, N, Ntheta)
    end
end

# =============================================================================
# 最適化コマンド
# =============================================================================

"""最適化コマンドの基底型"""
abstract type OptimizationCommand <: ComputationCommand end

"""
    OptimizeLinearStability

線形安定性最適化コマンド

目的: Λ = -3/2(Γ₁ + Γ₂) を最大化（最も負にする）

# Fields
- `individual_pcfs::Vector{IndividualPCF}`: 個別位相結合関数
- `target_q::Float64`: 目標値 q (Γ₁ + Γ₂ = q を達成)
- `N::Int`: 振動子数
- `Ntheta::Int`: 位相分割数
"""
struct OptimizeLinearStability <: OptimizationCommand
    individual_pcfs::Vector{IndividualPCF}
    target_q::Float64
    N::Int
    Ntheta::Int

    function OptimizeLinearStability(individual_pcfs::Vector{IndividualPCF},
                                      target_q::Float64, N::Int, Ntheta::Int)
        !isempty(individual_pcfs) || throw(ArgumentError("individual_pcfs cannot be empty"))
        new(individual_pcfs, target_q, N, Ntheta)
    end
end

"""
    OptimizeRotation

回転特性最適化コマンド

目的: R = √3/2 |Γ₁ - Γ₂| を最大化

# Fields
- `individual_pcfs::Vector{IndividualPCF}`: 個別位相結合関数
- `target_mu::Float64`: 目標値 μ (|Γ₁ - Γ₂| = μ を達成)
- `N::Int`: 振動子数
- `Ntheta::Int`: 位相分割数
"""
struct OptimizeRotation <: OptimizationCommand
    individual_pcfs::Vector{IndividualPCF}
    target_mu::Float64
    N::Int
    Ntheta::Int

    function OptimizeRotation(individual_pcfs::Vector{IndividualPCF},
                               target_mu::Float64, N::Int, Ntheta::Int)
        !isempty(individual_pcfs) || throw(ArgumentError("individual_pcfs cannot be empty"))
        target_mu >= 0 || throw(DomainError(target_mu, "target_mu must be non-negative"))
        new(individual_pcfs, target_mu, N, Ntheta)
    end
end

# =============================================================================
# シミュレーションコマンド
# =============================================================================

"""
    SimulatePhaseDynamics

縮約位相方程式のシミュレーションコマンド

# Fields
- `total_pcf::TotalPCF`: 全体位相結合関数
- `initial_phi1::Float64`: φ₁の初期値
- `initial_phi2::Float64`: φ₂の初期値
- `epsilon::Float64`: 結合強度
- `t_start::Float64`: シミュレーション開始時刻
- `t_end::Float64`: シミュレーション終了時刻
"""
struct SimulatePhaseDynamics <: ComputationCommand
    total_pcf::TotalPCF
    initial_phi1::Float64
    initial_phi2::Float64
    epsilon::Float64
    t_start::Float64
    t_end::Float64

    function SimulatePhaseDynamics(total_pcf::TotalPCF, initial_phi1::Float64,
                                    initial_phi2::Float64, epsilon::Float64,
                                    t_start::Float64, t_end::Float64)
        0 < epsilon < 0.1 || @warn "epsilon should satisfy 0 < ε << 1"
        t_end > t_start || throw(DomainError((t_start, t_end), "t_end must be > t_start"))
        new(total_pcf, initial_phi1, initial_phi2, epsilon, t_start, t_end)
    end
end

"""
    SimulateFullDynamics

元力学系の3ネットワーク結合シミュレーションコマンド

# Fields
- `three_network_params::ThreeNetworkParams`: 3ネットワークパラメータ
- `initial_states::Vector{Vector{Float64}}`: 初期状態 (3つのネットワーク)
- `t_start::Float64`: シミュレーション開始時刻
- `t_end::Float64`: シミュレーション終了時刻
"""
struct SimulateFullDynamics <: ComputationCommand
    three_network_params::ThreeNetworkParams
    initial_states::Vector{Vector{Float64}}
    t_start::Float64
    t_end::Float64

    function SimulateFullDynamics(params::ThreeNetworkParams,
                                   initial_states::Vector{Vector{Float64}},
                                   t_start::Float64, t_end::Float64)
        N = params.network.N
        length(initial_states) == N_NET || throw(DimensionMismatch("Must have $N_NET initial states"))
        for (i, state) in enumerate(initial_states)
            length(state) == Ds * N ||
                throw(DimensionMismatch("Initial state $i must have length $(Ds*N)"))
        end
        t_end > t_start || throw(DomainError((t_start, t_end), "t_end must be > t_start"))
        new(params, initial_states, t_start, t_end)
    end
end

# =============================================================================
# データ保存・読み込みコマンド
# =============================================================================

"""IO操作コマンドの基底型"""
abstract type IOCommand <: ComputationCommand end

"""
    SaveData

データ保存コマンド

# Fields
- `data`: 保存するデータ
- `filename::String`: ファイル名
- `format::Symbol`: フォーマット (:jld2, :csv)
"""
struct SaveData <: IOCommand
    data::Any
    filename::String
    format::Symbol

    function SaveData(data::Any, filename::String, format::Symbol=:jld2)
        format in (:jld2, :csv) || throw(ArgumentError("format must be :jld2 or :csv"))
        new(data, filename, format)
    end
end

"""
    LoadData

データ読み込みコマンド

# Fields
- `filename::String`: ファイル名
- `expected_type::Type`: 期待される型
"""
struct LoadData <: IOCommand
    filename::String
    expected_type::Type
end

# =============================================================================
# 計算スクリプト（コマンドのシーケンス）
# =============================================================================

"""計算スクリプト = コマンドのシーケンス"""
const ComputationScript = Vector{ComputationCommand}

# =============================================================================
# 計算結果型（Either風）
# =============================================================================

"""計算結果を表す型（Either風）"""
abstract type ComputationResult{T} end

"""成功結果"""
struct Success{T} <: ComputationResult{T}
    value::T
end

"""失敗結果"""
struct Failure{T} <: ComputationResult{T}
    error::String
    context::Dict{String, Any}

    Failure{T}(error::String) where T = new{T}(error, Dict{String, Any}())
    Failure{T}(error::String, context::Dict{String, Any}) where T = new{T}(error, context)
end

# ユーティリティ関数
"""結果が成功かどうかを判定"""
is_success(::Success) = true
is_success(::Failure) = false

"""成功結果から値を取り出す（失敗時はエラー）"""
function unwrap(r::Success{T}) where T
    return r.value
end

function unwrap(r::Failure{T}) where T
    throw(ErrorException("Unwrap failed: $(r.error)"))
end

"""成功結果から値を取り出す（失敗時はデフォルト値）"""
function unwrap_or(r::Success{T}, default::T) where T
    return r.value
end

function unwrap_or(r::Failure{T}, default::T) where T
    return default
end

"""結果に関数を適用（成功時のみ）"""
function map_result(f, r::Success{T}) where T
    return Success(f(r.value))
end

function map_result(f, r::Failure{T}) where T
    return r
end

"""bind演算（モナド風チェイン）"""
function bind_result(r::Success{T}, f) where T
    try
        return f(r.value)
    catch e
        return Failure{Any}(string(e), Dict("exception" => e))
    end
end

function bind_result(r::Failure{T}, f) where T
    return r
end

# |> 演算子でのチェインをサポート
Base.:>>(r::ComputationResult, f) = bind_result(r, f)

# =============================================================================
# スクリプトビルダー（Pure関数）
# =============================================================================

"""
    build_standard_pipeline(params; Ntheta, t_relax, n_iterations) -> ComputationScript

標準的な計算パイプラインを生成

SPS → PSF → IndividualPCFs の計算フローを定義
"""
function build_standard_pipeline(params::NetworkParams;
                                  Ntheta::Int=Nθ_DEFAULT,
                                  t_relax::Float64=TREL,
                                  n_iterations::Int=REP)
    return ComputationScript([
        ComputeSPS(params, Ntheta, t_relax)
        # 注: PSF, IndividualPCF計算は前のコマンドの結果に依存するため、
        # インタープリターで動的に生成される
    ])
end

"""
    build_optimization_pipeline(params, target_q; ...) -> ComputationScript

最適化パイプラインを生成

SPS → PSF → IndividualPCFs → Optimization の計算フローを定義
"""
function build_optimization_pipeline(params::NetworkParams,
                                      target_q::Float64;
                                      Ntheta::Int=Nθ_DEFAULT,
                                      optimization_type::Symbol=:linear)
    optimization_type in (:linear, :rotation) ||
        throw(ArgumentError("optimization_type must be :linear or :rotation"))

    return ComputationScript([
        ComputeSPS(params, Ntheta, TREL)
        # 最適化コマンドは動的に生成
    ])
end

"""
    build_comparison_pipeline(params; ...) -> ComputationScript

比較実験用パイプラインを生成

一様結合、線形安定性最適化、回転特性最適化を比較するためのスクリプト
"""
function build_comparison_pipeline(params::NetworkParams;
                                    Ntheta::Int=Nθ_DEFAULT,
                                    target_q::Float64=0.1,
                                    target_mu::Float64=0.5)
    return ComputationScript([
        ComputeSPS(params, Ntheta, TREL)
        # 3種類の結合条件での計算は動的に生成
    ])
end

# =============================================================================
# コマンド表示（デバッグ用）
# =============================================================================

Base.show(io::IO, cmd::ComputeSPS) =
    print(io, "ComputeSPS(N=$(cmd.params.N), Ntheta=$(cmd.Ntheta))")

Base.show(io::IO, cmd::ComputePSF) =
    print(io, "ComputePSF(iterations=$(cmd.n_iterations))")

Base.show(io::IO, cmd::ComputeIndividualPCF) =
    print(io, "ComputeIndividualPCF(i=$(cmd.i), j=$(cmd.j), k=$(cmd.k))")

Base.show(io::IO, cmd::ComputeTotalPCF) =
    print(io, "ComputeTotalPCF(N=$(cmd.N), Ntheta=$(cmd.Ntheta))")

Base.show(io::IO, cmd::OptimizeLinearStability) =
    print(io, "OptimizeLinearStability(target_q=$(cmd.target_q))")

Base.show(io::IO, cmd::OptimizeRotation) =
    print(io, "OptimizeRotation(target_mu=$(cmd.target_mu))")

Base.show(io::IO, cmd::SimulatePhaseDynamics) =
    print(io, "SimulatePhaseDynamics(ε=$(cmd.epsilon), t=$(cmd.t_start):$(cmd.t_end))")

Base.show(io::IO, cmd::SimulateFullDynamics) =
    print(io, "SimulateFullDynamics(t=$(cmd.t_start):$(cmd.t_end))")

Base.show(io::IO, r::Success{T}) where T = print(io, "Success{$T}(...)")
Base.show(io::IO, r::Failure{T}) where T = print(io, "Failure{$T}($(r.error))")

