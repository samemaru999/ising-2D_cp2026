"""
    Domain/Validation.jl - バリデーションロジック [Pure]

FDDのドメインバリデーション
型レベルでの不変条件検証と計算結果の検証
"""

using LinearAlgebra

# =============================================================================
# バリデーションエラー型
# =============================================================================

"""
    ValidationError

バリデーションエラーを表す型

# Fields
- `field::Symbol`: 問題のあるフィールド名
- `message::String`: エラーメッセージ
- `value::Any`: 問題の値
"""
struct ValidationError <: Exception
    field::Symbol
    message::String
    value::Any
end

Base.showerror(io::IO, e::ValidationError) =
    print(io, "ValidationError: $(e.field) - $(e.message) (value: $(e.value))")

"""バリデーション結果型"""
const ValidationResult{T} = Union{T,Vector{ValidationError}}

# =============================================================================
# パラメータバリデーション
# =============================================================================

"""
    validate_fhn_params(delta, a, b) -> Union{Nothing, Vector{ValidationError}}

FHNパラメータの検証

# Invariants
- δ > 0
- 0 < b < 1
"""
function validate_fhn_params(delta::Float64, a::Float64, b::Float64)
    errors = ValidationError[]

    delta > 0 || push!(errors,
        ValidationError(:delta, "must be positive", delta))

    0 < b < 1 || push!(errors,
        ValidationError(:b, "must be in (0, 1)", b))

    # aに対する制約は特になし（任意の実数）

    isempty(errors) ? nothing : errors
end

"""
    validate_network_params(N, I, K) -> Union{Nothing, Vector{ValidationError}}

ネットワークパラメータの検証

# Invariants
- N > 0
- length(I) == N
- size(K) == (N, N)
- diag(K) == zeros(N)
"""
function validate_network_params(N::Int, I::Vector{Float64}, K::Matrix{Float64})
    errors = ValidationError[]

    N > 0 || push!(errors,
        ValidationError(:N, "must be positive", N))

    length(I) == N || push!(errors,
        ValidationError(:I, "must have length N=$N", length(I)))

    size(K) == (N, N) || push!(errors,
        ValidationError(:K, "must be $N x $N", size(K)))

    if size(K) == (N, N)
        all(diag(K) .== 0) || push!(errors,
            ValidationError(:K, "diagonal must be zero (no self-coupling)", diag(K)))
    end

    isempty(errors) ? nothing : errors
end

"""
    validate_coupling_tensor(C, N) -> Union{Nothing, Vector{ValidationError}}

結合テンソルの検証

# Invariants
- size(C) == (N, N, N)
- !any(isnan, C)
- !any(isinf, C)
"""
function validate_coupling_tensor(C::Array{Float64,3}, N::Int)
    errors = ValidationError[]

    size(C) == (N, N, N) || push!(errors,
        ValidationError(:C, "must be $N x $N x $N", size(C)))

    any(isnan, C) && push!(errors,
        ValidationError(:C, "must not contain NaN", "NaN detected"))

    any(isinf, C) && push!(errors,
        ValidationError(:C, "must not contain Inf", "Inf detected"))

    isempty(errors) ? nothing : errors
end

"""
    validate_epsilon(epsilon) -> Union{Nothing, Vector{ValidationError}}

結合強度εの検証

# Invariants
- 0 < ε < 0.1 (位相縮約の妥当性条件)
"""
function validate_epsilon(epsilon::Float64)
    errors = ValidationError[]

    epsilon > 0 || push!(errors,
        ValidationError(:epsilon, "must be positive", epsilon))

    epsilon < 0.1 || push!(errors,
        ValidationError(:epsilon, "should be << 1 for phase reduction validity", epsilon))

    isempty(errors) ? nothing : errors
end

# =============================================================================
# 計算結果バリデーション
# =============================================================================

"""
    validate_sps(sps; tol) -> ComputationResult{StablePeriodicSolution}

安定周期解の検証

# Checks
- 周期性: X(0) ≈ X(T)
- 角周波数の整合性: ω = 2π/T
- NaN/Inf の不在
"""
function validate_sps(sps::StablePeriodicSolution; tol::Float64=1e-6)
    errors = ValidationError[]

    # 周期性の確認
    endpoint_diff = norm(sps.Xs[:, 1] - sps.Xs[:, end])
    endpoint_diff < tol || push!(errors,
        ValidationError(:Xs, "must be periodic (endpoint diff = $endpoint_diff)", endpoint_diff))

    # 角周波数の整合性
    expected_omega = 2π / sps.T
    isapprox(sps.omega, expected_omega; rtol=1e-10) || push!(errors,
        ValidationError(:omega, "must equal 2π/T", (omega=sps.omega, expected=expected_omega)))

    # NaN/Inf チェック
    any(isnan, sps.Xs) && push!(errors,
        ValidationError(:Xs, "must not contain NaN", "NaN detected"))

    any(isinf, sps.Xs) && push!(errors,
        ValidationError(:Xs, "must not contain Inf", "Inf detected"))

    # 周期が正であること
    sps.T > 0 || push!(errors,
        ValidationError(:T, "must be positive", sps.T))

    if isempty(errors)
        Success(sps)
    else
        Failure{StablePeriodicSolution}("SPS validation failed",
            Dict{String,Any}("errors" => errors))
    end
end

"""
    validate_psf(psf, sps, F; tol) -> ComputationResult{PhaseSensitivityFunction}

位相感受関数の検証

# Checks
- 正規化条件: Q(θ) · F(X(θ)) = ω
- 周期性: Q(0) ≈ Q(T)
- NaN/Inf の不在

# Arguments
- `psf::PhaseSensitivityFunction`: 検証対象
- `sps::StablePeriodicSolution`: 対応する周期解
- `F::Function`: 力学系のベクトル場 (X -> dX/dt)
- `tol::Float64`: 許容誤差
"""
function validate_psf(psf::PhaseSensitivityFunction, sps::StablePeriodicSolution,
    F::Function; tol::Float64=1e-4)
    errors = ValidationError[]

    # 正規化条件の確認
    max_normalization_error = 0.0
    for k in 1:psf.Ntheta
        X = sps.Xs[:, k]
        Q = psf.Q[:, k]
        normalization = dot(Q, F(X))
        error = abs(normalization - psf.omega)
        max_normalization_error = max(max_normalization_error, error)
    end

    max_normalization_error < tol || push!(errors,
        ValidationError(:Q, "normalization Q·F = ω not satisfied (max error = $max_normalization_error)",
            max_normalization_error))

    # 周期性の確認
    endpoint_diff = norm(psf.Q[:, 1] - psf.Q[:, end])
    endpoint_diff < tol || push!(errors,
        ValidationError(:Q, "must be periodic (endpoint diff = $endpoint_diff)", endpoint_diff))

    # NaN/Inf チェック
    any(isnan, psf.Q) && push!(errors,
        ValidationError(:Q, "must not contain NaN", "NaN detected"))

    any(isinf, psf.Q) && push!(errors,
        ValidationError(:Q, "must not contain Inf", "Inf detected"))

    if isempty(errors)
        Success(psf)
    else
        Failure{PhaseSensitivityFunction}("PSF validation failed",
            Dict{String,Any}("errors" => errors, "max_normalization_error" => max_normalization_error))
    end
end

"""
    validate_psf_simple(psf; tol) -> ComputationResult{PhaseSensitivityFunction}

位相感受関数の簡易検証（ベクトル場なし）

# Checks
- 周期性: Q(0) ≈ Q(T)
- NaN/Inf の不在
"""
function validate_psf_simple(psf::PhaseSensitivityFunction; tol::Float64=1e-6)
    errors = ValidationError[]

    # 周期性の確認
    endpoint_diff = norm(psf.Q[:, 1] - psf.Q[:, end])
    endpoint_diff < tol || push!(errors,
        ValidationError(:Q, "must be periodic (endpoint diff = $endpoint_diff)", endpoint_diff))

    # NaN/Inf チェック
    any(isnan, psf.Q) && push!(errors,
        ValidationError(:Q, "must not contain NaN", "NaN detected"))

    any(isinf, psf.Q) && push!(errors,
        ValidationError(:Q, "must not contain Inf", "Inf detected"))

    if isempty(errors)
        Success(psf)
    else
        Failure{PhaseSensitivityFunction}("PSF validation failed", Dict{String,Any}("errors" => errors))
    end
end

"""
    validate_pcf(pcf; tol) -> ComputationResult{TotalPCF}

位相結合関数の検証

# Checks
- Γ(0, 0) = 0 (同期点での値)
- NaN/Inf の不在
- 対称性（オプション）
"""
function validate_pcf(pcf::TotalPCF; tol::Float64=1e-6)
    errors = ValidationError[]

    # 同期点での値
    gamma_at_origin = pcf.Gamma[1, 1]
    abs(gamma_at_origin) < tol || push!(errors,
        ValidationError(:Gamma, "Γ(0,0) should be approximately 0", gamma_at_origin))

    # NaN/Inf チェック
    any(isnan, pcf.Gamma) && push!(errors,
        ValidationError(:Gamma, "must not contain NaN", "NaN detected"))

    any(isinf, pcf.Gamma) && push!(errors,
        ValidationError(:Gamma, "must not contain Inf", "Inf detected"))

    # 導関数のチェック
    isnan(pcf.dGamma_dphi) && push!(errors,
        ValidationError(:dGamma_dphi, "must not be NaN", pcf.dGamma_dphi))

    isnan(pcf.dGamma_dpsi) && push!(errors,
        ValidationError(:dGamma_dpsi, "must not be NaN", pcf.dGamma_dpsi))

    if isempty(errors)
        Success(pcf)
    else
        Failure{TotalPCF}("PCF validation failed", Dict{String,Any}("errors" => errors))
    end
end

"""
    validate_optimization_result(result; tol) -> ComputationResult{OptimizationResult}

最適化結果の検証

# Checks
- 達成値が目標値に近いか
- 結合テンソルがNaN/Infを含まないか
"""
function validate_optimization_result(result::OptimizationResult; tol::Float64=1e-3)
    errors = ValidationError[]

    # 目標達成の確認
    isapprox(result.achieved, result.target; atol=tol) || push!(errors,
        ValidationError(:achieved,
            "target not achieved (target=$(result.target), achieved=$(result.achieved))",
            (target=result.target, achieved=result.achieved)))

    # 最適化テンソルのチェック
    any(isnan, result.C_opt) && push!(errors,
        ValidationError(:C_opt, "must not contain NaN", "NaN detected"))

    any(isinf, result.C_opt) && push!(errors,
        ValidationError(:C_opt, "must not contain Inf", "Inf detected"))

    # フロベニウスノルムが正であること
    result.frobenius_norm >= 0 || push!(errors,
        ValidationError(:frobenius_norm, "must be non-negative", result.frobenius_norm))

    if isempty(errors)
        Success(result)
    else
        Failure{OptimizationResult}("Optimization validation failed", Dict{String,Any}("errors" => errors))
    end
end

"""
    validate_phase_time_series(pts; tol) -> ComputationResult{PhaseTimeSeries}

位相時系列の検証

# Checks
- 時間が単調増加
- 位相が有限値
"""
function validate_phase_time_series(pts::PhaseTimeSeries; tol::Float64=1e-10)
    errors = ValidationError[]

    # 時間の単調増加
    all(diff(pts.time) .> 0) || push!(errors,
        ValidationError(:time, "must be strictly increasing", "non-monotonic"))

    # NaN/Inf チェック
    any(isnan, pts.phases) && push!(errors,
        ValidationError(:phases, "must not contain NaN", "NaN detected"))

    any(isinf, pts.phases) && push!(errors,
        ValidationError(:phases, "must not contain Inf", "Inf detected"))

    any(isnan, pts.phase_diffs) && push!(errors,
        ValidationError(:phase_diffs, "must not contain NaN", "NaN detected"))

    any(isinf, pts.phase_diffs) && push!(errors,
        ValidationError(:phase_diffs, "must not contain Inf", "Inf detected"))

    if isempty(errors)
        Success(pts)
    else
        Failure{PhaseTimeSeries}("PhaseTimeSeries validation failed", Dict{String,Any}("errors" => errors))
    end
end

# =============================================================================
# ユーティリティ関数
# =============================================================================

"""
    collect_errors(validations...) -> Union{Nothing, Vector{ValidationError}}

複数のバリデーション結果を集約
"""
function collect_errors(validations...)
    all_errors = ValidationError[]
    for v in validations
        if v !== nothing
            append!(all_errors, v)
        end
    end
    isempty(all_errors) ? nothing : all_errors
end

"""
    is_valid(validation_result) -> Bool

バリデーション結果が有効かどうかを判定
"""
is_valid(::Nothing) = true
is_valid(::Vector{ValidationError}) = false

"""
    assert_valid(validation_result, context::String)

バリデーション結果を検証し、エラーがあれば例外をスロー
"""
function assert_valid(validation_result, context::String="Validation")
    if validation_result !== nothing
        error_msgs = ["$context failed:"]
        for err in validation_result
            push!(error_msgs, "  - $(err.field): $(err.message)")
        end
        throw(DomainError(validation_result, join(error_msgs, "\n")))
    end
end

"""
    check_dimensions(expected, actual, name::String) -> Union{Nothing, ValidationError}

次元のチェック
"""
function check_dimensions(expected, actual, name::String)
    if expected != actual
        return ValidationError(Symbol(name), "dimension mismatch: expected $expected, got $actual", actual)
    end
    return nothing
end

"""
    check_positive(value, name::Symbol) -> Union{Nothing, ValidationError}

正値のチェック
"""
function check_positive(value::Number, name::Symbol)
    if value <= 0
        return ValidationError(name, "must be positive", value)
    end
    return nothing
end

"""
    check_in_range(value, lo, hi, name::Symbol) -> Union{Nothing, ValidationError}

範囲チェック
"""
function check_in_range(value::Number, lo::Number, hi::Number, name::Symbol)
    if !(lo < value < hi)
        return ValidationError(name, "must be in ($lo, $hi)", value)
    end
    return nothing
end

"""
    check_no_nan(array, name::Symbol) -> Union{Nothing, ValidationError}

NaNチェック
"""
function check_no_nan(array::AbstractArray, name::Symbol)
    if any(isnan, array)
        return ValidationError(name, "must not contain NaN", "NaN detected")
    end
    return nothing
end

"""
    check_no_inf(array, name::Symbol) -> Union{Nothing, ValidationError}

Infチェック
"""
function check_no_inf(array::AbstractArray, name::Symbol)
    if any(isinf, array)
        return ValidationError(name, "must not contain Inf", "Inf detected")
    end
    return nothing
end

"""
    check_finite(array, name::Symbol) -> Union{Nothing, Vector{ValidationError}}

有限値チェック（NaN と Inf 両方）
"""
function check_finite(array::AbstractArray, name::Symbol)
    errors = ValidationError[]
    nan_err = check_no_nan(array, name)
    inf_err = check_no_inf(array, name)
    nan_err !== nothing && push!(errors, nan_err)
    inf_err !== nothing && push!(errors, inf_err)
    isempty(errors) ? nothing : errors
end

