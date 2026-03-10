"""
    Interpreters/PhaseInterpreter.jl - 位相計算インタープリター [Impure]

位相結合関数計算と最適化コマンドを実行するインタープリター
"""

using LinearAlgebra
using Printf

# =============================================================================
# 個別位相結合関数（IndividualPCF）計算
# =============================================================================

"""
    interpret(cmd::ComputeIndividualPCF) -> ComputationResult{IndividualPCF}

個別位相結合関数 γ_ijk(φ, ψ) 計算コマンドを実行
"""
function interpret(cmd::ComputeIndividualPCF)
    try
        sps = cmd.sps
        psf = cmd.psf
        i, j, k = cmd.i, cmd.j, cmd.k

        # Pure関数を呼び出し
        pcf = compute_individual_pcf(sps, psf, i, j, k)

        return Success(pcf)

    catch e
        return Failure{IndividualPCF}("IndividualPCF computation failed for ($i,$j,$k): $e",
            Dict("exception" => e, "indices" => (i, j, k)))
    end
end

"""
    compute_all_individual_pcfs_parallel(sps, psf; show_progress=true)
        -> Vector{IndividualPCF}

全ての個別PCFを計算（進捗表示付き）
"""
function compute_all_individual_pcfs_parallel(sps::StablePeriodicSolution,
                                               psf::PhaseSensitivityFunction;
                                               show_progress::Bool=true)
    N = sps.N
    total = N^3
    pcfs = Vector{IndividualPCF}(undef, total)

    idx = 1
    for i in 1:N
        for j in 1:N
            for k in 1:N
                if show_progress && idx % 100 == 0
                    @printf("Computing PCF %d/%d (%.1f%%)\r", idx, total, 100*idx/total)
                end

                pcfs[idx] = compute_individual_pcf(sps, psf, i, j, k)
                idx += 1
            end
        end
    end

    if show_progress
        println("\nCompleted computing $total individual PCFs")
    end

    return pcfs
end

# =============================================================================
# 全体位相結合関数（TotalPCF）計算
# =============================================================================

"""
    interpret(cmd::ComputeTotalPCF) -> ComputationResult{TotalPCF}

全体位相結合関数 Γ(φ, ψ) 計算コマンドを実行
"""
function interpret(cmd::ComputeTotalPCF)
    try
        individual_pcfs = cmd.individual_pcfs
        C = cmd.C
        N = cmd.N
        Ntheta = cmd.Ntheta

        # Pure関数を呼び出し
        total_pcf = compute_total_pcf(individual_pcfs, C, N, Ntheta)

        # 検証
        validation_result = validate_pcf(total_pcf)
        if !is_success(validation_result)
            @warn "TotalPCF validation warning: $(validation_result.error)"
        end

        return Success(total_pcf)

    catch e
        return Failure{TotalPCF}("TotalPCF computation failed: $e",
            Dict("exception" => e))
    end
end

# =============================================================================
# 最適化コマンド実行
# =============================================================================

"""
    interpret(cmd::OptimizeLinearStability) -> ComputationResult{OptimizationResult}

線形安定性最適化コマンドを実行
"""
function interpret(cmd::OptimizeLinearStability)
    try
        individual_pcfs = cmd.individual_pcfs
        target_q = cmd.target_q
        N = cmd.N
        Ntheta = cmd.Ntheta

        # Pure関数を呼び出し
        result = optimize_linear_stability(individual_pcfs, target_q, N, Ntheta)

        # 検証
        validation_result = validate_optimization_result(result)
        if !is_success(validation_result)
            @warn "Optimization validation warning: $(validation_result.error)"
        end

        # 性能評価をログ出力
        performance = evaluate_coupling_performance(result.C_opt, individual_pcfs, N, Ntheta)
        @info "Linear stability optimization completed" target_q=target_q achieved=result.achieved frobenius_norm=result.frobenius_norm Lambda=performance.Lambda

        return Success(result)

    catch e
        return Failure{OptimizationResult}("Linear stability optimization failed: $e",
            Dict("exception" => e))
    end
end

"""
    interpret(cmd::OptimizeRotation) -> ComputationResult{OptimizationResult}

回転特性最適化コマンドを実行
"""
function interpret(cmd::OptimizeRotation)
    try
        individual_pcfs = cmd.individual_pcfs
        target_mu = cmd.target_mu
        N = cmd.N
        Ntheta = cmd.Ntheta

        # Pure関数を呼び出し
        result = optimize_rotation(individual_pcfs, target_mu, N, Ntheta)

        # 検証
        validation_result = validate_optimization_result(result)
        if !is_success(validation_result)
            @warn "Optimization validation warning: $(validation_result.error)"
        end

        # 性能評価をログ出力
        performance = evaluate_coupling_performance(result.C_opt, individual_pcfs, N, Ntheta)
        @info "Rotation optimization completed" target_mu=target_mu achieved=result.achieved frobenius_norm=result.frobenius_norm R=performance.R

        return Success(result)

    catch e
        return Failure{OptimizationResult}("Rotation optimization failed: $e",
            Dict("exception" => e))
    end
end

# =============================================================================
# パイプライン実行
# =============================================================================

"""
    run_full_pipeline(params; Ntheta, t_relax, n_iterations, coupling_type, target_value)
        -> NamedTuple

SPS → PSF → IndividualPCFs → Optimization の完全パイプラインを実行
"""
function run_full_pipeline(params::NetworkParams;
                            Ntheta::Int=Nθ_DEFAULT,
                            t_relax::Float64=TREL,
                            n_iterations::Int=REP,
                            coupling_type::Symbol=:uniform,
                            target_value::Float64=0.0)
    results = Dict{Symbol, Any}()

    # Step 1: SPS計算
    @info "Step 1: Computing Stable Periodic Solution..."
    sps_cmd = ComputeSPS(params, Ntheta, t_relax)
    sps_result = interpret(sps_cmd)

    if !is_success(sps_result)
        error("SPS computation failed: $(sps_result.error)")
    end
    sps = unwrap(sps_result)
    results[:sps] = sps
    @info "SPS computed" T=sps.T omega=sps.omega

    # Step 2: PSF計算
    @info "Step 2: Computing Phase Sensitivity Function..."
    psf_cmd = ComputePSF(sps, n_iterations)
    psf_result = interpret(psf_cmd)

    if !is_success(psf_result)
        error("PSF computation failed: $(psf_result.error)")
    end
    psf = unwrap(psf_result)
    results[:psf] = psf
    @info "PSF computed"

    # Step 3: 全IndividualPCF計算
    @info "Step 3: Computing Individual Phase Coupling Functions..."
    individual_pcfs = compute_all_individual_pcfs_parallel(sps, psf)
    results[:individual_pcfs] = individual_pcfs

    # Step 4: 結合テンソルの選択と最適化
    N = params.N

    if coupling_type == :uniform
        @info "Using uniform coupling tensor"
        C = create_diagonal_coupling_tensor(N)
        results[:C] = C

    elseif coupling_type == :linear
        @info "Optimizing for linear stability with target_q = $target_value"
        opt_cmd = OptimizeLinearStability(individual_pcfs, target_value, N, Ntheta)
        opt_result = interpret(opt_cmd)

        if !is_success(opt_result)
            error("Linear optimization failed: $(opt_result.error)")
        end
        opt = unwrap(opt_result)
        results[:optimization] = opt
        C = opt.C_opt

    elseif coupling_type == :rotation
        @info "Optimizing for rotation with target_mu = $target_value"
        opt_cmd = OptimizeRotation(individual_pcfs, target_value, N, Ntheta)
        opt_result = interpret(opt_cmd)

        if !is_success(opt_result)
            error("Rotation optimization failed: $(opt_result.error)")
        end
        opt = unwrap(opt_result)
        results[:optimization] = opt
        C = opt.C_opt

    else
        error("Unknown coupling_type: $coupling_type")
    end

    # Step 5: TotalPCF計算
    @info "Step 5: Computing Total Phase Coupling Function..."
    pcf_cmd = ComputeTotalPCF(individual_pcfs, C, N, Ntheta)
    pcf_result = interpret(pcf_cmd)

    if !is_success(pcf_result)
        error("TotalPCF computation failed: $(pcf_result.error)")
    end
    total_pcf = unwrap(pcf_result)
    results[:total_pcf] = total_pcf

    # 安定性指標を計算
    metrics = compute_stability_metrics(total_pcf)
    results[:metrics] = metrics

    @info "Pipeline completed" Lambda=metrics.Lambda R=metrics.R is_stable=metrics.is_stable

    return (
        sps = sps,
        psf = psf,
        individual_pcfs = individual_pcfs,
        total_pcf = total_pcf,
        metrics = metrics,
        coupling_type = coupling_type,
        optimization = get(results, :optimization, nothing)
    )
end

"""
    run_comparison_experiment(params; Ntheta, target_q, target_mu)
        -> Dict{Symbol, NamedTuple}

一様結合、線形最適化、回転最適化の3条件を比較
"""
function run_comparison_experiment(params::NetworkParams;
                                    Ntheta::Int=Nθ_DEFAULT,
                                    target_q::Float64=-0.1,
                                    target_mu::Float64=0.5)
    N = params.N
    results = Dict{Symbol, NamedTuple}()

    # 共通のSPS, PSF, IndividualPCFsを計算
    @info "Computing common components (SPS, PSF, IndividualPCFs)..."

    sps_result = interpret(ComputeSPS(params, Ntheta, TREL))
    sps = unwrap(sps_result)

    psf_result = interpret(ComputePSF(sps, REP))
    psf = unwrap(psf_result)

    individual_pcfs = compute_all_individual_pcfs_parallel(sps, psf; show_progress=false)

    # 1. 一様結合
    @info "Evaluating uniform coupling..."
    C_uniform = create_diagonal_coupling_tensor(N)
    pcf_uniform = compute_total_pcf(individual_pcfs, C_uniform, N, Ntheta)
    metrics_uniform = compute_stability_metrics(pcf_uniform)
    results[:uniform] = (
        C = C_uniform,
        total_pcf = pcf_uniform,
        metrics = metrics_uniform
    )

    # 2. 線形安定性最適化
    @info "Optimizing for linear stability..."
    opt_linear = optimize_linear_stability(individual_pcfs, target_q, N, Ntheta)
    pcf_linear = compute_total_pcf(individual_pcfs, opt_linear.C_opt, N, Ntheta)
    metrics_linear = compute_stability_metrics(pcf_linear)
    results[:linear_opt] = (
        optimization = opt_linear,
        total_pcf = pcf_linear,
        metrics = metrics_linear
    )

    # 3. 回転特性最適化
    @info "Optimizing for rotation..."
    opt_rotation = optimize_rotation(individual_pcfs, target_mu, N, Ntheta)
    pcf_rotation = compute_total_pcf(individual_pcfs, opt_rotation.C_opt, N, Ntheta)
    metrics_rotation = compute_stability_metrics(pcf_rotation)
    results[:rotation_opt] = (
        optimization = opt_rotation,
        total_pcf = pcf_rotation,
        metrics = metrics_rotation
    )

    # 結果サマリー
    @info "=== Comparison Results ===" uniform_Lambda=metrics_uniform.Lambda uniform_R=metrics_uniform.R linear_Lambda=metrics_linear.Lambda linear_R=metrics_linear.R rotation_Lambda=metrics_rotation.Lambda rotation_R=metrics_rotation.R

    return results
end

# =============================================================================
# ユーティリティ
# =============================================================================

"""
    interpret_script(script::ComputationScript) -> Vector{ComputationResult}

計算スクリプト全体を実行
"""
function interpret_script(script::ComputationScript)
    results = ComputationResult[]

    for (i, cmd) in enumerate(script)
        @info "Executing command $i: $(typeof(cmd))"
        result = interpret(cmd)
        push!(results, result)

        if !is_success(result)
            @warn "Command $i failed: $(result.error)"
            break
        end
    end

    return results
end

"""
    get_predefined_coupling(N::Int, type::Symbol) -> Array{Float64,3}

事前定義された結合テンソルを取得
"""
function get_predefined_coupling(N::Int, type::Symbol)
    if N == 4
        if type == :uniform
            return copy(HO_INTER_C_UNI_4)
        elseif type == :linear
            return copy(HO_LOPT_C_4)
        elseif type == :rotation
            return copy(HO_ROPT_C_4)
        end
    end

    # 事前定義がない場合は生成
    if type == :uniform
        return create_diagonal_coupling_tensor(N)
    else
        error("No predefined $type coupling for N=$N. Run optimization instead.")
    end
end

