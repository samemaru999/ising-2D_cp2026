#!/usr/bin/env julia
"""
    main.jl - HONetSync 実行スクリプト

高次結合を持つ集団振動ネットワーク同期最適化シミュレーション

使用方法:
    julia --project=. scripts/main.jl [options]

オプション:
    --N=4           振動子数 (4 or 10)
    --mode=full     実行モード (full, linear, rotation, compare)
    --save          結果をファイルに保存
    --output=output 出力ディレクトリ

例:
    julia --project=. scripts/main.jl --N=4 --mode=full --save
    julia --project=. scripts/main.jl --mode=compare
"""

using HONetSync
using Printf
using Dates

# =============================================================================
# コマンドライン引数のパース
# =============================================================================

function parse_args(args)
    config = Dict{Symbol, Any}(
        :N => 4,
        :mode => :full,
        :save => false,
        :output => "output",
        :Ntheta => 101,
        :t_relax => 1000.0,
        :target_q => -0.1,
        :target_mu => 0.5,
        :epsilon => 2e-4,
        :t_end => 10000.0,
        :suffix => ""  # ファイル名に追加するサフィックス
    )

    for arg in args
        if startswith(arg, "--N=")
            config[:N] = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--mode=")
            config[:mode] = Symbol(split(arg, "=")[2])
        elseif arg == "--save"
            config[:save] = true
        elseif startswith(arg, "--output=")
            config[:output] = split(arg, "=")[2]
        elseif startswith(arg, "--Ntheta=")
            config[:Ntheta] = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--target-q=")
            config[:target_q] = parse(Float64, split(arg, "=")[2])
        elseif startswith(arg, "--target-mu=")
            config[:target_mu] = parse(Float64, split(arg, "=")[2])
        elseif startswith(arg, "--epsilon=")
            config[:epsilon] = parse(Float64, split(arg, "=")[2])
        elseif startswith(arg, "--suffix=")
            suffix_value = String(split(arg, "=")[2])
            if suffix_value == "timestamp" || suffix_value == "time"
                # 日本標準時(UTC+9)のタイムスタンプを生成
                jst = now(UTC) + Hour(9)
                config[:suffix] = Dates.format(jst, "yyyymmdd_HHMMSS")
            else
                config[:suffix] = suffix_value
            end
        elseif arg == "--help" || arg == "-h"
            println(@doc Main)
            exit(0)
        end
    end

    return config
end

# =============================================================================
# ステップ1: 安定周期解 (SPS) の計算
# =============================================================================

function compute_sps_step(params::NetworkParams, Ntheta::Int, t_relax::Float64)
    println("\n" * "="^60)
    println("Step 1: Computing Stable Periodic Solution (SPS)")
    println("="^60)
    println("  N = $(params.N), Ntheta = $Ntheta, t_relax = $t_relax")

    # eDSLコマンド作成
    cmd = ComputeSPS(params, Ntheta, t_relax)
    println("  Command: $cmd")

    # インタープリター実行
    @time result = interpret(cmd)

    if is_success(result)
        sps = unwrap(result)
        println("  Period T = $(sps.T)")
        println("  Angular frequency ω = $(sps.omega)")
        return sps
    else
        error("SPS computation failed: $(result.error)")
    end
end

# =============================================================================
# ステップ2: 位相感受関数 (PSF) の計算
# =============================================================================

function compute_psf_step(sps::StablePeriodicSolution)
    println("\n" * "="^60)
    println("Step 2: Computing Phase Sensitivity Function (PSF)")
    println("="^60)

    # eDSLコマンド作成
    cmd = ComputePSF(sps)
    println("  Command: $cmd")

    # インタープリター実行
    @time result = interpret(cmd)

    if is_success(result)
        psf = unwrap(result)
        println("  PSF computed successfully")
        println("  Q matrix size: $(size(psf.Q))")
        return psf
    else
        error("PSF computation failed: $(result.error)")
    end
end

# =============================================================================
# ステップ3: 個別位相結合関数 (Individual PCF) の計算
# =============================================================================

function compute_individual_pcfs_step(sps::StablePeriodicSolution, psf::PhaseSensitivityFunction)
    println("\n" * "="^60)
    println("Step 3: Computing Individual Phase Coupling Functions")
    println("="^60)

    N = sps.N
    total_combinations = N^3
    println("  Computing $total_combinations combinations (N=$N)")

    individual_pcfs = IndividualPCF[]

    @time for i in 1:N, j in 1:N, k in 1:N
        cmd = ComputeIndividualPCF(sps, psf, i, j, k)
        result = interpret(cmd)

        if is_success(result)
            push!(individual_pcfs, unwrap(result))
        else
            @warn "Failed to compute PCF for ($i, $j, $k): $(result.error)"
        end
    end

    println("  Computed $(length(individual_pcfs)) individual PCFs")
    return individual_pcfs
end

# =============================================================================
# ステップ4: 最適化
# =============================================================================

function optimize_linear_step(individual_pcfs::Vector{IndividualPCF}, target_q::Float64, N::Int, Ntheta::Int)
    println("\n" * "="^60)
    println("Step 4a: Linear Stability Optimization")
    println("="^60)
    println("  Target q = $target_q")

    # eDSLコマンド作成
    cmd = OptimizeLinearStability(individual_pcfs, target_q, N, Ntheta)
    println("  Command: $cmd")

    # インタープリター実行
    @time result = interpret(cmd)

    if is_success(result)
        opt = unwrap(result)
        println("  Achieved Λ = $(opt.achieved)")
        println("  Frobenius norm ||C|| = $(opt.frobenius_norm)")
        return opt
    else
        error("Linear optimization failed: $(result.error)")
    end
end

function optimize_rotation_step(individual_pcfs::Vector{IndividualPCF}, target_mu::Float64, N::Int, Ntheta::Int)
    println("\n" * "="^60)
    println("Step 4b: Rotation Optimization")
    println("="^60)
    println("  Target μ = $target_mu")

    # eDSLコマンド作成
    cmd = OptimizeRotation(individual_pcfs, target_mu, N, Ntheta)
    println("  Command: $cmd")

    # インタープリター実行
    @time result = interpret(cmd)

    if is_success(result)
        opt = unwrap(result)
        println("  Achieved R = $(opt.achieved)")
        println("  Frobenius norm ||C|| = $(opt.frobenius_norm)")
        return opt
    else
        error("Rotation optimization failed: $(result.error)")
    end
end

# =============================================================================
# ステップ5: Total PCF の計算
# =============================================================================

function compute_total_pcf_step(individual_pcfs::Vector{IndividualPCF}, C::Array{Float64,3}, N::Int, Ntheta::Int)
    println("\n" * "="^60)
    println("Step 5: Computing Total Phase Coupling Function")
    println("="^60)

    # eDSLコマンド作成
    cmd = ComputeTotalPCF(individual_pcfs, C, N, Ntheta)
    println("  Command: $cmd")

    # インタープリター実行
    @time result = interpret(cmd)

    if is_success(result)
        total_pcf = unwrap(result)
        println("  Γ₁ (∂Γ/∂φ at origin) = $(total_pcf.dGamma_dphi)")
        println("  Γ₂ (∂Γ/∂ψ at origin) = $(total_pcf.dGamma_dpsi)")

        # 安定性指標の計算
        Lambda = compute_linear_stability(total_pcf.dGamma_dphi, total_pcf.dGamma_dpsi)
        R = compute_rotation_characteristic(total_pcf.dGamma_dphi, total_pcf.dGamma_dpsi)
        println("  Linear stability Λ = $Lambda")
        println("  Rotation characteristic R = $R")

        return total_pcf
    else
        error("Total PCF computation failed: $(result.error)")
    end
end

# =============================================================================
# ステップ6: 位相ダイナミクスシミュレーション
# =============================================================================

function simulate_phase_dynamics_step(total_pcf::TotalPCF, epsilon::Float64, t_end::Float64)
    println("\n" * "="^60)
    println("Step 6: Simulating Phase Dynamics")
    println("="^60)
    println("  ε = $epsilon, t_end = $t_end")

    # 初期位相差
    phi1_0 = 0.1
    phi2_0 = -0.1

    # eDSLコマンド作成
    cmd = SimulatePhaseDynamics(total_pcf, phi1_0, phi2_0, epsilon, 0.0, t_end)
    println("  Command: $cmd")
    println("  Initial phases: φ₁ = $phi1_0, φ₂ = $phi2_0")

    # インタープリター実行
    @time result = interpret(cmd)

    if is_success(result)
        pts = unwrap(result)
        println("  Simulation completed")
        println("  Time points: $(length(pts.time))")
        println("  Final phases: φ₁ = $(pts.phase_diffs[1, end]), φ₂ = $(pts.phase_diffs[2, end])")
        return pts
    else
        error("Phase dynamics simulation failed: $(result.error)")
    end
end

# =============================================================================
# メイン実行関数
# =============================================================================

function run_full_pipeline(config)
    println("\n" * "#"^60)
    println("# HONetSync - Full Pipeline Execution")
    println("#"^60)
    println("Configuration:")
    for (k, v) in config
        println("  $k = $v")
    end

    N = config[:N]
    Ntheta = config[:Ntheta]
    t_relax = config[:t_relax]
    target_q = config[:target_q]
    target_mu = config[:target_mu]
    epsilon = config[:epsilon]
    t_end = config[:t_end]

    # ネットワークパラメータ作成
    params = NetworkParams(N)

    # パイプライン実行
    sps = compute_sps_step(params, Ntheta, t_relax)
    psf = compute_psf_step(sps)
    individual_pcfs = compute_individual_pcfs_step(sps, psf)

    results = Dict{Symbol, Any}(
        :sps => sps,
        :psf => psf,
        :individual_pcfs => individual_pcfs
    )

    # モードに応じた最適化
    if config[:mode] == :full || config[:mode] == :linear
        opt_linear = optimize_linear_step(individual_pcfs, target_q, N, Ntheta)
        results[:opt_linear] = opt_linear

        total_pcf_linear = compute_total_pcf_step(individual_pcfs, opt_linear.C_opt, N, Ntheta)
        results[:total_pcf_linear] = total_pcf_linear

        pts_linear = simulate_phase_dynamics_step(total_pcf_linear, epsilon, t_end)
        results[:pts_linear] = pts_linear
    end

    if config[:mode] == :full || config[:mode] == :rotation
        opt_rotation = optimize_rotation_step(individual_pcfs, target_mu, N, Ntheta)
        results[:opt_rotation] = opt_rotation

        total_pcf_rotation = compute_total_pcf_step(individual_pcfs, opt_rotation.C_opt, N, Ntheta)
        results[:total_pcf_rotation] = total_pcf_rotation

        pts_rotation = simulate_phase_dynamics_step(total_pcf_rotation, epsilon, t_end)
        results[:pts_rotation] = pts_rotation
    end

    # 結果の保存
    if config[:save]
        save_results(config[:output], results, config)
    end

    println("\n" * "#"^60)
    println("# Pipeline completed successfully!")
    println("#"^60)

    return results
end

function run_comparison(config)
    println("\n" * "#"^60)
    println("# HONetSync - Strategy Comparison")
    println("#"^60)

    N = config[:N]
    params = NetworkParams(N)

    # 比較実験実行
    results = run_comparison_experiment(params;
        target_q=config[:target_q],
        target_mu=config[:target_mu])

    # 結果表示
    println("\n" * "="^60)
    println("Comparison Results")
    println("="^60)

    strategies = [:uniform, :linear, :rotation]
    for strategy in strategies
        if haskey(results, strategy)
            r = results[strategy]
            println("\n[$strategy]")
            println("  Λ = $(r.Lambda)")
            println("  R = $(r.R)")
        end
    end

    return results
end

# =============================================================================
# 結果保存
# =============================================================================

function save_results(output_dir::AbstractString, results::Dict, config)
    println("\n" * "="^60)
    println("Saving Results")
    println("="^60)

    N = config[:N]
    Ntheta = config[:Ntheta]
    suffix = get(config, :suffix, "")

    # サブディレクトリの作成
    sps_dir = joinpath(output_dir, "sps")
    psf_dir = joinpath(output_dir, "psf")
    opt_dir = joinpath(output_dir, "opt")
    pcf_dir = joinpath(output_dir, "pcf")
    sim_dir = joinpath(output_dir, "sim")

    mkpath(sps_dir)
    mkpath(psf_dir)
    mkpath(opt_dir)
    mkpath(pcf_dir)
    mkpath(sim_dir)

    # SPS保存
    if haskey(results, :sps)
        filepath = joinpath(sps_dir, generate_filename(:sps, N, Ntheta; suffix=suffix))
        save_sps(filepath, results[:sps])
        println("  Saved: $filepath")
    end

    # PSF保存
    if haskey(results, :psf)
        filepath = joinpath(psf_dir, generate_filename(:psf, N, Ntheta; suffix=suffix))
        save_psf(filepath, results[:psf])
        println("  Saved: $filepath")
    end

    # 最適化結果保存
    if haskey(results, :opt_linear)
        filepath = joinpath(opt_dir, generate_filename(:optimization, N, Ntheta; coupling_type=:linear, suffix=suffix))
        save_optimization_result(filepath, results[:opt_linear])
        println("  Saved: $filepath")
    end

    if haskey(results, :opt_rotation)
        filepath = joinpath(opt_dir, generate_filename(:optimization, N, Ntheta; coupling_type=:rotation, suffix=suffix))
        save_optimization_result(filepath, results[:opt_rotation])
        println("  Saved: $filepath")
    end

    # Total PCF保存 (suffixを追加)
    suffix_str = isempty(suffix) ? "" : "_$(suffix)"
    if haskey(results, :total_pcf_linear)
        filepath = joinpath(pcf_dir, "total_pcf_linear_N$(N)_Nth$(Ntheta)$(suffix_str).jld2")
        save_total_pcf(filepath, results[:total_pcf_linear])
        println("  Saved: $filepath")
    end

    if haskey(results, :total_pcf_rotation)
        filepath = joinpath(pcf_dir, "total_pcf_rotation_N$(N)_Nth$(Ntheta)$(suffix_str).jld2")
        save_total_pcf(filepath, results[:total_pcf_rotation])
        println("  Saved: $filepath")
    end

    # 位相ダイナミクス保存
    if haskey(results, :pts_linear)
        filepath = joinpath(sim_dir, "phase_dynamics_linear_N$(N)_Nth$(Ntheta)$(suffix_str).jld2")
        save_phase_time_series(filepath, results[:pts_linear])
        println("  Saved: $filepath")
    end

    if haskey(results, :pts_rotation)
        filepath = joinpath(sim_dir, "phase_dynamics_rotation_N$(N)_Nth$(Ntheta)$(suffix_str).jld2")
        save_phase_time_series(filepath, results[:pts_rotation])
        println("  Saved: $filepath")
    end

    println("  All results saved to: $output_dir")
end

# =============================================================================
# エントリーポイント
# =============================================================================

function main(args=ARGS)
    config = parse_args(args)

    if config[:mode] == :compare
        return run_comparison(config)
    else
        return run_full_pipeline(config)
    end
end

# スクリプトとして実行された場合
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
