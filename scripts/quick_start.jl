#!/usr/bin/env julia
"""
    quick_start.jl - HONetSync クイックスタート

最小限のコードでeDSLパイプラインを実行する例

使用方法:
    julia --project=. scripts/quick_start.jl
"""

using HONetSync

println("HONetSync Quick Start Example")
println("="^50)

# =============================================================================
# Step 1: パラメータ設定
# =============================================================================
println("\n[Step 1] Creating network parameters...")

N = 4  # 振動子数
params = NetworkParams(N)

println("  N = $N oscillators")
println("  δ = $(params.fhn.delta), a = $(params.fhn.a), b = $(params.fhn.b)")

# =============================================================================
# Step 2: 安定周期解 (SPS) の計算
# =============================================================================
println("\n[Step 2] Computing Stable Periodic Solution...")

# eDSLコマンド: ComputeSPS
cmd_sps = ComputeSPS(params, 51, 100.0)  # 短いパラメータで高速化
result_sps = interpret(cmd_sps)

if is_success(result_sps)
    sps = unwrap(result_sps)
    println("  Success! Period T = $(round(sps.T, digits=4))")
else
    error("Failed: $(result_sps.error)")
end

# =============================================================================
# Step 3: 位相感受関数 (PSF) の計算
# =============================================================================
println("\n[Step 3] Computing Phase Sensitivity Function...")

# eDSLコマンド: ComputePSF
cmd_psf = ComputePSF(sps, 5)  # 少ない反復回数で高速化
result_psf = interpret(cmd_psf)

if is_success(result_psf)
    psf = unwrap(result_psf)
    println("  Success! Q matrix size = $(size(psf.Q))")
else
    error("Failed: $(result_psf.error)")
end

# =============================================================================
# Step 4: 個別PCFの計算 (一部のみ)
# =============================================================================
println("\n[Step 4] Computing Individual PCFs (subset)...")

individual_pcfs = IndividualPCF[]

# デモ用に一部の組み合わせのみ計算
for i in 1:2, j in 1:2, k in 1:2
    cmd = ComputeIndividualPCF(sps, psf, i, j, k)
    result = interpret(cmd)
    if is_success(result)
        push!(individual_pcfs, unwrap(result))
    end
end

println("  Computed $(length(individual_pcfs)) individual PCFs")

# =============================================================================
# Step 5: 事前定義された結合を使用
# =============================================================================
println("\n[Step 5] Using predefined coupling tensor...")

# 事前定義された線形最適化結合を取得
C_linear = get_predefined_coupling(N, :linear)
println("  Using linear-optimized coupling (N=$N)")
println("  ||C|| = $(round(compute_coupling_frobenius_norm(C_linear), digits=6))")

# =============================================================================
# 結果サマリー
# =============================================================================
println("\n" * "="^50)
println("Quick Start Completed!")
println("="^50)
println("""
Next steps:
  1. Run full pipeline: julia --project=. scripts/main.jl
  2. Run tests: julia --project=. -e 'using Pkg; Pkg.test()'
  3. See README.md for detailed documentation
""")
