#!/usr/bin/env julia
"""
    example_edsl.jl - eDSLパターンの使用例

eDSL (embedded Domain Specific Language) パターンの理解を助けるサンプルコード
実際のODE計算は行わず、APIの使用方法を示す

使用方法:
    julia --project=. scripts/example_edsl.jl
"""

using HONetSync

println("="^60)
println("HONetSync eDSL Pattern Examples")
println("="^60)

# =============================================================================
# 1. ドメイン型の作成
# =============================================================================
println("\n[1] Domain Types (Pure Layer)")
println("-"^40)

# FHNパラメータ
fhn = FHNParams()  # デフォルト値
println("FHNParams: δ=$(fhn.delta), a=$(fhn.a), b=$(fhn.b)")

# ネットワークパラメータ
params = NetworkParams(4)
println("NetworkParams: N=$(params.N), topology=$(params.topology)")

# 3ネットワークパラメータ
three_params = ThreeNetworkParams(params)
println("ThreeNetworkParams: ε=$(three_params.epsilon)")

# =============================================================================
# 2. eDSLコマンドの作成
# =============================================================================
println("\n[2] eDSL Commands (Pure Layer)")
println("-"^40)

# 各コマンドは計算の「記述」であり、実行ではない
cmd1 = ComputeSPS(params)
cmd2 = ComputeSPS(params, 101, 1000.0)  # カスタムパラメータ

println("Command 1: $cmd1")
println("Command 2: $cmd2")

# コマンドはイミュータブルなデータ構造
println("\nCommand properties:")
println("  params.N = $(cmd1.params.N)")
println("  Ntheta = $(cmd1.Ntheta)")
println("  t_relax = $(cmd1.t_relax)")

# =============================================================================
# 3. コマンドスクリプト（コマンドのリスト）
# =============================================================================
println("\n[3] Command Scripts")
println("-"^40)

# スクリプトはコマンドのベクトル
script = ComputationScript([
    ComputeSPS(params),
    # ComputePSF, ComputeIndividualPCF などを追加可能
])

println("Script with $(length(script.commands)) command(s)")

# =============================================================================
# 4. 結果型 (Success / Failure)
# =============================================================================
println("\n[4] Result Types (Either Pattern)")
println("-"^40)

# Success の例
success_result = Success(42)
println("Success result: $(unwrap(success_result))")
println("  is_success: $(is_success(success_result))")

# Failure の例
failure_result = Failure{Int}("Something went wrong", Dict("code" => 500))
println("Failure result: $(failure_result.error)")
println("  is_success: $(is_success(failure_result))")

# unwrap_or でデフォルト値を指定
println("  unwrap_or: $(unwrap_or(failure_result, -1))")

# map_result で変換
mapped = map_result(x -> x * 2, success_result)
println("  mapped success: $(unwrap(mapped))")

# =============================================================================
# 5. インタープリターパターン
# =============================================================================
println("\n[5] Interpreter Pattern (Impure Layer)")
println("-"^40)

println("""
eDSLの実行フロー:

  [Pure Layer]                    [Impure Layer]
  ─────────────                   ──────────────
  コマンド作成                     interpret(cmd)
       │                              │
       ▼                              ▼
  ComputeSPS(params)  ──────▶  ODE Solver実行
       │                              │
       ▼                              ▼
  (データ構造)                   (実際の計算)
                                      │
                                      ▼
                              ComputationResult
                              Success or Failure
""")

# interpret関数の多重ディスパッチ
println("interpret() dispatches on command type:")
println("  interpret(::ComputeSPS) → SimulationInterpreter.jl")
println("  interpret(::ComputePSF) → SimulationInterpreter.jl")
println("  interpret(::ComputeIndividualPCF) → PhaseInterpreter.jl")
println("  interpret(::OptimizeLinearStability) → PhaseInterpreter.jl")

# =============================================================================
# 6. 事前定義された結合行列
# =============================================================================
println("\n[6] Predefined Coupling Matrices")
println("-"^40)

for coupling_type in [:uniform, :linear, :rotation]
    C = get_predefined_coupling(coupling_type, 4)
    norm = compute_coupling_frobenius_norm(C)
    println("$coupling_type (N=4): ||C|| = $(round(norm, digits=6))")
end

# =============================================================================
# 7. Pure関数の例
# =============================================================================
println("\n[7] Pure Functions (Logic Layer)")
println("-"^40)

# FHNダイナミクス（純粋関数）
u, v = 0.5, 0.3
I = 0.4
fhn_params = FHNParams()

du, dv = fhn_single_dynamics(u, v, I, fhn_params)
println("FHN dynamics at (u=$u, v=$v):")
println("  du/dt = $du")
println("  dv/dt = $dv")

# 位相縮約の安定性計算（純粋関数）
Gamma1, Gamma2 = 0.1, 0.05
Lambda = compute_linear_stability(Gamma1, Gamma2)
R = compute_rotation_characteristic(Gamma1, Gamma2)
println("\nPhase reduction (Γ₁=$Gamma1, Γ₂=$Gamma2):")
println("  Linear stability Λ = $Lambda")
println("  Rotation R = $R")

# =============================================================================
# 完了
# =============================================================================
println("\n" * "="^60)
println("eDSL Pattern Summary")
println("="^60)
println("""
Key Concepts:
  1. Commands are data structures describing computations
  2. interpret() executes commands (impure)
  3. Results are wrapped in Success/Failure
  4. Pure functions have no side effects
  5. Separation of concerns: Domain → Logic → Interpreters → Runtime
""")
