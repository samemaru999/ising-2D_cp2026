#!/usr/bin/env julia
"""
    plot_3network_dynamics.jl - 3ネットワーク結合系の時系列プロット

各ネットワーク（A, B, C）のv成分を縦に並べて表示

使用方法:
    julia --project=. scripts/plot_3network_dynamics.jl
"""

using HONetSync
using GLMakie
using DifferentialEquations

println("3ネットワーク結合系シミュレーション")
println("="^50)

# パラメータ設定
N = 4  # 各ネットワークの振動子数
Ds = 2  # 状態次元 (u, v)

params = NetworkParams(N)
println("N = $N, FHN params: δ=$(params.fhn.delta), a=$(params.fhn.a), b=$(params.fhn.b)")

# 3ネットワークパラメータ
C = get_predefined_coupling(N, :linear)
epsilon = 2e-4
three_params = ThreeNetworkParams(params, C, epsilon)

# 初期条件（各ネットワークに少しずれた初期値）
function create_initial_state_3net(N, Ds)
    X_A = zeros(Ds * N)
    X_B = zeros(Ds * N)
    X_C = zeros(Ds * N)

    for i in 1:N
        # ネットワークA: 基準
        X_A[2i-1] = 1.0 + 0.1 * randn()
        X_A[2i] = 1.0 + 0.1 * randn()

        # ネットワークB: 少しずれた位相
        X_B[2i-1] = 0.8 + 0.1 * randn()
        X_B[2i] = 0.5 + 0.1 * randn()

        # ネットワークC: さらにずれた位相
        X_C[2i-1] = 0.5 + 0.1 * randn()
        X_C[2i] = -0.5 + 0.1 * randn()
    end

    return X_A, X_B, X_C
end

X_A0, X_B0, X_C0 = create_initial_state_3net(N, Ds)
X0 = vcat(X_A0, X_B0, X_C0)

# ODEシミュレーション
function dynamics!(dX, X, p, t)
    N = p.network.N
    Ds = 2

    X_A = X[1:Ds*N]
    X_B = X[Ds*N+1:2*Ds*N]
    X_C = X[2*Ds*N+1:3*Ds*N]

    state = ThreeNetworkState(X_A, X_B, X_C)
    dX_A, dX_B, dX_C = compute_three_network_dynamics(state, p)

    dX[1:Ds*N] = dX_A
    dX[Ds*N+1:2*Ds*N] = dX_B
    dX[2*Ds*N+1:3*Ds*N] = dX_C
end

println("シミュレーション実行中...")
t_end = 250.0
prob = ODEProblem(dynamics!, X0, (0.0, t_end), three_params)
sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-8, saveat=0.5)

println("シミュレーション完了: $(length(sol.t)) points")

# プロット作成
println("プロット作成中...")

# 色の定義
colors = [:black, :red, :green, :blue, :cyan, :magenta, :orange, :purple, :brown, :gray]

fig = Figure(size=(900, 700), backgroundcolor=:white)

network_labels = ["Network A", "Network B", "Network C"]

for (net_idx, label) in enumerate(network_labels)
    ax = Axis(fig[net_idx, 1],
        title = label,
        xlabel = net_idx == 3 ? "t" : "",
        ylabel = "vᵢ",
        xticklabelsvisible = net_idx == 3,
    )

    # 各振動子のv成分をプロット
    for i in 1:N
        # v成分のインデックス
        v_idx = Ds * N * (net_idx - 1) + 2 * i

        v_series = [u[v_idx] for u in sol.u]

        lines!(ax, sol.t, v_series,
            color = colors[i],
            linewidth = 1.0,
            label = "#$i")
    end

    ylims!(ax, -4, 4)
    hlines!(ax, [0], color=(:black, 0.3), linestyle=:dash, linewidth=0.5)
end

# 凡例
Legend(fig[1:3, 2], fig.content[1], framevisible=false, labelsize=12)

display(fig)

println("\n3ネットワーク v-component プロット完了!")
println("Press Enter to close...")
readline()
