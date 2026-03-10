#!/usr/bin/env julia
"""
    visualize_coupling_tensor.jl - 高次結合テンソルの可視化

3ネットワーク間の高次結合を三角形で表現
配列形式: C[i,j,k] = [networkB_node, networkC_node, networkA_node]
"""

using GLMakie
using Statistics

"""
    visualize_higher_order_coupling(C::Array{Float64,3};
                                     threshold::Float64=1.0,
                                     filename::String="coupling_tensor.png")

高次結合テンソルを三角形で可視化

# Arguments
- `C`: 3次元配列 [networkB, networkC, networkA]
- `threshold`: この値以上の結合のみ表示
- `filename`: 出力ファイル名
"""
function visualize_higher_order_coupling(C::Array{Float64,3};
                                          threshold::Float64=1.0,
                                          filename::String="coupling_tensor.png",
                                          figsize=(1200, 900))
    N = size(C, 1)

    # ノード位置の定義
    # Network A: 上部（横一列）
    # Network B: 左側（45度右回り = 左上→右下）
    # Network C: 右側（45度左回り = 右上→左下）

    spacing = 1.5  # ノード間隔
    angle = π / 4  # 45度

    # Network A: 上部（横一列、中央揃え）
    pos_A = [(i * spacing, 10.0) for i in 1:N]

    # Network B: 左側、45度時計回り（左上→右下）
    B_start_x = 0.5
    B_start_y = 6.0
    dx_B = spacing * cos(angle)  # 右方向
    dy_B = spacing * sin(angle)  # 下方向
    pos_B = [(B_start_x + (i-1)*dx_B, B_start_y - (i-1)*dy_B) for i in 1:N]

    # Network C: 右側、45度反時計回り（右上→左下）
    C_start_x = (N + 1) * spacing + 1.0
    C_start_y = 6.0
    dx_C = spacing * cos(angle)  # 左方向
    dy_C = spacing * sin(angle)  # 下方向
    pos_C = [(C_start_x - (i-1)*dx_C, C_start_y - (i-1)*dy_C) for i in 1:N]

    # Figure作成
    fig = Figure(size=figsize, backgroundcolor=:white)
    ax = Axis(fig[1, 1],
              aspect=DataAspect(),
              backgroundcolor=:white)
    hidedecorations!(ax)
    hidespines!(ax)

    # カラーマップ（正負で色分け）
    max_val = maximum(abs.(C))

    # 三角形を描画（0でない値のみ）
    triangles_pos = []  # 正の値
    triangles_neg = []  # 負の値
    values_pos = Float64[]
    values_neg = Float64[]

    for k in 1:N  # Network A
        for i in 1:N  # Network B
            for j in 1:N  # Network C
                val = C[i, j, k]
                if abs(val) >= threshold
                    # 三角形の頂点
                    pA = pos_A[k]
                    pB = pos_B[i]
                    pC = pos_C[j]

                    if val > 0
                        push!(triangles_pos, (pA, pB, pC))
                        push!(values_pos, val)
                    else
                        push!(triangles_neg, (pA, pB, pC))
                        push!(values_neg, abs(val))
                    end
                end
            end
        end
    end

    # 三角形を描画（値が小さいものから先に描画）
    # 正の値（青系）- 塗りつぶし
    if !isempty(triangles_pos)
        sorted_idx = sortperm(values_pos)
        for idx in sorted_idx
            tri = triangles_pos[idx]
            val = values_pos[idx]
            alpha = 0.15 + 0.35 * (val / max_val)  # 塗りつぶし用（薄め）

            # 三角形を塗りつぶし
            vertices = Point2f[(tri[1][1], tri[1][2]), (tri[2][1], tri[2][2]), (tri[3][1], tri[3][2])]
            poly!(ax, vertices, color=(:royalblue, alpha), strokecolor=(:royalblue, alpha + 0.2), strokewidth=1)
        end
    end

    # 負の値（緑系）- 塗りつぶし
    if !isempty(triangles_neg)
        sorted_idx = sortperm(values_neg)
        for idx in sorted_idx
            tri = triangles_neg[idx]
            val = values_neg[idx]
            alpha = 0.15 + 0.35 * (val / max_val)  # 塗りつぶし用（薄め）

            # 三角形を塗りつぶし
            vertices = Point2f[(tri[1][1], tri[1][2]), (tri[2][1], tri[2][2]), (tri[3][1], tri[3][2])]
            poly!(ax, vertices, color=(:seagreen, alpha), strokecolor=(:seagreen, alpha + 0.2), strokewidth=1)
        end
    end

    # ノードを描画
    node_size = 30

    # Network A（上部）
    for (i, p) in enumerate(pos_A)
        scatter!(ax, [p[1]], [p[2]], markersize=node_size, color=:white,
                 strokecolor=:royalblue, strokewidth=2)
        text!(ax, p[1], p[2], text=string(i), align=(:center, :center),
              fontsize=16, color=:royalblue)
    end

    # Network B（左下）
    for (i, p) in enumerate(pos_B)
        scatter!(ax, [p[1]], [p[2]], markersize=node_size, color=:white,
                 strokecolor=:royalblue, strokewidth=2)
        text!(ax, p[1], p[2], text=string(i), align=(:center, :center),
              fontsize=16, color=:royalblue)
    end

    # Network C（右下）
    for (i, p) in enumerate(pos_C)
        scatter!(ax, [p[1]], [p[2]], markersize=node_size, color=:white,
                 strokecolor=:royalblue, strokewidth=2)
        text!(ax, p[1], p[2], text=string(i), align=(:center, :center),
              fontsize=16, color=:royalblue)
    end

    # ネットワークラベル
    text!(ax, mean([p[1] for p in pos_A]), 11.5, text="Network A",
          align=(:center, :bottom), fontsize=20, color=:black, font=:bold)
    text!(ax, pos_B[1][1] - 1.5, pos_B[1][2], text="Network B",
          align=(:right, :center), fontsize=20, color=:black, font=:bold)
    text!(ax, pos_C[1][1] + 1.5, pos_C[1][2], text="Network C",
          align=(:left, :center), fontsize=20, color=:black, font=:bold)

    # 軸の余白を設定
    xlims!(ax, -3.0, (N + 1) * spacing + 4.0)
    ylims!(ax, pos_B[N][2] - 2.0, 13.0)

    # 凡例
    Legend(fig[1, 2],
           [PolyElement(color=(:royalblue, 0.5), strokecolor=:royalblue, strokewidth=2),
            PolyElement(color=(:seagreen, 0.5), strokecolor=:seagreen, strokewidth=2)],
           ["Positive coupling", "Negative coupling"],
           framevisible=false)

    # 保存
    save(filename, fig, px_per_unit=2)
    println("Saved: $filename")

    return fig
end

# 使用例
function demo()
    # サンプルデータ（ユーザー提供のテンソル）
    C = zeros(4, 4, 4)

    # [:, :, 1]
    C[:, :, 1] = [
        0.0      -8.4675  30.7509   13.3751;
        8.4675    0.0     73.3927   63.1924;
       -30.7509  -73.3927   0.0       9.58208;
       -13.3751  -63.1924  -9.58208   0.0
    ]

    # [:, :, 2]
    C[:, :, 2] = [
        0.0       4.96969  -44.2086  -74.4692;
       -4.96969   0.0        4.6082   -8.64349;
       44.2086   -4.6082     0.0      10.2945;
       74.4692    8.64349  -10.2945    0.0
    ]

    # [:, :, 3]
    C[:, :, 3] = [
        0.0     -11.5582   -12.4097   -11.4672;
       11.5582    0.0        5.40493    5.0409;
       12.4097   -5.40493    0.0       -7.73085;
       11.4672   -5.0409     7.73085    0.0
    ]

    # [:, :, 4]
    C[:, :, 4] = [
        0.0       -7.53523  12.1312   2.03443;
        7.53523    0.0      11.7236   2.85113;
       -12.1312   -11.7236    0.0      9.61348;
       -2.03443   -2.85113  -9.61348  0.0
    ]

    # 閾値以上の結合のみ表示
    visualize_higher_order_coupling(C; threshold=10.0,
                                     filename="data/coupling_tensor_viz.png")
end

# スクリプトとして実行
if abspath(PROGRAM_FILE) == @__FILE__
    demo()
end
