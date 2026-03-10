#!/usr/bin/env julia
"""
    visualize_user_tensor.jl - ユーザー提供テンソルの可視化
"""

# visualize_coupling_tensor.jl から関数をインクルード
include("visualize_coupling_tensor.jl")

function main()
    # ユーザー提供の4×4×4テンソル
    C = zeros(4, 4, 4)

    # [:, :, 1]
    C[:, :, 1] = [
        110.201   66.5065   24.9034   19.1023;
         66.5065  36.3752    9.22956   1.6673;
         24.9034   9.22956   2.26415  -3.03168;
         19.1023   1.6673   -3.03168  -7.65055
    ]

    # [:, :, 2]
    C[:, :, 2] = [
          7.48441  28.9589  34.6101   27.4172;
         28.9589   54.8922  18.8754   16.4867;
         34.6101   18.8754   1.21353  -1.82615;
         27.4172   16.4867  -1.82615  -1.84915
    ]

    # [:, :, 3]
    C[:, :, 3] = [
         -3.27374  -11.1674     8.84357    2.35428;
        -11.1674    -9.39606   16.866     10.0738;
          8.84357   16.866    -24.4425   -13.8933;
          2.35428   10.0738   -13.8933    -4.62289
    ]

    # [:, :, 4]
    C[:, :, 4] = [
          6.90885   10.1421    0.577644   0.463651;
         10.1421     3.52756  -5.74043   -6.09644;
          0.577644  -5.74043   7.11396    0.254868;
          0.463651  -6.09644   0.254868  -4.55575
    ]

    # テンソルの統計情報を表示
    println("=== Tensor Statistics ===")
    println("Size: ", size(C))
    println("Min: ", minimum(C))
    println("Max: ", maximum(C))
    println("Mean: ", sum(C) / length(C))
    println("Non-zero elements: ", count(x -> abs(x) > 1e-6, C))
    println()

    # 値の分布を確認
    println("=== Value Distribution ===")
    for thresh in [5.0, 10.0, 20.0, 50.0]
        n_above = count(x -> abs(x) >= thresh, C)
        println("  |C| >= $thresh: $n_above elements")
    end
    println()

    # 可視化（閾値=10.0）
    println("Generating visualization with threshold=10.0...")
    visualize_higher_order_coupling(C;
        threshold=10.0,
        filename="data/user_tensor_viz.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
