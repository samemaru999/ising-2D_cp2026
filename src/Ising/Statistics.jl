"""
MCサンプリングの蓄積と熱力学量の導出を行う純粋関数群。
§2.4 蓄積・導出（Logic, Pure）
"""

"""
    compute_thermodynamics(avg::ThermalAverages, N::Int, T::Float64) -> ThermodynamicQuantities

ThermalAverages から熱力学量を計算する純粋関数。
- N: 系のサイト数 (L×L)
- T: 温度
"""
function compute_thermodynamics(avg::ThermalAverages, N::Int, T::Float64)
    avg.n_samples > 0 || throw(ArgumentError("n_samples must be > 0"))
    N > 0 || throw(ArgumentError("N must be > 0"))
    T > 0 || throw(ArgumentError("T must be > 0"))

    n = avg.n_samples
    mean_e = avg.sum_e / n
    mean_e2 = avg.sum_e2 / n
    mean_abs_m = avg.sum_abs_m / n
    mean_m2 = avg.sum_m2 / n
    mean_m4 = avg.sum_m4 / n

    C = N / T^2 * (mean_e2 - mean_e^2)
    chi = N / T * (mean_m2 - mean_abs_m^2)
    U_L = 1.0 - mean_m4 / (3.0 * mean_m2^2)

    return ThermodynamicQuantities(mean_e, mean_abs_m, C, chi, U_L)
end
