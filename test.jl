# %%
using Pkg
Pkg.status()
# %%
#using ModelingToolkit: t_nounits as t, D_nounits as D
#import ModelingToolkit as MTK
import DifferentialEquations as DiEq
using Enzyme, LinearAlgebra, Statistics
using CairoMakie, ProgressMeter, JLD2, StaticArrays, Graphs
# %%
include("src/myconst.jl")
include("src/mydyn.jl")

abstract type Result{T,E} end
struct Ok{T,E} <: Result{T,E}
    value::T
end
struct Err{T,E} <: Result{T,E}
    error::E
end

# Result型のパイプライン演算子対応
Base.:|>(r::Ok, f::Function) = f(r.value)
Base.:|>(r::Err, f::Function) = r
# %%
# Mylib
#include("src/find_sps_net.jl")
#include("src/phase_sensitivity.jl")
#include("src/phase_coupling_func.jl")
#include("src/differential_evolution.jl")
#using .PhaseCouplingFunc
#
# %%
evol(x::AllNetworkState, p, t) = begin
    #du = AllNetworkState(undef) → いける

    dx = map(1:length(x)) do target_net_idx
        # 各ネットワークの状態 du[target_net_idx]を更新

        # 個々の振動子の固有ダイナミクスを計算
        dyn_f = map(enumerate(x[target_net_idx])) do (i, state_i)
            oscillator_dynamics(state_i, p.oscillator_params, i)
        end

        #同一ネットワーク内相互作用
        intra_coup_g = intraNet_dynamics(x[target_net_idx],
            p.intra_network_topologies)
        #network_dynamics .+= p.intra_coupling_strength .* intra_coupling

        # ネットワーク間相互作用
        inter_coup_h = map(1:length(x)) do source_net_idx
            if source_net_idx == target_net_idx
                return [zeros(2) for i = 1:N_osc]
            else
                compute_inter_network_coupling(
                    x[source_net_idx],
                    x[target_net_idx],
                    p.inter_network_topologies
                )
            end
            #network_dynamics .+= p.inter_coupling_strength .* inter_coupling
        end |> sum
        network_dynamics = dyn_f + intra_coup_g + inter_coup_h
        return network_dynamics
    end
    # 各ネットワークの状態を更新
    return reduce(vcat, dx)
end

# %%
