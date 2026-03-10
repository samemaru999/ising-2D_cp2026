"""
    Runtime/IO/IOHandler.jl - データ入出力ハンドラー [Impure]

JLD2形式でのデータ保存・読み込みを行う

ファイル命名規則:
{データ種別}_N{振動子数}_Nθ{位相分割数}_{オプション}_{日付}.jld2
"""

using JLD2
using Dates
using Printf

# =============================================================================
# ファイル命名
# =============================================================================

"""
    generate_filename(data_type::Symbol, N::Int, Ntheta::Int;
                       topology_type::Symbol=:none,
                       coupling_type::Symbol=:none, target_value::Float64=0.0,
                       suffix::String="")
        -> String

命名規則に従ったファイル名を生成

# Arguments
- `topology_type::Symbol`: ネットワーク結合構造 (:full, :tree, :none)
- `coupling_type::Symbol`: 最適化タイプ (:linear, :rotation, :none)

suffixがタイムスタンプ形式（yyyymmdd_HHMMSS）の場合は日付を追加しない
"""
function generate_filename(data_type::Symbol, N::Int, Ntheta::Int;
                            topology_type::Symbol=:none,
                            coupling_type::Symbol=:none,
                            target_value::Float64=0.0,
                            suffix::String="")
    base = "$(data_type)_N$(N)_Nth$(Ntheta)"

    if topology_type != :none
        base *= "_$(topology_type)"
    end

    if coupling_type != :none
        base *= "_$(coupling_type)"
    end

    if target_value != 0.0
        base *= @sprintf("_%.2f", abs(target_value))
    end

    if !isempty(suffix)
        base *= "_$(suffix)"
        # suffixがタイムスタンプ形式（8桁の数字で始まる）場合は日付を追加しない
        if occursin(r"^\d{8}", suffix)
            return "$(base).jld2"
        end
    end

    # suffixがない or タイムスタンプ形式でない場合は日付を追加
    date_str = Dates.format(today(), "yyyymmdd")
    return "$(base)_$(date_str).jld2"
end

# =============================================================================
# SaveData コマンド実行
# =============================================================================

"""
    interpret(cmd::SaveData) -> ComputationResult{Nothing}

データ保存コマンドを実行
"""
function interpret(cmd::SaveData)
    try
        data = cmd.data
        filename = cmd.filename
        format = cmd.format

        if format == :jld2
            save_jld2(filename, data)
        elseif format == :csv
            save_csv(filename, data)
        else
            return Failure{Nothing}("Unsupported format: $format")
        end

        @info "Data saved to $filename"
        return Success(nothing)

    catch e
        return Failure{Nothing}("Save failed: $e", Dict("exception" => e))
    end
end

"""
    interpret(cmd::LoadData) -> ComputationResult

データ読み込みコマンドを実行
"""
function interpret(cmd::LoadData)
    try
        filename = cmd.filename
        expected_type = cmd.expected_type

        if !isfile(filename)
            return Failure{expected_type}("File not found: $filename")
        end

        data = load_jld2(filename, expected_type)
        return Success(data)

    catch e
        return Failure{cmd.expected_type}("Load failed: $e", Dict("exception" => e))
    end
end

# =============================================================================
# JLD2 保存関数
# =============================================================================

"""
    save_jld2(filename, data)

JLD2形式でデータを保存
"""
function save_jld2(filename::String, data)
    # ディレクトリが存在しなければ作成
    dir = dirname(filename)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end

    jldsave(filename; data=data, metadata=create_metadata())
end

"""
    save_sps(filename, sps::StablePeriodicSolution)

安定周期解を保存
"""
function save_sps(filename::String, sps::StablePeriodicSolution)
    jldsave(filename;
        type = "StablePeriodicSolution",
        version = "2.1",
        N = sps.N,
        Ntheta = sps.Ntheta,
        T = sps.T,
        omega = sps.omega,
        Xs = sps.Xs,
        params = serialize_params(sps.params),
        X0_initial = sps.X0_initial,  # 緩和前の初期値
        metadata = create_metadata()
    )
    @info "SPS saved to $filename"
end

"""
    save_psf(filename, psf::PhaseSensitivityFunction)

位相感受関数を保存
"""
function save_psf(filename::String, psf::PhaseSensitivityFunction)
    jldsave(filename;
        type = "PhaseSensitivityFunction",
        version = "2.0",
        N = psf.N,
        Ntheta = psf.Ntheta,
        Q = psf.Q,
        T = psf.T,
        omega = psf.omega,
        metadata = create_metadata()
    )
    @info "PSF saved to $filename"
end

"""
    save_total_pcf(filename, pcf::TotalPCF)

全体位相結合関数を保存
"""
function save_total_pcf(filename::String, pcf::TotalPCF)
    jldsave(filename;
        type = "TotalPCF",
        version = "2.0",
        N = pcf.N,
        Ntheta = pcf.Ntheta,
        Gamma = pcf.Gamma,
        C = pcf.C,
        dGamma_dphi = pcf.dGamma_dphi,
        dGamma_dpsi = pcf.dGamma_dpsi,
        metadata = create_metadata()
    )
    @info "TotalPCF saved to $filename"
end

"""
    save_optimization_result(filename, result::OptimizationResult)

最適化結果を保存
"""
function save_optimization_result(filename::String, result::OptimizationResult)
    jldsave(filename;
        type = "OptimizationResult",
        version = "2.0",
        N = result.N,
        C_opt = result.C_opt,
        opt_type = result.opt_type,
        target = result.target,
        achieved = result.achieved,
        frobenius_norm = result.frobenius_norm,
        metadata = create_metadata()
    )
    @info "OptimizationResult saved to $filename"
end

"""
    save_phase_time_series(filename, pts::PhaseTimeSeries)

位相時系列を保存
"""
function save_phase_time_series(filename::String, pts::PhaseTimeSeries)
    jldsave(filename;
        type = "PhaseTimeSeries",
        version = "2.0",
        time = pts.time,
        phases = pts.phases,
        phase_diffs = pts.phase_diffs,
        metadata = create_metadata()
    )
    @info "PhaseTimeSeries saved to $filename"
end

"""
    save_pipeline_results(dirname, results::NamedTuple)

パイプライン実行結果を一括保存
"""
function save_pipeline_results(dirname::String, results::NamedTuple)
    if !isdir(dirname)
        mkpath(dirname)
    end

    N = results.sps.N
    Ntheta = results.sps.Ntheta

    # 各結果を保存
    save_sps(joinpath(dirname, generate_filename(:sps, N, Ntheta)), results.sps)
    save_psf(joinpath(dirname, generate_filename(:psf, N, Ntheta)), results.psf)
    save_total_pcf(joinpath(dirname, generate_filename(:pcf, N, Ntheta;
        coupling_type=results.coupling_type)), results.total_pcf)

    if results.optimization !== nothing
        save_optimization_result(joinpath(dirname,
            generate_filename(:opt, N, Ntheta; coupling_type=results.coupling_type)),
            results.optimization)
    end

    @info "Pipeline results saved to $dirname"
end

# =============================================================================
# JLD2 読み込み関数
# =============================================================================

"""
    load_jld2(filename, ::Type{T}) where T

JLD2形式からデータを読み込み
"""
function load_jld2(filename::String, ::Type{T}) where T
    data = JLD2.load(filename)
    return reconstruct(T, data)
end

"""
    load_sps(filename) -> StablePeriodicSolution

安定周期解を読み込み
"""
function load_sps(filename::String)
    data = JLD2.load(filename)

    params = deserialize_params(data["params"])

    return StablePeriodicSolution(
        data["N"],
        data["Ntheta"],
        data["T"],
        data["omega"],
        data["Xs"],
        params
    )
end

"""
    load_psf(filename) -> PhaseSensitivityFunction

位相感受関数を読み込み
"""
function load_psf(filename::String)
    data = JLD2.load(filename)

    return PhaseSensitivityFunction(
        data["N"],
        data["Ntheta"],
        data["Q"],
        data["T"],
        data["omega"]
    )
end

"""
    load_total_pcf(filename) -> TotalPCF

全体位相結合関数を読み込み
"""
function load_total_pcf(filename::String)
    data = JLD2.load(filename)

    return TotalPCF(
        data["N"],
        data["Ntheta"],
        data["Gamma"],
        data["C"],
        data["dGamma_dphi"],
        data["dGamma_dpsi"]
    )
end

"""
    load_optimization_result(filename) -> OptimizationResult

最適化結果を読み込み
"""
function load_optimization_result(filename::String)
    data = JLD2.load(filename)

    return OptimizationResult(
        data["N"],
        data["C_opt"],
        data["opt_type"],
        data["target"],
        data["achieved"],
        data["frobenius_norm"]
    )
end

# =============================================================================
# シリアライゼーション補助関数
# =============================================================================

"""
    create_metadata() -> Dict{String, Any}

メタデータを生成
"""
function create_metadata()
    return Dict{String, Any}(
        "created_at" => Dates.now(),
        "julia_version" => string(VERSION),
        "package_version" => "2.0.0"
    )
end

"""
    serialize_params(params::NetworkParams) -> Dict

NetworkParamsをシリアライズ可能な辞書に変換
"""
function serialize_params(params::NetworkParams)
    return Dict{String, Any}(
        "N" => params.N,
        "fhn" => Dict(
            "delta" => params.fhn.delta,
            "a" => params.fhn.a,
            "b" => params.fhn.b
        ),
        "I" => params.I,
        "K" => params.K,
        "topology_type" => string(typeof(params.topology))
    )
end

"""
    deserialize_params(data::Dict) -> NetworkParams

辞書からNetworkParamsを再構築
"""
function deserialize_params(data::Dict)
    N = data["N"]
    fhn_data = data["fhn"]
    fhn = FHNParams(fhn_data["delta"], fhn_data["a"], fhn_data["b"])

    return NetworkParams(N, fhn, data["I"], data["K"], FullCoupling())
end

"""
    reconstruct(::Type{T}, data::Dict) where T

辞書からドメイン型を再構築（汎用）
"""
function reconstruct(::Type{StablePeriodicSolution}, data::Dict)
    return load_sps_from_dict(data)
end

function reconstruct(::Type{PhaseSensitivityFunction}, data::Dict)
    return PhaseSensitivityFunction(
        data["N"], data["Ntheta"], data["Q"], data["T"], data["omega"]
    )
end

function reconstruct(::Type{TotalPCF}, data::Dict)
    return TotalPCF(
        data["N"], data["Ntheta"], data["Gamma"], data["C"],
        data["dGamma_dphi"], data["dGamma_dpsi"]
    )
end

function reconstruct(::Type{OptimizationResult}, data::Dict)
    return OptimizationResult(
        data["N"], data["C_opt"], data["opt_type"],
        data["target"], data["achieved"], data["frobenius_norm"]
    )
end

function load_sps_from_dict(data::Dict)
    params = deserialize_params(data["params"])
    return StablePeriodicSolution(
        data["N"], data["Ntheta"], data["T"], data["omega"],
        data["Xs"], params
    )
end

# =============================================================================
# CSV出力（簡易データ用）
# =============================================================================

"""
    save_csv(filename, data)

CSV形式でデータを保存（主に位相時系列用）
"""
function save_csv(filename::String, data)
    dir = dirname(filename)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end

    if data isa PhaseTimeSeries
        save_phase_time_series_csv(filename, data)
    else
        error("CSV export not implemented for $(typeof(data))")
    end
end

"""
    save_phase_time_series_csv(filename, pts::PhaseTimeSeries)

位相時系列をCSV形式で保存
"""
function save_phase_time_series_csv(filename::String, pts::PhaseTimeSeries)
    open(filename, "w") do io
        # ヘッダー
        println(io, "time,phase_A,phase_B,phase_C,phi1,phi2")

        # データ
        for t in 1:length(pts.time)
            @printf(io, "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                pts.time[t],
                pts.phases[1, t], pts.phases[2, t], pts.phases[3, t],
                pts.phase_diffs[1, t], pts.phase_diffs[2, t])
        end
    end
end

# =============================================================================
# ファイル一覧・管理
# =============================================================================

"""
    list_data_files(directory::String; pattern::Regex=r"\\.jld2\$")
        -> Vector{String}

指定ディレクトリ内のデータファイル一覧を取得
"""
function list_data_files(directory::String; pattern::Regex=r"\.jld2$")
    if !isdir(directory)
        return String[]
    end

    files = readdir(directory; join=true)
    return filter(f -> occursin(pattern, f), files)
end

"""
    get_file_info(filename::String) -> NamedTuple

ファイルのメタデータを取得
"""
function get_file_info(filename::String)
    if !isfile(filename)
        error("File not found: $filename")
    end

    data = JLD2.load(filename)

    return (
        filename = basename(filename),
        type = get(data, "type", "unknown"),
        version = get(data, "version", "unknown"),
        created_at = get(get(data, "metadata", Dict()), "created_at", nothing),
        N = get(data, "N", nothing),
        Ntheta = get(data, "Ntheta", nothing)
    )
end

