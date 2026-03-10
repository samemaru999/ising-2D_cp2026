"""
シミュレーション結果のJLD2保存・読み込み（Impure）。
"""

using JLD2

"""
    save_result(filename::String, result::TemperatureSweepResult) -> Nothing

温度掃引結果をJLD2形式で保存する。
"""
function save_result(filename::String, result::TemperatureSweepResult)
    jldsave(filename;
        params = result.params,
        config = result.config,
        temperatures = result.temperatures,
        results = result.results
    )
    return nothing
end

"""
    save_result(filename::String, result::SingleTemperatureResult) -> Nothing

単一温度結果をJLD2形式で保存する。
"""
function save_result(filename::String, result::SingleTemperatureResult)
    jldsave(filename;
        params = result.params,
        temperature = result.temperature,
        averages = result.averages,
        thermodynamics = result.thermodynamics,
        acceptance_rate = result.acceptance_rate
    )
    return nothing
end

"""
    load_result(filename::String) -> Dict{String, Any}

JLD2ファイルからデータを読み込み、Dictとして返す。
"""
function load_result(filename::String)
    return load(filename)
end
