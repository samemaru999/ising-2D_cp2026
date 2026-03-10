"""
    Mocks/MockIO.jl - IOハンドラーのモック実装

ファイルシステムを使用せずにメモリ内でデータ保存をシミュレート
テストの独立性と再現性のために使用
"""

using HONetSync

# =============================================================================
# In-Memory Storage
# =============================================================================

"""
    MockStorage

メモリ内ストレージ - ファイルシステムの代替
"""
mutable struct MockStorage
    data::Dict{String, Any}
    metadata::Dict{String, Dict{String, Any}}
    write_count::Int
    read_count::Int
    last_written_key::Union{String, Nothing}
    last_read_key::Union{String, Nothing}
end

"""
    create_mock_storage() -> MockStorage

空のモックストレージを作成
"""
function create_mock_storage()
    return MockStorage(
        Dict{String, Any}(),
        Dict{String, Dict{String, Any}}(),
        0,
        0,
        nothing,
        nothing
    )
end

"""
    create_prefilled_storage() -> MockStorage

テストデータで事前に埋められたストレージを作成
"""
function create_prefilled_storage()
    storage = create_mock_storage()

    # SPSデータ
    mock_sps = create_mock_sps(N=4, Ntheta=101)
    storage.data["sps_N4_Nth101.jld2"] = mock_sps
    storage.metadata["sps_N4_Nth101.jld2"] = Dict(
        "created_at" => "2024-01-01T00:00:00",
        "type" => :sps,
        "N" => 4,
        "Ntheta" => 101
    )

    # PCFデータ
    mock_pcfs = create_mock_individual_pcfs(N=4, Ntheta=101)
    storage.data["pcf_N4_Nth101.jld2"] = mock_pcfs
    storage.metadata["pcf_N4_Nth101.jld2"] = Dict(
        "created_at" => "2024-01-01T00:00:00",
        "type" => :pcf,
        "N" => 4,
        "Ntheta" => 101
    )

    return storage
end

# =============================================================================
# Mock IO Operations
# =============================================================================

"""
    mock_save!(storage::MockStorage, key::String, data::Any; metadata=Dict())

モックストレージにデータを保存
"""
function mock_save!(storage::MockStorage, key::String, data::Any; metadata=Dict())
    storage.data[key] = deepcopy(data)
    storage.metadata[key] = merge(
        Dict("saved_at" => string(now())),
        metadata
    )
    storage.write_count += 1
    storage.last_written_key = key
    return key
end

"""
    mock_load(storage::MockStorage, key::String) -> Union{Any, Nothing}

モックストレージからデータを読み込み
"""
function mock_load(storage::MockStorage, key::String)
    storage.read_count += 1
    storage.last_read_key = key

    if haskey(storage.data, key)
        return deepcopy(storage.data[key])
    else
        return nothing
    end
end

"""
    mock_exists(storage::MockStorage, key::String) -> Bool

キーの存在確認
"""
function mock_exists(storage::MockStorage, key::String)
    return haskey(storage.data, key)
end

"""
    mock_delete!(storage::MockStorage, key::String) -> Bool

データの削除
"""
function mock_delete!(storage::MockStorage, key::String)
    if haskey(storage.data, key)
        delete!(storage.data, key)
        delete!(storage.metadata, key)
        return true
    end
    return false
end

"""
    mock_list_keys(storage::MockStorage) -> Vector{String}

全キーのリスト
"""
function mock_list_keys(storage::MockStorage)
    return collect(keys(storage.data))
end

"""
    mock_get_metadata(storage::MockStorage, key::String) -> Union{Dict, Nothing}

メタデータの取得
"""
function mock_get_metadata(storage::MockStorage, key::String)
    if haskey(storage.metadata, key)
        return storage.metadata[key]
    end
    return nothing
end

# =============================================================================
# Mock IO Service (Dependency Injection Pattern)
# =============================================================================

"""
    MockIOService

モックIOサービス構造体（依存性注入用）
"""
struct MockIOService
    storage::MockStorage
    save::Function
    load::Function
    exists::Function
    delete::Function
    list_keys::Function
end

"""
    create_mock_io_service(; prefilled=false) -> MockIOService

モックIOサービスを作成
"""
function create_mock_io_service(; prefilled=false)
    storage = prefilled ? create_prefilled_storage() : create_mock_storage()

    return MockIOService(
        storage,
        (key, data; metadata=Dict()) -> mock_save!(storage, key, data; metadata=metadata),
        (key) -> mock_load(storage, key),
        (key) -> mock_exists(storage, key),
        (key) -> mock_delete!(storage, key),
        () -> mock_list_keys(storage)
    )
end

# =============================================================================
# Failure Injection for IO Testing
# =============================================================================

"""
    create_failing_io_service(fail_on::Symbol) -> MockIOService

指定した操作で失敗するIOサービスを作成

# Arguments
- `fail_on::Symbol`: `:save`, `:load`, `:delete` のいずれか
"""
function create_failing_io_service(fail_on::Symbol)
    storage = create_mock_storage()

    function failing_save(key, data; metadata=Dict())
        if fail_on == :save
            error("Simulated IO error on save")
        end
        return mock_save!(storage, key, data; metadata=metadata)
    end

    function failing_load(key)
        if fail_on == :load
            error("Simulated IO error on load")
        end
        return mock_load(storage, key)
    end

    function failing_delete(key)
        if fail_on == :delete
            error("Simulated IO error on delete")
        end
        return mock_delete!(storage, key)
    end

    return MockIOService(
        storage,
        failing_save,
        failing_load,
        (key) -> mock_exists(storage, key),
        failing_delete,
        () -> mock_list_keys(storage)
    )
end

"""
    create_slow_io_service(delay_seconds::Float64) -> MockIOService

意図的に遅いIOサービスを作成（タイムアウトテスト用）
"""
function create_slow_io_service(delay_seconds::Float64)
    storage = create_mock_storage()

    function slow_save(key, data; metadata=Dict())
        sleep(delay_seconds)
        return mock_save!(storage, key, data; metadata=metadata)
    end

    function slow_load(key)
        sleep(delay_seconds)
        return mock_load(storage, key)
    end

    return MockIOService(
        storage,
        slow_save,
        slow_load,
        (key) -> mock_exists(storage, key),
        (key) -> mock_delete!(storage, key),
        () -> mock_list_keys(storage)
    )
end

# =============================================================================
# Statistics and Verification
# =============================================================================

"""
    get_io_stats(storage::MockStorage) -> NamedTuple

IOの統計情報を取得
"""
function get_io_stats(storage::MockStorage)
    return (
        total_items = length(storage.data),
        write_count = storage.write_count,
        read_count = storage.read_count,
        last_written = storage.last_written_key,
        last_read = storage.last_read_key
    )
end

"""
    reset_io_stats!(storage::MockStorage)

IO統計をリセット
"""
function reset_io_stats!(storage::MockStorage)
    storage.write_count = 0
    storage.read_count = 0
    storage.last_written_key = nothing
    storage.last_read_key = nothing
end

"""
    verify_io_sequence(storage::MockStorage, expected_writes::Int, expected_reads::Int) -> Bool

期待されるIO操作回数を検証
"""
function verify_io_sequence(storage::MockStorage, expected_writes::Int, expected_reads::Int)
    return storage.write_count == expected_writes && storage.read_count == expected_reads
end

