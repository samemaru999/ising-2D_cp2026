# Runtime/IO Layer Documentation

`src/Runtime/IO/` ディレクトリは、計算結果の永続化、ファイル入出力、シリアライゼーションを担当します。これは Impure なレイヤーであり、ファイルシステムへの副作用を持ちます。

## ファイル一覧

- `IOHandler.jl`: JLD2形式でのデータ保存・読み込みとファイル管理。

---

## 1. IOHandler.jl

`JLD2.jl` パッケージを使用して、Juliaのデータ構造をバイナリ形式で効率的に保存・読み込みします。

### ファイル命名 (`generate_filename`)
以下の形式で一貫したファイル名を生成します。

`{Type}_N{N}_Nth{Ntheta}_{Option}_{Date}.jld2`

- 例: `sps_N4_Nth101_20240101.jld2`
- オプション: `_linear`, `_rotation`, `_0.10` など

### コマンド処理
- **`interpret(cmd::SaveData)`**: 汎用的なデータ保存を実行します。フォーマットとして `:jld2` または `:csv` をサポートします。
- **`interpret(cmd::LoadData)`**: ファイルからデータを読み込み、指定された型として返します。

### JLD2 保存関数
各データ型ごとに特化した保存関数が用意されています。これらはメタデータ（作成日時、バージョン情報）を自動的に付与します。

- **`save_sps`**: 安定周期解 (`StablePeriodicSolution`) を保存。
- **`save_psf`**: 位相感受関数 (`PhaseSensitivityFunction`) を保存。
- **`save_total_pcf`**: 全体位相結合関数 (`TotalPCF`) を保存。
- **`save_optimization_result`**: 最適化結果 (`OptimizationResult`) を保存。
- **`save_pipeline_results`**: パイプライン実行結果全体をディレクトリに一括保存します。

### JLD2 読み込み関数
保存されたファイルからデータを復元します。

- **`load_sps`**: `StablePeriodicSolution` を復元（`NetworkParams` のデシリアライズを含む）。
- **`load_psf`**: `PhaseSensitivityFunction` を復元。
- **`load_total_pcf`**: `TotalPCF` を復元。
- **`load_optimization_result`**: `OptimizationResult` を復元。

### CSV出力
位相ダイナミクスなどの時系列データは、解析のためにCSV形式でも出力可能です。

- **`save_phase_time_series_csv`**: 位相時系列データ (`PhaseTimeSeries`) をCSVとして保存します。

### ユーティリティ
- **`list_data_files`**: 指定ディレクトリ内の `.jld2` ファイル一覧を取得。
- **`get_file_info`**: ファイルのヘッダー情報を読み取り、メタデータを返します。
