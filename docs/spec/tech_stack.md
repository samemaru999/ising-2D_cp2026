# 技術スタック仕様書

**プロジェクト名**: 2次元イジングモデル モンテカルロシミュレーション (Ising2D)
**バージョン**: 1.0
**作成日**: 2026-03-10

---

## 1. 開発言語

### 1.1 Julia

| 項目 | 仕様 |
|------|------|
| **言語** | Julia |
| **必須バージョン** | 1.10 以上 |
| **推奨バージョン** | 1.11.x (最新安定版) |

### 1.2 Julia選定理由

| 理由 | FDDとの関連 |
|------|------------|
| 高速な数値計算（JITコンパイル） | 本質的複雑さへの集中 |
| 科学計算エコシステム | 偶有的複雑さの削減 |
| 多重ディスパッチによる拡張性 | eDSLパターンの実現 |
| REPL/Jupyterでの対話的開発 | 反復的開発の支援 |
| 型システム（パラメトリック型） | 型による設計原則 |
| 効率的な配列操作 | 2次元スピン格子の操作 |

---

## 2. 主要パッケージ

### 2.1 必須パッケージ

| パッケージ | 用途 | FDDレイヤー |
|------------|------|------------|
| Random (stdlib) | モンテカルロ乱数生成 | Interpreters (Impure) |
| Statistics (stdlib) | 統計量計算（平均、分散） | Logic (Pure) |
| LinearAlgebra (stdlib) | 線形代数演算 | Domain (Pure) |
| JLD2.jl | データ保存・読み込み | Runtime (Impure) |
| GLMakie.jl | 可視化（ヒートマップ、プロット） | App (Impure) |
| ProgressMeter.jl | 進捗表示 | App (Impure) |

### 2.2 補助パッケージ

| パッケージ | 用途 | 重要度 |
|------------|------|--------|
| StaticArrays.jl | 高速な固定長配列（必要に応じて） | 任意 |
| LsqFit.jl | 最小二乗フィッティング（臨界指数推定） | 任意 |
| Test (stdlib) | ユニットテスト | 推奨 |
| BenchmarkTools.jl | 性能測定 | 任意 |

### 2.3 モンテカルロ法に関する注意事項

```
重要: 乱数シード管理について
----------------------------
再現可能なシミュレーション結果のため、乱数シードを必ず設定可能にすること。
乱数生成はImpure処理としてInterpreters層に限定する。

例:
  rng = MersenneTwister(seed)
  rand(rng) < exp(-ΔE / T)
```

---

## 3. 開発環境

### 3.1 推奨環境

```
OS: macOS / Linux / Windows (WSL2推奨)
CPU: 2コア以上
メモリ: 8GB以上
ストレージ: SSD推奨
```

### 3.2 開発ツール

| ツール | 用途 | 推奨 |
|--------|------|------|
| エディタ/IDE | コード編集 | VS Code + Julia拡張 |
| Jupyter | 対話的開発 | IJulia |
| REPL | 動作確認 | Julia標準REPL |
| バージョン管理 | ソースコード管理 | Git |

---

## 4. 依存関係管理

### 4.1 Project.toml

```toml
name = "Ising2D"
uuid = "..."
authors = ["S. Kuroiwa"]
version = "1.0.0"

[deps]
GLMakie = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
julia = "1.10"
JLD2 = "0.4"
GLMakie = "0.10"
ProgressMeter = "1"
```

---

## 5. FDDレイヤーとパッケージの対応

```
+-------------------------------------------------------------+
|                    Application Layer (Impure)                |
|  GLMakie.jl, ProgressMeter.jl                                |
|  - 可視化、進捗表示                                           |
+-------------------------------------------------------------+
|                    Business Logic Layer (Pure)               |
|  Statistics (stdlib)                                         |
|  - 統計量計算                                                 |
+-------------------------------------------------------------+
|                    Domain Model Layer (Pure)                 |
|  LinearAlgebra (stdlib)                                      |
|  - ハミルトニアン計算、物理量計算                              |
+-------------------------------------------------------------+
|                    Interoperability Layer (Impure)           |
|  Random (stdlib)                                              |
|  - メトロポリス法の乱数生成                                    |
+-------------------------------------------------------------+
|                    Persistence Layer (Impure)                |
|  JLD2.jl                                                      |
|  - データの保存・読み込み                                      |
+-------------------------------------------------------------+
```

### Pure/Impureパッケージ分類

| 分類 | パッケージ | 理由 |
|------|-----------|------|
| **Pure** | LinearAlgebra | 数学演算、副作用なし |
| **Pure** | Statistics | 統計計算、副作用なし |
| **Impure** | Random | 乱数生成（内部状態変更） |
| **Impure** | JLD2.jl | ファイルIO |
| **Impure** | GLMakie.jl | 描画（副作用） |
| **Impure** | ProgressMeter.jl | 標準出力への書き込み |

---

## 付録

### A. 関連ドキュメント

- [システムアーキテクチャ設計書](./system_architecture.md)
- [データ構造定義書](./data_structure.md)
- [ディレクトリ構造](./directory_structure.md)
- [要件定義書](./requirement.md)
