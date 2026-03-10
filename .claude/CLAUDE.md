# CLAUDE.md - FDD (Functional Declarative Design) Workflow for Julia

## Overview

このプロジェクトはFunctional Declarative Design (FDD)の方法論に従って、Juliaで2次元イジングモデルのモンテカルロシミュレーションソフトウェアを開発する。FDDはオブジェクト指向設計(OOD)に対応する関数型言語向けの体系的アプローチであり、本質的複雑さと偶有的複雑さを明確に分離することを目指す。

---

# Code Style

## Fundamental Philosophy (関数型プログラミングの基本思想)

- プログラムを「計算機への命令」ではなく「関数を適用して値を得る計算」の定義と捉える
- 関数fは値を写すものである。入力の集合(定義域)から出力の集合(値域)への写像
- どの集合からどの集合に写されるか、その集合を**型**で表現する
- 集合の表現、つまり「型と型を組み合わせた構造化」がより良くできれば、プログラムが堅牢になる

この思想に基づき、設計・実装のあらゆる場面で「関数の型シグネチャ（何を受け取り何を返すか）」を起点に考える。

## Design Principles

- **本質的複雑さ**(ドメイン固有)と**偶有的複雑さ**(実装選択由来)を分離し、後者を最小化する
- **Pure/Impure分離**: ビジネスロジックを純粋に保ち、副作用はインタープリター層に集約する
- **型による設計**: 不正な状態を表現不可能にする。ドメインモデルは型で表現する
- **eDSLによる抽象化**: 問題記述と実行詳細を分離する
- **Service Handleパターン**: テスト可能な依存性注入

## Architecture Layers

**Pure/Impure分離**:
- Pure Layer: Domain Types, eDSL定義, Business Logic, Pure Functions
- Impure Layer: Interpreters, IO Operations, State Management, External APIs

**FDDアーキテクチャレイヤー** (上から下へ):
- Application Layer (設定, ロギング, スレッド管理)
- Business Logic Layer (シナリオ, スクリプト - Pure)
- Service Layer (機能インターフェース, eDSL)
- Domain Model Layer (データ型, DSL - 完全にPure)
- Interoperability Layer (リアクティブストリーム, パイプライン)
- Persistence Layer (データ保存抽象化)

## Julia Implementation Patterns

1. **代数的データ型 (ADT)** - abstract typeによる直和型、structによる直積型
2. **eDSL設計** - コマンドを値として表現し、スクリプトを純粋なデータ構造にする
3. **型クラス風パターン（多重ディスパッチ）** - interpret/validate等のインターフェース定義と具体実装
4. **Service Handle Pattern** - 関数の集合としてのインターフェース、依存性注入、テスト用モック
5. **純粋関数を優先** - 入力のみに依存、副作用なし。副作用はインタープリター層に集約

---

# Workflow

## Development Phases

FDDは反復的トップダウンアプローチを採用する。各フェーズは前のフェーズに戻ることができる。

1. **Requirements Analysis** - ドメイン理解、マインドマップ、ユーザーシナリオ、ユースケース、Q&A
2. **Architecture Design** - レイヤー構造、サブシステム、依存関係の設計
3. **Domain Model Design** - ドメインを型とeDSLで表現
4. **Interface Design** - サブシステム間の境界定義
5. **Implementation** - 設計を実行可能コードに変換
6. **Testing** - ユニット / プロパティ / インテグレーション / E2Eテスト

## Iteration Guidelines

1. **フェーズ間の戻り**: 問題発見時は前のフェーズに戻って修正
2. **ドキュメント更新**: コード変更時は関連ドキュメントも更新
3. **テスト優先**: 実装前にテストケースを考慮
4. **レビューポイント**: 各フェーズ完了時に設計を再検討

## Git
### Commit
- 関数を1個書いたら、その都度 `git commit -m` を行う
- 細かい粒度でコミットすることで、変更の追跡と問題発生時のロールバックを容易にする

### Push
- 小さい機能: 実装完了時に `git push` する
- 大きい機能: いくつかの段階（a1, a2, ..., an）に分けて、各段階の完了時に `git push` する

適切なタイミングで`/git-gh_executor`を実行せよ.
## Notes for Claude Code

- 各フェーズは順序どおりに進行するが、前のフェーズへの反復は推奨される
- テストはモックを活用してユニット分離を確保する
- ドキュメントはコードと同期を維持する
- コードを書く際、必要になったらその都度/docs/内のファイルを確認する
    - 確認する際は, `knowledge-researcher`を用いて、そこから要約情報を受けとる
    - 外部情報を参照する必要が出た際も、`knowledge-researcher`エージェントを用いて調査をして、そのエージェントから要約情報を受け取る
- コーディングの際、プロジェクトによってはFDDの原則を全て必ずしも守らなくて良い(src/ではできる限りFDDの原則を使用)
- 実装作業フェーズのどこまで終了したかを細かいタスクが完了するごとに逐次報告する
    - 各jlファイル内のどの部分の実装が終わったかを✅で報告する
      - (例: ✅ src/Domain/Types.jl 〇〇完了)
- errorとその訂正情報は, `docs-updater`を使用して, /docs/error_log.mdに記載する
- 物理や数学の専門的な数式の計算やその実装計画には `/algo-forge` スキルを使用する

## Test
- TDDを採用.
- Unit test (小さいテスト) を全体のテストの80%にする
- 全テスト: `julia --project -e 'using Pkg; Pkg.test()'`
- 個別実行: `julia --project -e 'using Test, Ising2D; @testset "名前" begin include("test/対象ファイル.jl") end'`
