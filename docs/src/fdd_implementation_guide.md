# FDD (Functional Declarative Design) と Ising2D の実装

FDDの思想が Ising2D プロジェクトにどのように適用されているかを解説する。

---

## 1. FDDとは何か

### 基本理念

FDD (Functional Declarative Design) は関数型プログラミングの原則に基づく設計手法:

1. **本質的複雑さと偶有的複雑さの分離**: ドメイン固有の物理を純粋に保ち、実装都合の複雑さを最小化
2. **Pure/Impure分離**: 副作用のないコードと副作用を持つコードの明確な分離
3. **型による設計**: 不正な状態を型で表現不可能にする

### Ising2Dでの適用

物理シミュレーションの本質的複雑さ（イジングモデルの物理、統計力学）と偶有的複雑さ（乱数生成、ファイルIO、プロット）を分離する。

---

## 2. Pure/Impure分離

### Pure Layer（副作用なし）

```
src/Types.jl              # 型定義
src/Ising/Lattice.jl      # 格子生成
src/Ising/Physics.jl      # エネルギー・磁化・ΔE
src/Ising/Statistics.jl   # 熱力学量の導出
```

特徴:
- 同じ入力 → 常に同じ出力
- 乱数を使わない
- ファイルIOなし
- 単独でテスト可能

### Impure Layer（副作用あり）

```
src/Runtime/MonteCarlo.jl    # メトロポリス更新（乱数使用）
src/Runtime/IO.jl            # JLD2ファイル保存・読み込み
src/Runtime/Visualization.jl # GLMakieプロット
```

特徴:
- 乱数生成、ファイルシステム、GUIに依存
- Pure層の関数を呼び出す

---

## 3. 型による設計

### 不正な状態の排除

内部コンストラクタでバリデーションを行い、不正な値を型レベルで排除:

```julia
struct IsingParams
    L::Int
    J::Float64
    function IsingParams(L::Int, J::Float64)
        L >= 2 || throw(ArgumentError("L must be ≥ 2"))
        J != 0.0 || throw(ArgumentError("J must be non-zero"))
        new(L, J)
    end
end
```

これにより `IsingParams` が存在すれば、`L >= 2` かつ `J != 0` が保証される。

### 結果型の階層

```
IsingParams + SimulationConfig    # 入力パラメータ
       ↓
ThermalAverages                    # MC測定の累積（mutable）
       ↓  compute_thermodynamics (Pure)
ThermodynamicQuantities            # 導出された物理量（immutable）
       ↓
SingleTemperatureResult            # 単一温度の完全な結果
       ↓
TemperatureSweepResult             # 温度掃引の完全な結果
```

データは下流に流れ、各段階で型が保証する不変条件が増える。

---

## 4. 具体例: 物理計算の Pure/Impure 分離

### Pure: `delta_energy`（Physics.jl）

```julia
function delta_energy(lattice::Matrix{Int}, i::Int, j::Int; J::Float64=1.0)
    L = size(lattice, 1)
    s = lattice[i, j]
    neighbors = lattice[mod1(i-1, L), j] + lattice[mod1(i+1, L), j] +
                lattice[i, mod1(j-1, L)] + lattice[i, mod1(j+1, L)]
    return 2.0 * J * s * neighbors
end
```

- 入力（格子、座標）のみに依存、常に同じ結果
- テスト容易: 既知の格子配置で期待値と比較するだけ

### Impure: `metropolis_step!`（MonteCarlo.jl）

```julia
function metropolis_step!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0)
    L = size(lattice, 1)
    i = rand(1:L)           # Impure: 乱数
    j = rand(1:L)           # Impure: 乱数
    dE = delta_energy(lattice, i, j; J=J)  # Pure関数の呼び出し
    if dE <= 0 || rand() < exp(-beta * dE)  # Impure: 乱数
        lattice[i, j] *= -1  # Impure: 状態変更
        return true
    end
    return false
end
```

- `rand()` と `lattice` の変更が副作用
- Pure な `delta_energy` を内部で呼び出す

### Pure: `compute_thermodynamics`（Statistics.jl）

```julia
function compute_thermodynamics(avg::ThermalAverages, N::Int, T::Float64)
    # 累積値 → 熱力学量への純粋な変換
    mean_e = avg.sum_e / avg.n_samples
    C = N / T^2 * (mean_e2 - mean_e^2)
    # ...
    return ThermodynamicQuantities(mean_e, mean_abs_m, C, chi, U_L)
end
```

- `ThermalAverages`（Impureな測定ループで蓄積）を受け取り、Pure な計算で `ThermodynamicQuantities` を返す
- Impure → Pure の境界

---

## 5. テスタビリティ

Pure/Impure 分離の最大の利点はテストの容易さ:

```julia
# Pure関数は単独でテスト可能（モック不要）
@testset "delta_energy" begin
    lattice = [1 1; 1 1]  # 全スピン+1
    @test delta_energy(lattice, 1, 1; J=1.0) == 8.0
end

@testset "compute_thermodynamics" begin
    avg = ThermalAverages(...)  # 既知の値を設定
    thermo = compute_thermodynamics(avg, 16, 2.0)
    @test thermo.specific_heat ≈ expected_C
end
```

テスト構成（TDD採用）:
- **Unit tests (80%)**: Pure関数中心。高速・決定的
- **Integration / E2E (20%)**: Impure層を含むシナリオテスト

---

## 6. ディレクトリ構造

```
src/
├── Ising2D.jl              # モジュール定義 + exports
├── Types.jl                # Domain Layer [Pure] - 全型定義
├── Ising/                  # Logic Layer [Pure]
│   ├── Lattice.jl          # 格子生成
│   ├── Physics.jl          # 物理量計算
│   └── Statistics.jl       # 熱力学量の導出
└── Runtime/                # Runtime Layer [Impure]
    ├── MonteCarlo.jl       # メトロポリスMC、温度掃引
    ├── IO.jl               # JLD2保存・読み込み
    └── Visualization.jl    # GLMakieプロット

scripts/                    # 実行スクリプト
├── quick_start.jl          # 単一温度MC
├── temperature_sweep.jl    # 温度掃引（固定パラメータ）
├── run_simulation.jl       # 温度掃引（CLI引数対応）
├── benchmark.jl            # パフォーマンス計測
└── plot_results.jl         # 結果の再プロット

test/
├── runtests.jl
├── Unit/                   # 単体テスト (80%)
│   ├── test_types.jl
│   ├── test_lattice.jl
│   ├── test_physics.jl
│   ├── test_statistics.jl
│   ├── test_observables.jl
│   ├── test_results.jl
│   └── test_montecarlo.jl
├── Ising/
│   └── test_integration.jl
└── E2E/
    └── test_simulation.jl
```
