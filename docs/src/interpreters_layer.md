# Runtime/MonteCarlo Layer Documentation

`src/Runtime/MonteCarlo.jl` は、メトロポリス法によるモンテカルロ更新と温度掃引シミュレーションを担当する。乱数を使用するため **Impure** レイヤー。

## ファイル

- `src/Runtime/MonteCarlo.jl`: メトロポリス更新、スイープ、温度掃引

---

## 関数一覧

### 低レベル関数

#### `metropolis_step!(lattice, beta; J=1.0) -> Bool`

ランダムに1サイトを選び、メトロポリス法でスピン反転を試みる。

1. サイト $(i,j)$ をランダムに選択
2. $\Delta E = $ `delta_energy(lattice, i, j; J)` を計算
3. 受容判定: $\Delta E \leq 0$ または $\text{rand}() < \exp(-\beta \Delta E)$
4. 受容時にスピン反転: `lattice[i,j] *= -1`

- **戻り値**: 反転が受容されたら `true`
- **ベンチマーク**: ~16.6 ns, 0 allocations

#### `sweep!(lattice, beta; J=1.0) -> Int`

$L \times L$ 回の `metropolis_step!` を実行（1モンテカルロスイープ）。

- **戻り値**: 受容された回数
- **ベンチマーク** (L=32): ~15.8 μs, 0 allocations
- **1スピン更新あたり**: ~15.4 ns

### 高レベル関数

#### `run_single_temperature(lattice, beta, N, config; J=1.0) -> (ThermalAverages, Int)`

1つの温度で平衡化 + 測定を実行。

1. **平衡化**: `config.n_equilibration` スイープ（結果は捨てる）
2. **測定ループ**: `config.n_sweeps` スイープ
   - `config.sample_interval` ごとに物理量を測定
   - `energy(lattice; J) / N`, `magnetization(lattice)` を計算
   - `ThermalAverages` に累積

- **戻り値**: `(ThermalAverages, total_accepted)`
- `lattice` は in-place で更新される

#### `run_temperature_sweep(params, config, temperatures) -> TemperatureSweepResult`

温度配列をループし、各温度でシミュレーションを実行。

1. 温度を**高温側からソート**（降順）
2. `random_lattice(L)` で初期化
3. 各温度で `run_single_temperature` → `compute_thermodynamics`
4. 前の温度の格子状態を引き継ぐ（高温から徐冷）
5. 結果を温度昇順に並べ替えて返す

---

## メトロポリス法の物理

受容確率:
$$P(\text{accept}) = \min\left(1, \exp(-\beta \Delta E)\right)$$

2次元正方格子では $\Delta E \in \{-8J, -4J, 0, +4J, +8J\}$ の5値のみ。これを利用した事前テーブル化 (`BoltzmannTable`) は将来の最適化候補。

## パフォーマンス特性

| 関数 | L=32 時間 | アロケーション |
|------|-----------|-------------|
| `metropolis_step!` | ~17 ns | 0 |
| `sweep!` | ~16 μs | 0 |
| 1スピン更新 | ~15 ns | 0 |
