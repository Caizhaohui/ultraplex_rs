# Ultraplex-RS

Ultraplex-RS 是使用Rust对开源项目 Ultraplex(https://github.com/ulelab/ultraplex) 重现与改写版本。目标是在保持核心功能与行为一致的同时，利用 Rust 的类型安全与并行生态，提供一个独立的 CLI 工具，覆盖质量修剪、UMI 抽取、5'/3' 条码匹配与分流、以及 `.fastq.gz` I/O 支持。

## 功能特性

- 质量修剪（3' 端优先）与 NextSeq 特性支持（`nextseq`）：`src/trim.rs`
- 将 UMI（条码中的 N 位）抽取并写入 read header 的 `rbc:` 字段：`src/align.rs:26`
- 5' 前缀条码匹配与（可选）3' 末端条码匹配，支持 `N` 通配与错配阈值（`threeprimemismatches`）：`src/align.rs:3`, `src/align.rs:14`
- 组合分流写出：为每个匹配条码（或条码组合）写出独立的 FASTQ 文件（支持 `.fastq.gz`）：`src/demux.rs:36`, `src/demux.rs:104`
- three_prime_only 模式：以 5' 前缀匹配为入口，在 3' 末端匹配条码并抽取 UMI，输出按样本名或组合键命名：`src/cli.rs:133`
- 条码一致性校验：3' 条码的非 N 位置对齐校验（保证 UMI 与固定位一致性）：`src/demux.rs:116`
- 并行 Reader-Workers：按批次并行处理，聚合写出（`rayon`）：`src/cli.rs:120`
- `.fastq` 与 `.fastq.gz` 输入/输出支持：`src/demux.rs:94`

> Ultraplex 的目标与行为简介：移除低质碱基、移除测序接头、将 UMI 移至 read header、检测 5'/3' 条码进行（组合）分流，并以高性能完成整 lane 的处理。

## 安装与构建

### 先决条件

- Rust（stable）与 Cargo
- Linux/macOS 环境（Windows 未经充分测试）

### 构建项目

```bash
cargo build --release
```

### 运行测试

```bash
cargo test
```

## 快速开始

### 最简用法（单端，5' 条码）

```bash
cargo run --bin ultraplex_rs -- \
  -i your.fastq.gz \
  -b barcodes.csv \
  -d out_dir \
  -o run1 \
  --gzip
```

输出文件将命名为 `ultraplex_<prefix>_<key>.fastq.gz`，无匹配的读写入 `ultraplex_<prefix>_no_match.fastq.gz`。

### three_prime_only 模式（5' + 3'，UMI 抽取）

```bash
cargo run --bin ultraplex_rs -- \
  -i reads1.fastq.gz \
  -b barcodes_5_and_3.csv \
  -d out_dir \
  -o run3p \
  --three_prime_only \
  -M 1 \
  --gzip
```

当匹配到 3' 条码时，将从条码中的 `N` 位抽取 UMI 并追加到 read header 的 `rbc:` 字段，输出按样本名或 `5bc_<5'>_3bc_<3'>` 组合命名：`src/cli.rs:166`。

### 常用参数

- `-i, --inputfastq <path>`：输入 FASTQ（支持 `.fastq.gz`）
- `-b, --barcodes <csv>`：条码 CSV（见下方格式说明）
- `-d, --directory <out>`：输出目录
- `-o, --outputprefix <prefix>`：输出前缀（用于文件名）
- `--nextseq`：启用 NextSeq 风格的质量修剪：`src/cli.rs:57`
- `--gzip`：以 `.gz` 压缩写出
- `--three_prime_only`：启用 3' 条码末端匹配 + UMI 抽取：`src/cli.rs:24`
- `-M, --threeprimemismatches <n>`：3' 匹配可允许错配数：`src/cli.rs:28`
- `-t, --threads <n>`：并行线程数：`src/cli.rs:30`
- `--keep_barcode`：保留条码本体在序列中（默认匹配到 3' 时剪去条码）：`src/cli.rs:140`
- `-l, --final_min_length <n>`：长度过滤，短于阈值的读将跳过写出：`src/cli.rs:180`
- `--ignore_no_match`：忽略无匹配的读（不写入 `no_match`）：`src/cli.rs:171`

> 目前 `-i2/--input_2` 参数已预留，但完整的成对读处理（mate adapter 移除等）仍在路线图中。

## CLI 帮助

<!-- BEGIN:CLI_HELP -->
````
Usage: ultraplex_rs [OPTIONS] --inputfastq <INPUTFASTQ>

Options:
  -i, --inputfastq <INPUTFASTQ>
          输入 FASTQ 文件路径，支持 .fastq 与 .fastq.gz；推荐 gzip 压缩
  -d, --directory <DIRECTORY>
          输出目录，仅用于写出结果文件；不批量读取目录内文件 [default: ]
  -b, --barcodes <BARCODES>
          条码 CSV。首列 5’ 条码，后续列为链接的 3’ 条码；支持 :样本名 [default: ]
  -o, --outputprefix <OUTPUTPREFIX>
          输出前缀，用于命名 ultraplex_<prefix>_<key>.fastq[.gz] [default: demux]
      --nextseq
          启用 NextSeq 风格的质量修剪（主要针对 3’ 端低质位）
      --gzip
          以 .fastq.gz 格式写出结果文件
      --three_prime_only
          启用 3’ 条码末端匹配与 UMI 抽取（结合 5’ 前缀）
  -I, --input_2 <INPUT_2>
          成对测序第二个 FASTQ 路径（预留；后续扩展） [default: ]
  -M, --threeprimemismatches <THREEPRIMEMISMATCHES>
          3’ 条码末端匹配允许的错配数 [default: 0]
  -t, --threads <THREADS>
          并行处理线程数 [default: 4]
      --keep_barcode
          匹配到 3’ 条码后是否保留条码本体在序列中
  -l, --final_min_length <FINAL_MIN_LENGTH>
          长度过滤阈值，短于该长度的读将跳过写出 [default: 0]
      --ignore_no_match
          忽略无匹配的读（不写入 no_match 文件）
  -q, --phredquality <PHREDQUALITY>
          质量修剪的 Phred 阈值（默认 30，ASCII 偏移 33） [default: 30]
  -h, --help
          Print help

````
<!-- END:CLI_HELP -->

## 条码 CSV 格式

CSV 第一列为所有 5' 条码，可选的后续列为与该 5' 条码链接的 3' 条码。条目可在尾部用冒号指定样本名，样本名应唯一。

示例：

```
NNNATGNN,
NNNCCGNN,ATG:sample2,TCA:sample3,
NNNCACNN,
```

约束：

- 5' 条码之间非 N 字符的位置必须一致（长度与具体字符可不同，N 位个数可不同）。位置相对于读的 5' 端定义。
- 3' 条码相对于读尾端（3' 端）定义，连接到同一个 5' 条码的所有 3' 条码必须具有一致的非 N 位置。
- 5' 条码与 3' 条码之间不要求一致性。

一致性检查实现：`src/demux.rs:116`。

## 行为与实现概览（代码参考）

- 质量修剪：`quality_trim_index`、`nextseq_trim_index`（`src/trim.rs`）
- 5' 前缀匹配：`prefix_match`（忽略 `N` 并支持错配）：`src/align.rs:3`
- 3' 末端匹配：`suffix_match`（忽略 `N` 并支持错配）：`src/align.rs:14`
- UMI 抽取：`extract_umi_from_suffix`（将条码中的 `N` 位对应的读端碱基收集为 UMI）：`src/align.rs:26`
- 输出写出：`write_fastq_record` 与动态 `get_writer`（根据键创建 writer）：`src/demux.rs:36`, `src/demux.rs:104`
- CLI 主流程（并行批处理）：`src/cli.rs:120`

## 与原 Python 版的差异与当前状态

- 已覆盖：质量修剪、`.fastq.gz` I/O、5'/3' 条码匹配（含 `N` 与错配）、UMI 抽取并写入 `rbc:`、按组合或样本名分流写出、并行批处理。
- 进行中：完整的适配器切除（anywhere/linked）、成对读（`-i2`）管线的适配器逻辑、更多参数对齐（如 `-m5`、TSO 参数等）。
- 文件命名与样本名映射规则已对齐原版的约定（`ultraplex_<prefix>_<key>.fastq(.gz)`），对齐逻辑见：`src/cli.rs:166` 与 `src/demux.rs:36`。

## 测试与示例数据

- 单元与端到端测试：`cargo test`
- 内置 e2e 测试覆盖：
  - 基本修剪与 5' 分流：`tests/e2e_trim.rs:28`
  - three_prime_only + UMI + 样本命名：`tests/e2e_trim.rs:55`
  - 若存在 `tests/test_simple`，自动进行若干对齐场景验证（命名映射与长度过滤）：`tests/e2e_trim.rs:94`, `tests/e2e_trim.rs:106`

## 性能与并行

- 采用 `rayon` 进行批处理并行，默认批大小 1024，`-t/--threads` 控制并行度。
- Writer 按需创建，减少文件句柄与 I/O 抖动。

## 依赖

- `clap`（CLI 参数解析）
- `needletail`（FASTQ 读入）
- `flate2`（`.gz` 压缩写出）
- `csv`（条码 CSV 解析）
- `rayon`（并行）

详见 `Cargo.toml`。

## 许可证与致谢

- 本项目遵循与原项目相同的开源精神，重现并扩展其功能以服务更广泛的 Rust 生态使用场景。
- 致谢 Ultraplex 原作者及 Cutadapt 项目对核心匹配逻辑的启发。

## 路线图

- 适配器切除与 anywhere/linked 匹配逻辑完备
- 成对读（`-i2`）适配器处理与双端一致性验证
- 更多 `tests/test_simple` 场景自动化对齐
- 统计与报告（匹配率、UMI 分布、长度分布）
