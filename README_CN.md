# nf-core/sarek 中文指南

## 目录
1. [什么是 sarek](#什么是-sarek)
2. [WGS 分析流程总览](#wgs-分析流程总览)
3. [代码结构导读](#代码结构导读)
4. [关键文件速查表](#关键文件速查表)
5. [核心概念解释](#核心概念解释)
6. [生物学背景对照表](#生物学背景对照表)
7. [常用参数说明](#常用参数说明)
8. [三种变异检测模式](#三种变异检测模式)

---

## 什么是 sarek

**一句话介绍**: nf-core/sarek 是一个完整的全基因组测序 (WGS) / 全外显子组测序 (WES) 变异检测流程工具，用 Nextflow 编写，可以从原始测序数据 (FASTQ) 一步步处理到最终的变异检测结果 (VCF) 和质量报告。

**核心特点**:
- 📊 **端到端自动化**: 从 FASTQ 到注释的 VCF，一条命令搞定
- 🧬 **支持多种变异检测工具**: HaplotypeCaller、DeepVariant、Manta 等
- 🚀 **高效并行处理**: 将基因组分块处理，大幅加速分析
- 📈 **生产级质量**: nf-core 社区维护，经过严格测试
- 💻 **可移植性强**: 支持本地运行、集群、云平台 (AWS/GCP) 等

---

## WGS 分析流程总览

### 完整流程图

```
┌─────────────────────────────────────────────────────────────────────┐
│                       FASTQ 原始测序数据                              │
│              (来自测序仪，包含碱基质量信息)                          │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                         ▼
            ┌────────────────────────────┐
            │  1️⃣  质控 (FastQC)          │
            │  评估测序质量和偏差          │
            │  Output: HTML 质控报告      │
            └────────────┬────────────────┘
                         │
                         ▼
            ┌────────────────────────────┐
            │  2️⃣  去接头 (FASTP)         │
            │  移除测序接头和低质量 reads │
            │  Output: 清洁 FASTQ        │
            └────────────┬────────────────┘
                         │
                         ▼
            ┌────────────────────────────┐
            │  3️⃣  比对 (BWA-MEM/Bowtie2)│
            │  将 reads 比对到参考基因组  │
            │  Output: SAM 文件 → BAM    │
            └────────────┬────────────────┘
                         │
                         ▼
            ┌────────────────────────────┐
            │  4️⃣  标记重复 (MarkDup)     │
            │  找出测序过程产生的重复 read│
            │  Output: 标记重复的 BAM    │
            └────────────┬────────────────┘
                         │
                         ▼
            ┌────────────────────────────┐
            │  5️⃣  碱基质量校正 (BQSR)   │
            │  校正测序机产生的系统偏差   │
            │  Output: 校正后 BAM        │
            └────────────┬────────────────┘
                         │
                         ▼
         ┌──────────────────────────────────────┐
         │    6️⃣  变异检测 (多工具可选)          │
         ├──────────────────────────────────────┤
         │  A) HaplotypeCaller (推荐 SNP/InDel)│
         │     局部重组装 + 贝叶斯推断          │
         │                                      │
         │  B) DeepVariant (深度学习，高准确)  │
         │     CNN 从 reads 堆叠图像推断        │
         │                                      │
         │  C) Manta (结构变异)                │
         │     大型重排、倒位等                │
         │                                      │
         │  Output: gVCF 和 VCF 文件            │
         └──────────────┬───────────────────────┘
                        │
                        ▼
         ┌──────────────────────────────────────┐
         │  7️⃣  联合分型 (GenotypeGVCFs)       │
         │  多样本一起重新分型，提高准确度      │
         │  (可选，仅多样本时)                  │
         │  Output: 整合的 VCF                 │
         └──────────────┬───────────────────────┘
                        │
                        ▼
         ┌──────────────────────────────────────┐
         │  8️⃣  变异注释 (VEP/SnpEff)          │
         │  标注变异的生物学意义                 │
         │  - 编码区/非编码区                   │
         │  - 已知数据库 (gnomAD, ClinVar)     │
         │  - 预测影响 (missense, frameshift)   │
         │  Output: 注释的 VCF                 │
         └──────────────┬───────────────────────┘
                        │
                        ▼
         ┌──────────────────────────────────────┐
         │  9️⃣  质控报告 (MultiQC)             │
         │  汇总所有分析步骤的统计               │
         │  - 比对率、覆盖深度                   │
         │  - 变异数量、分布                     │
         │  Output: HTML 交互式报告             │
         └──────────────────────────────────────┘
                        │
                        ▼
            ┌────────────────────────────┐
            │  ✅ 最终输出: 注释的 VCF   │
            │      + HTML 质控报告        │
            │      + 详细日志              │
            └────────────────────────────┘
```

### 流程说明

| 步骤 | 工具 | 输入 | 输出 | 生物学意义 |
|------|------|------|------|----------|
| 1️⃣ 质控 | FastQC | FASTQ | HTML 报告 | 评估测序是否成功，数据质量 |
| 2️⃣ 去接头 | FASTP | FASTQ | 清洁 FASTQ | 移除人工引入的接头 |
| 3️⃣ 比对 | BWA-MEM | FASTQ + 参考 | BAM | reads 在基因组中的位置 |
| 4️⃣ 标记重复 | MarkDuplicates | BAM | BAM | 移除 PCR 重复 |
| 5️⃣ 碱基校正 | BQSR | BAM | BAM | 校正测序机系统误差 |
| 6️⃣ 变异检测 | HC/DV/Manta | BAM | VCF | 找出 SNP/InDel/SV |
| 7️⃣ 联合分型 | GenotypeGVCFs | gVCF | VCF | 多样本共同分型 |
| 8️⃣ 注释 | VEP | VCF | VCF | 变异功能预测 |
| 9️⃣ 报告 | MultiQC | 所有日志 | HTML | 质量总结 |

---

## 代码结构导读

### 目录树概览

```
sarek/
├── README.md                          # 官方英文文档
├── nextflow.config                    # 主配置文件
├── conf/                              # 配置目录
│   ├── modules.config                 # 每个进程的参数配置
│   ├── genomes.config                 # 参考基因组配置 (hg38, hg37 等)
│   └── test.config                    # 测试数据配置
│
├── workflows/                         # 主工作流目录
│   ├── sarek.nf                       # 主工作流 (从这里开始)
│   └── *.nf                           # 不同分析模式的工作流
│
├── subworkflows/                      # 子流程目录 (可复用的步骤)
│   ├── local/
│   │   ├── fastq_align/               # 🔹 reads 比对
│   │   ├── bam_qc/                    # 🔹 BAM 质控
│   │   ├── bam_variant_calling_*      # 🔹 变异检测 (多种工具)
│   │   └── *.nf                       # 其他小流程
│   └── nf-core/
│       └── ...                        # 来自 nf-core 社区的子流程
│
├── modules/                           # 模块目录 (执行单个工具)
│   ├── nf-core/                       # nf-core 社区的模块
│   │   ├── gatk4/                     # GATK4 工具
│   │   │   ├── haplotypecaller/       # 🔹 HaplotypeCaller 模块
│   │   │   │   └── main.nf           # 实际执行代码
│   │   │   ├── mergevcfs/             # VCF 合并
│   │   │   └── ...
│   │   ├── bwa/                       # BWA 比对工具
│   │   ├── samtools/                  # SAMtools (BAM 处理)
│   │   ├── deepvariant/               # DeepVariant
│   │   └── ...
│   └── local/
│       └── ...                        # 自定义本地模块
│
└── assets/
    ├── schema_input.json              # 输入参数规范
    └── ...
```

### 看代码的顺序

**新手推荐路线**:

1. **第一步**: 打开 `workflows/sarek.nf`
   - 这是整个流程的"导演脚本"
   - 看里面的 `include` 语句，了解用了哪些子流程
   - 看 `workflow` 块，了解执行顺序

2. **第二步**: 理解 4 个关键子流程
   - `subworkflows/local/fastq_align/main.nf` — reads 比对
   - `subworkflows/local/bam_qc/main.nf` — BAM 质控
   - `subworkflows/local/bam_variant_calling_haplotypecaller/main.nf` — 变异检测
   - `subworkflows/local/bam_variant_calling_deepvariant/main.nf` — DeepVariant 变异检测

3. **第三步**: 深入感兴趣的模块
   - 如 `modules/nf-core/gatk4/haplotypecaller/main.nf` — 实际执行 HaplotypeCaller
   - 这里才是真正的命令行代码

4. **第四步**: 参考配置
   - `conf/modules.config` — 看每个工具的参数设定
   - `conf/genomes.config` — 看参考基因组和索引位置

**关键概念**:
- **include**: 导入其他 Nextflow 文件中的工作流或模块
- **workflow**: 定义一个工作流程，包含多个步骤
- **process**: 执行单个工具的最小单位
- **channel**: 进程间传递数据的"管道"

---

## 关键文件速查表

### 最重要的 5 个文件

| 文件路径 | 作用 | 何时查看 |
|---------|------|---------|
| `workflows/sarek.nf` | 整个流程的大脑 | 想了解整体流程 |
| `subworkflows/local/fastq_align/main.nf` | reads 如何比对 | 不清楚比对步骤 |
| `subworkflows/local/bam_variant_calling_haplotypecaller/main.nf` | 如何检测变异 | 想理解变异检测原理 |
| `modules/nf-core/gatk4/haplotypecaller/main.nf` | HaplotypeCaller 实际命令 | 需要调参或排故 |
| `conf/modules.config` | 每个工具的参数 | 想改某工具的参数 |

### 其他常看文件

| 文件 | 用途 |
|------|------|
| `conf/genomes.config` | 参考基因组位置和索引文件 |
| `nextflow.config` | 全局配置 (执行方式、内存等) |
| `conf/test.config` | 用小数据测试流程 |
| `assets/schema_input.json` | 输入参数规范文档 |

---

## 核心概念解释

### 1️⃣ Channel (数据流)

**类比**: 想象一条流水线上的传送带，samples 就是流水线上的工件。

```
Channel 的三种形式:
- emit:   输出数据的地方
- take:   接收输入数据的地方
- main:   处理数据的逻辑
```

**例子**:
```
workflow FASTQ_ALIGN {
    take:
    reads                    // 输入: FASTQ 文件流
    index                    // 输入: 参考基因组索引

    main:
    BWAMEM1_MEM(reads, index, sort)  // 处理数据: 比对

    emit:
    bam = BWAMEM1_MEM.out.bam        // 输出: BAM 文件流
}
```

**关键点**:
- Channel 是**异步**传递数据的 — 不是等一个完成后才开始下一个
- Channel 支持**并行处理** — 1000 个样本可以同时进行

### 2️⃣ Process (执行步骤)

**类比**: 一个独立的任务工作室，每次可以处理一个样本或多个样本。

```
process 结构:
┌─────────────────────────────────────┐
│ process 进程名 {                    │
│   tag "$meta.id"         # 标签    │
│   label 'process_*'      # 资源    │
│   container "..."        # 环境    │
│   input:  ...            # 输入    │
│   output: ...            # 输出    │
│   script: { ... }        # 命令    │
│ }                         │
└─────────────────────────────────────┘
```

**例子** — 运行 HaplotypeCaller:
```
process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"                    # 标签: 样本 ID
    label 'process_low'               # 资源: 低 CPU 低内存
    container 'gatk4:...'             # 环境: 包含 GATK4 的容器

    input:
    tuple val(meta), path(bam)        # 输入: 样本 ID + BAM 文件

    output:
    tuple val(meta), path("*.vcf.gz") # 输出: 样本 ID + VCF 文件

    script:
    """
    gatk HaplotypeCaller --input $bam ...
    """
}
```

**关键点**:
- 每个 process 在**隔离环境**中运行 (容器或 conda)
- 可以设置 **CPU/内存/时间**限制
- Process 是**可重用**的，不同工作流可以共用

### 3️⃣ Subworkflow (子流程)

**类比**: 一个可复用的"小工厂"，包含多个 process 组成的工作步骤。

```
Subworkflow 结构:
┌────────────────────────────────────┐
│ workflow 子流程名 {                │
│   take:     ...                    │
│   main:                            │
│     process1 | process2 | process3 │
│   emit:     ...                    │
│ }                                  │
└────────────────────────────────────┘
```

**例子** — 比对子流程:
```
workflow FASTQ_ALIGN {
    take:
    reads, index, sort

    main:
    // 选择一个比对工具运行
    BWAMEM1_MEM(reads, index, sort)

    emit:
    bam = BWAMEM1_MEM.out.bam
}
```

**关键点**:
- Subworkflow = N 个 process 的组合
- 可以嵌套调用其他 subworkflow
- 使代码模块化、可重用

### 4️⃣ Include (导入工具)

**类比**: 在 Python 中 `import numpy as np`，引入一个库或函数。

```
include { 工作流或模块名 } from '文件路径'
```

**例子**:
```
include { BWAMEM1_MEM }        from '../modules/bwa/mem/main'
include { GATK4_HAPLOTYPECALLER } from '../modules/gatk4/haplotypecaller/main'

workflow {
    BWAMEM1_MEM(...)           // 使用引入的模块
}
```

**关键点**:
- `include` 把外部的 process/workflow 导入当前文件
- 支持别名: `as BWAMEM1_MEM` 重命名
- 使代码分离，便于维护

### 5️⃣ Meta Map (样本元信息)

**类比**: 每个样本附带的"身份证"，记录样本的所有属性。

```
Meta map 例子:
[
    id: 'sample1',                    # 样本 ID
    single_end: false,                # 是否为单端测序
    read_group: 'rg1',                # 测序库 ID
    interval_name: 'chr1',            # 当前处理的基因组区间
    num_intervals: 100,               # 总区间数
    variantcaller: 'haplotypecaller'  # 使用的变异检测工具
]
```

**为什么重要**:
- 数据和元信息总是"成对"传递
- 帮助 process 追踪是哪个样本的数据
- 用于后续的分组、合并操作

**在代码中的用法**:
```
tuple val(meta), path(bam)    // meta 是字典，bam 是文件路径

// 在 process 中使用 meta:
tag "$meta.id"                // 用样本 ID 标记 log 输出
output: tuple val(meta), ...  // 传递 meta 到下一步
```

---

## 生物学背景对照表

### 每个分析步骤的生物学意义

| 步骤 | 生物学问题 | 工具 | 原理 | 输出的意义 |
|------|----------|------|------|----------|
| **FastQC** | 测序质量如何? | FastQC | 统计碱基质量分布 | 是否需要过滤低质 reads |
| **FASTP** | 哪些 reads 质量太差? | FASTP | 剪除低质段和接头 | 提升下游分析准确度 |
| **BWA-MEM** | reads 来自基因组哪里? | BWA-MEM | seed-extend 算法 | reads 在基因组上的位置 |
| **MarkDup** | 这些 reads 是 PCR 重复吗? | MarkDuplicates | 按序列相同性判断 | 避免 PCR 偏差的遗传影响 |
| **BQSR** | 碱基质量值准确吗? | GATK BQSR | 用已知变异校正偏差 | 校正测序机系统误差 |
| **HaplotypeCaller** | 这里有变异吗? 是什么类型? | GATK HC | 局部重组装 + 贝叶斯推断 | SNP/InDel 位置和基因型 |
| **DeepVariant** | 用机器学习看有什么变异? | DeepVariant | CNN 从图像推断 | 更准确的 SNP/InDel 检测 |
| **Manta** | 有大的基因组重排吗? | Manta | 读对异常和跳跃分析 | SV (>50bp) 的位置 |
| **GenotypeGVCFs** | 所有样本在这里的基因型是什么? | GATK GG | 贝叶斯多样本分型 | 样本间一致的基因型 |
| **VEP/SnpEff** | 这个变异有什么影响? | VEP/SnpEff | 对标注数据库 | 变异的功能预测 |

### 几个关键的生物学概念

#### SNP vs InDel
- **SNP** (Single Nucleotide Polymorphism): 单个碱基改变 — 最常见的遗传变异
  ```
  ref:  A T G C C G
  alt:  A T A C C G  ← C→A 替换 = SNP
  ```

- **InDel** (Insertion/Deletion): 插入或删除一段序列
  ```
  ref:  A T G C C G
  alt:  A T G - - G  ← 删除 CC = deletion
  alt:  A T G C T C C G  ← 插入 TC = insertion
  ```

#### 基因型 (Genotype)
- **0/0** = 杂合 (两条染色体都是参考碱基) — "野生型"
- **0/1** = 杂合 (一条参考，一条变异) — 携带变异
- **1/1** = 纯合 (两条都是变异碱基) — 同合变异

#### gVCF vs VCF
- **VCF**: 只包含变异位点
  ```
  #CHROM  POS    REF  ALT  GT:GQ
  chr1    1000   G    A    0/1:50
  chr1    2000   T    C    1/1:99
  ```

- **gVCF** (Genomic VCF): 包含所有位点 (包括没有变异的)
  ```
  #CHROM  POS    REF  ALT  GT:GQ
  chr1    1000   G    .    0/0:99  ← 没有变异
  chr1    2000   T    C    1/1:99  ← 有变异
  ```
  **为什么用 gVCF?** — 多样本联合分型时，需要知道"没有变异"的位点，才能准确判断样本间的一致性。

#### Depth (覆盖深度)
- **定义**: 某位点有多少 reads 覆盖
- **意义**:
  - 深度低 (<5×) — 该位点不可信
  - 深度正常 (20-100×) — 好的覆盖
  - 深度过高 (>200×) — 可能是重复区或污染
- **30×** = 常用的目标深度 (Whole Genome)

---

## 常用参数说明

### 运行流程时最常调整的参数

| 参数 | 默认值 | 说明 | 何时改 |
|------|-------|------|-------|
| `--input` | - | 样本输入表 (CSV) | 必填 |
| `--outdir` | `./results` | 输出目录 | 想改输出位置 |
| `--genome` | `GRCh38` | 参考基因组版本 | 用 hg37 时改为 `GRCh37` |
| `--tools` | `haplotypecaller` | 变异检测工具 | 想用 DeepVariant 改为 `deepvariant` |
| `--max_cpus` | 16 | 最大 CPU 数 | 服务器 CPU 不够时减少 |
| `--max_memory` | 128GB | 最大内存 | 服务器内存不足时减少 |
| `-resume` | - | 断点续跑 | 流程中断后重新开始 |
| `--step` | - | 从某步开始 | 只想跑变异检测，不重新比对 |

### 如何设置参数

**方式 1: 命令行指定**
```bash
nextflow run sarek.nf \
  --input samplesheet.csv \
  --genome GRCh37 \
  --tools deepvariant \
  --max_cpus 8
```

**方式 2: 配置文件**
```groovy
// nextflow.config
params {
    input = 'samplesheet.csv'
    genome = 'GRCh37'
    tools = 'deepvariant'
    max_cpus = 8
}
```

### 工具特定参数

#### HaplotypeCaller 参数
| 参数 | 说明 |
|------|------|
| `--standard-min-confidence-threshold-for-calling` | 变异置信度阈值 (默认 30) — 越低越敏感 |
| `--dbsnp` | 已知 SNP 数据库 (用于标注) |
| `--dragstr-params-path` | DragStr 模型 (优化短串联重复序列检测) |

#### DeepVariant 参数
| 参数 | 说明 |
|------|------|
| `--model-type` | `WGS` (全基因组) / `WES` (全外显子) / `PACBIO` |
| `--custom-model` | 用自定义的深度学习模型 |

---

## 三种变异检测模式

### 1️⃣ Germline (生殖系) — 最常用

**适用场景**:
- 单个个体的遗传变异检测
- 家系研究、群体遗传学
- 疾病相关遗传因子发现

**特点**:
- 检测个人**遗传**的变异 (存在于所有细胞)
- 预期杂合度: 50% (一条染色体 = 0.5 深度)
- 输出: VCF + gVCF

**工具选择**:
```
--tools haplotypecaller    # 推荐，准确率高
--tools deepvariant        # 推荐，特别适合有复杂区域的样本
```

**代码位置**:
```
subworkflows/local/bam_variant_calling_haplotypecaller/main.nf
subworkflows/local/bam_variant_calling_deepvariant/main.nf
```

**流程图**:
```
BAM 文件 (比对结果)
    ↓
HaplotypeCaller (局部重组装)
    ├→ 构建候选变异
    ├→ 贝叶斯推断
    └→ 基因型分型 (0/0, 0/1, 1/1)
    ↓
gVCF 文件 (所有位点)
    ↓
[可选] 多样本联合分型
    ↓
最终 VCF (变异注释)
```

### 2️⃣ Somatic (体细胞) — Paired 模式 — 癌症/肿瘤

**适用场景**:
- 配对样本: 正常组织 + 肿瘤组织
- 癌症驱动基因发现
- 体细胞突变检测

**特点**:
- 检测**仅在肿瘤细胞中**的变异 (正常组织没有)
- 预期杂合度: 取决于肿瘤纯度 (如 50% 纯度 = 0.25 深度)
- 需要两个样本: **normal** + **tumor**
- 输出: VCF (体细胞变异)

**工具选择**:
```
--tools strelka2           # 推荐 SNP/InDel
--tools manta              # 推荐 SV (大型重排)
--tools mutect2            # GATK 体细胞检测
```

**代码位置**:
```
subworkflows/local/bam_variant_calling_somatic/main.nf
```

**流程图**:
```
┌─────────────────────────┬────────────────────┐
│   正常样本 BAM          │   肿瘤样本 BAM      │
│  (背景/对照)             │  (关注)              │
└──────────┬──────────────┴────────────┬───────┘
           │                           │
           └───────────────┬───────────┘
                           ↓
                    Strelka2 / MuTect2
              (配对体细胞变异检测工具)
                    ↙         ↘
        在肿瘤中找到      在正常中找到
        但不在正常中      的位点(忽略)
                    ↓
            体细胞 VCF 文件
    (只包含肿瘤特异的变异)
                    ↓
            变异注释 + 驱动基因预测
```

**关键参数**:
```
--pair_samples true        # 启用配对分析
--tools strelka2,manta     # 可用逗号组合多个工具
```

### 3️⃣ Somatic (体细胞) — Tumor-only 模式 — 单肿瘤

**适用场景**:
- 只有肿瘤样本，没有正常对照
- 临床诊断情景
- 成本控制 (不需要测序正常样本)

**特点**:
- **无法区分**良性变异 vs 驱动变异 (没有对照)
- 预期杂合度: 任意 (取决于克隆性)
- 输出: VCF (但需要更严格的过滤)

**工具选择**:
```
--tools strelka2           # 支持 tumor-only
--tools mutect2            # 有特殊 tumor-only 模式
```

**流程图**:
```
    肿瘤样本 BAM
         ↓
  Strelka2 (Tumor-only)
    或 MuTect2 (TLOD)
         ↓
   VCF 文件 (包含所有变异)
         ↓
  [必须!] 过滤和注释
  - 群体频率过滤 (去除常见变异)
  - 驱动基因数据库标注
  - 临床相关性评分
         ↓
  高置信的肿瘤变异
```

**关键限制**:
- ⚠️ 无法检测 **clone 2** 的变异 (2 个独立克隆时)
- ⚠️ 需要**更严格的过滤**
- ✅ 成本低，快速诊断

---

## 工作流执行命令示例

### 基础 Germline 分析
```bash
nextflow run sarek/workflows/sarek.nf \
  --input samplesheet.csv \
  --genome GRCh38 \
  --tools haplotypecaller \
  --outdir results/ \
  -profile docker
```

### Germline + DeepVariant (推荐准确度高)
```bash
nextflow run sarek/workflows/sarek.nf \
  --input samplesheet.csv \
  --genome GRCh38 \
  --tools deepvariant \
  --outdir results/ \
  -profile docker
```

### 配对肿瘤分析 (Somatic)
```bash
nextflow run sarek/workflows/sarek.nf \
  --input samplesheet.csv \
  --genome GRCh38 \
  --tools strelka2,manta \
  --outdir results/ \
  -profile docker
```

### 只运行变异检测步骤 (跳过比对)
```bash
nextflow run sarek/workflows/sarek.nf \
  --input samplesheet.csv \
  --step variant_calling \
  --genome GRCh38 \
  --tools haplotypecaller
```

### 使用小数据快速测试
```bash
nextflow run sarek/workflows/sarek.nf \
  -profile docker,test
```

---

## 常见问题速查

| 问题 | 答案 | 相关文件 |
|------|------|---------|
| 如何修改 HaplotypeCaller 参数? | 编辑 `conf/modules.config` 的 `GATK4_HAPLOTYPECALLER` 部分 | `modules/nf-core/gatk4/haplotypecaller/main.nf` |
| reads 比对到哪里了? | 看 `subworkflows/local/fastq_align/main.nf` 的 bam.mix() 部分 | `subworkflows/local/fastq_align/main.nf` |
| DeepVariant 为什么比 HC 慢? | CNN 推断需要更多计算，但准确度更高 | `subworkflows/local/bam_variant_calling_deepvariant/main.nf` |
| 如何只跑某个样本? | 在输入 CSV 中只包含那个样本 | `assets/schema_input.json` |
| VCF 文件多大? | 取决于覆盖深度和基因组大小，通常 1-10 GB | - |

---

## 进阶阅读

### 深入理解 Nextflow 语法
- 官方文档: https://www.nextflow.io/docs/latest/
- nf-core 指南: https://nf-co.re/

### 理解变异检测原理
- HaplotypeCaller 论文: "From FastQ to High-Confidence Variant Calls: the Genome Analysis Toolkit Best Practices Pipeline"
- DeepVariant 论文: "A Universal SNP and Small-Indel Variant Caller Using Deep Convolutional Neural Networks"
- GATK 最佳实践: https://gatk.broadinstitute.org/hc/en-us

### Sarek 特定资源
- GitHub: https://github.com/nf-core/sarek
- 官方文档: https://nf-co.re/sarek
- 讨论区: https://github.com/nf-core/sarek/discussions

---

**祝你的基因组分析顺利!** 🧬

有问题时，记住三个关键文件:
1. 📄 `workflows/sarek.nf` — 了解全局流程
2. 📄 `subworkflows/local/*/main.nf` — 了解某个步骤
3. 📄 `conf/modules.config` — 了解参数配置
