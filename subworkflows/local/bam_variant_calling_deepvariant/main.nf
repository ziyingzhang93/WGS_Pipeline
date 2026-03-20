// ╔══════════════════════════════════════════════════════════════════╗
// ║  DeepVariant Variant Calling Subworkflow                        ║
// ║  DeepVariant 变异检测子流程                                       ║
// ╚══════════════════════════════════════════════════════════════════╝
// DeepVariant is a deep learning (CNN-based) variant caller developed by Google
// DeepVariant 是由 Google 开发的基于深度学习 (卷积神经网络) 的变异检测工具
// Methodology: Converts alignments to images → CNN inference → variant calls
// 方法: 将 reads 堆叠转换为图像 → CNN 推断 → 变异检测结果
// Advantage: High accuracy across different sequencing technologies and depths
// 优势: 对不同测序技术和深度都有很高的准确性
// Reference: https://github.com/google/deepvariant/issues/510

//
// DEEPVARIANT germline calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { DEEPVARIANT_RUNDEEPVARIANT                } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_GVCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_VCF  } from '../../../modules/nf-core/gatk4/mergevcfs/main'

// Deepvariant: https://github.com/google/deepvariant/issues/510
// ╔══════════════════════════════════════════════════════════════════╗
// ║  Main Workflow: BAM_VARIANT_CALLING_DEEPVARIANT                 ║
// ║  主工作流程: 用 DeepVariant 从 BAM 文件检测变异                  ║
// ╚══════════════════════════════════════════════════════════════════╝
// Takes aligned BAM/CRAM files and outputs variant calls using DeepVariant
// 接收比对后的 BAM/CRAM 文件，使用 DeepVariant 输出变异检测结果

workflow BAM_VARIANT_CALLING_DEEPVARIANT {
    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Input Channels (输入数据流)                                  ║
    // ╚══════════════════════════════════════════════════════════════╝
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
                  // 压缩比对文件 (CRAM) - 比 BAM 更节省空间
                  // Contains reads aligned to reference genome
    dict          // channel: [optional]  [ meta, dict ]
                  // GATK SAM 序列字典 (包含染色体信息)
    fasta         // channel: [mandatory] [ fasta ]
                  // 参考基因组 FASTA 文件 (核酸序列)
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
                  // FASTA 索引文件 (加速序列查询)
    intervals     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
                  // 基因组区间文件 - 用于并行化处理
                  // Interval splitting strategy: 将全基因组分成多个块，每块独立运行
                  // DeepVariant，最后合并结果

    main:
    versions = Channel.empty()

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Prepare Data for Parallel Processing (准备并行处理的数据)    ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Strategy: Spread-and-Gather Pattern (分散-聚合策略)          ║
    // ║                                                               ║
    // ║ 基本思想: 将基因组分成多个区间，每个区间独立运行             ║
    // ║ DeepVariant，最后合并结果                                    ║
    // ║ This enables parallel processing for faster analysis          ║
    // ║ 这种并行化策略可以大大加速分析                                ║
    // ║ (尤其是全基因组测序，通常有 50-100+ 个区间)                  ║
    // ╚══════════════════════════════════════════════════════════════╝

    // Combine cram and intervals for spread and gather strategy
    // 将 CRAM 文件和基因组区间组合在一起，实现 spread-and-gather 策略
    cram_intervals = cram.combine(intervals)
        // ╔══════════════════════════════════════════════════════════╗
        // ║  Update Metadata Map (更新样本元信息)                    ║
        // ╚══════════════════════════════════════════════════════════╝
        // Move num_intervals to meta map
        // 将"区间数量"(num_intervals) 字段添加到 meta 字典中
        // 用于后续判断是否需要合并多个区间的结果
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram, crai, intervals ]}

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Run DeepVariant (执行 DeepVariant)                          ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Process:                                                     ║
    // ║ 1. Make Examples: 将 reads 堆叠转换为 TensorFlow examples   ║
    // ║ 2. Call Variants: CNN 推断，预测变异位点                    ║
    // ║ 3. Post-processing: 后处理，生成 VCF/gVCF 文件              ║
    // ║                                                               ║
    // ║ Output: VCF (final variant calls) + gVCF (all sites)         ║
    // ║ 输出: VCF (最终变异结果) + gVCF (所有位点)                   ║
    // ║                                                               ║
    // ║ Why images? CNN can learn patterns in read pileups:          ║
    // ║ 为什么用图像? 卷积神经网络可以学习 reads 堆叠中的模式:      ║
    // ║   - Depth patterns (覆盖深度模式)                            ║
    // ║   - Base calling quality (碱基质量模式)                      ║
    // ║   - Read strand bias (reads 链偏差)                          ║
    // ║   - Flanking sequences (两侧序列)                            ║
    // ║ This makes DeepVariant more robust to different data types  ║
    // ║ 这使得 DeepVariant 对不同类型的数据更加稳健                  ║
    // ╚══════════════════════════════════════════════════════════════╝

    DEEPVARIANT_RUNDEEPVARIANT(cram_intervals, fasta, fasta_fai, [ [ id:'null' ], [] ], [ [ id:'null' ], [] ])

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Handling Interval Splitting - VCF (处理区间分割 - VCF)     ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Branch VCF outputs into two categories:                      ║
    // ║ 根据区间数量将 VCF 输出分为两类:                              ║
    // ║                                                               ║
    // ║ intervals:    num_intervals > 1   (多个区间 → 需要合并)       ║
    // ║ no_intervals: num_intervals <= 1  (单个区间 → 无需合并)       ║
    // ║                                                               ║
    // ║ This allows conditional processing: if multiple regions,    ║
    // ║ merge them; if single region, use directly                  ║
    // ║                                                               ║
    // ║ Example: WGS 测序通常将基因组分成 100+ 个区间                ║
    // ║ DeepVariant 为每个区间输出一个 VCF 文件                      ║
    // ║ → 需要合并成一个最终 VCF 文件                                ║
    // ╚══════════════════════════════════════════════════════════════╝

    // Figuring out if there is one or more vcf(s) from the same sample
    // 判断同一样本是否产生了多个 VCF (对应多个基因组区间) 或仅一个 VCF
    vcf_out = DEEPVARIANT_RUNDEEPVARIANT.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        // 根据 num_intervals 字段判断是否产生了多个区间的 VCF
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Handling Interval Splitting - gVCF (处理区间分割 - gVCF)   ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ gVCF (Genomic VCF): Includes ALL sites, even non-variants    ║
    // ║ gVCF 文件包含所有位点，包括没有变异的位点                    ║
    // ║                                                               ║
    // ║ Example gVCF content:                                        ║
    // ║   chr1 1000 . A . . PASS ... GT:DP:GQ 0/0:50:99             ║
    // ║   chr1 2000 . G A . PASS ... GT:DP:GQ 0/1:45:88             ║
    // ║ (first = no variant, second = heterozygous variant)         ║
    // ║                                                               ║
    // ║ Essential for multi-sample joint genotyping                  ║
    // ║ 多样本联合分型时必须使用 gVCF                                ║
    // ╚══════════════════════════════════════════════════════════════╝

    // Figuring out if there is one or more gvcf(s) from the same sample
    // 同样处理 gVCF 文件
    gvcf_out = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Prepare Channels for Merging (准备合并用的数据流)            ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ groupKey: 创建分组键，按样本分组所有的区间 VCF/gVCF         ║
    // ║ groupTuple: 将相同样本的所有 VCF 文件合并为一个列表          ║
    // ║                                                               ║
    // ║ Example processing for one sample with 3 intervals:          ║
    // ║   Input:  Sample1_interval1.vcf, Sample1_interval2.vcf...    ║
    // ║ → groupKey(Sample1, 3)                                       ║
    // ║ → groupTuple([Sample1_interval1, Sample1_interval2...])      ║
    // ║ → MERGE outputs single file: Sample1.vcf                     ║
    // ╚══════════════════════════════════════════════════════════════╝

    // Only when using intervals
    // 仅当划分了多个基因组区间时才执行合并步骤
    gvcf_to_merge = gvcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    vcf_to_merge = vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Merge gVCF and VCF Files (合并 gVCF 和 VCF 文件)            ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ GATK4_MERGEVCFS: 标准 VCF 合并工具                           ║
    // ║ Concatenates VCF records from multiple regions into single   ║
    // ║ 将多个区间的 VCF 记录连接成一个文件                          ║
    // ║                                                               ║
    // ║ Requirements:                                                ║
    // ║ - Input VCF must be sorted and indexed (bgzipped + tabix)   ║
    // ║ - Output is sorted and re-indexed                           ║
    // ║ - dict: 用于验证染色体顺序和长度                              ║
    // ╚══════════════════════════════════════════════════════════════╝

    MERGE_DEEPVARIANT_GVCF(gvcf_to_merge, dict)
    MERGE_DEEPVARIANT_VCF(vcf_to_merge, dict)

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Handling Interval Splitting - TBI (处理区间分割 - TBI 索引) ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ TBI (Tabix Index): Binary index for bgzipped VCF files       ║
    // ║ TBI 是 bgzip 压缩 VCF 文件的二进制索引                        ║
    // ║ Enables quick lookup of specific genomic regions             ║
    // ║ 加速对特定基因组区域的查询                                    ║
    // ║                                                               ║
    // ║ VCF 文件必须被 bgzip 压缩，然后 tabix 索引才能使用           ║
    // ║ 这个索引对 IGV、VEP 等下游工具是必需的                        ║
    // ╚══════════════════════════════════════════════════════════════╝

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_out = DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Consolidate Final Output (统一最终输出)                     ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Mix two sources:                                             ║
    // ║ 1. Merged VCF/gVCF (from multiple intervals)                 ║
    // ║ 2. Single VCF/gVCF (from no_intervals)                       ║
    // ║                                                               ║
    // ║ Result: Regardless of interval split, each sample has one   ║
    // ║ final VCF and gVCF file for downstream analysis             ║
    // ║ 结果: 不管是否分割基因组，每个样本都有一个最终 VCF 和 gVCF ║
    // ║                                                               ║
    // ║ Add variantcaller field to distinguish from other callers   ║
    // ║ 在 meta 中添加 'variantcaller:deepvariant' 标签              ║
    // ║ 用于后续区分不同的变异检测工具的结果                          ║
    // ╚══════════════════════════════════════════════════════════════╝

    // Mix intervals and no_intervals channels together
    // 混合多区间合并的结果和单区间的结果
    gvcf = Channel.empty().mix(MERGE_DEEPVARIANT_GVCF.out.vcf, gvcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        // 添加 variantcaller 标识，删除不再需要的 num_intervals 字段
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], vcf ] }

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MERGE_DEEPVARIANT_VCF.out.vcf, vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], vcf ] }

    tbi = Channel.empty().mix(MERGE_DEEPVARIANT_VCF.out.tbi, tbi_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, tbi -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], tbi ] }

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Collect Software Versions (记录软件版本)                    ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Track versions for reproducibility and troubleshooting       ║
    // ║ 记录所有使用的工具版本，确保结果可重复                        ║
    // ╚══════════════════════════════════════════════════════════════╝
    versions = versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions)
    versions = versions.mix(MERGE_DEEPVARIANT_GVCF.out.versions)
    versions = versions.mix(MERGE_DEEPVARIANT_VCF.out.versions)

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Emit Output Channels (输出数据流)                            ║
    // ╚══════════════════════════════════════════════════════════════╝
    emit:
    gvcf    // gVCF file with all sites (包含所有位点的 gVCF)
            // Essential for joint genotyping across multiple samples
            // 用于多样本联合分型
    vcf     // Final variant calls (最终变异检测结果 VCF)
            // Contains only variant sites for variant annotation and analysis
            // 只包含变异位点，用于变异注释和分析
    tbi     // VCF tabix index (VCF Tabix 索引)
            // Required for downstream tools like VEP, SnpEff
            // 下游工具 (如 VEP) 需要此索引

    versions // Software versions for reproducibility (软件版本)
             // All tool versions used in this workflow
}
