// ╔══════════════════════════════════════════════════════════════════╗
// ║  HaplotypeCaller Variant Calling Subworkflow                    ║
// ║  HaplotypeCaller 变异检测子流程                                   ║
// ╚══════════════════════════════════════════════════════════════════╝
// GATK4 HaplotypeCaller performs germline (non-somatic) variant calling
// GATK4 的 HaplotypeCaller 工具用于检测生殖系 (non-somatic) 变异
// Methodology: Local assembly + Bayesian inference for variant genotyping
// 方法: 局部重组装 reads + 贝叶斯推断进行变异分型
// Outputs both gVCF (for joint genotyping) and final VCF files
// 同时输出 gVCF (用于多样本联合分型) 和最终 VCF 文件

//
// GATK4 HAPLOTYPACALLER germline variant calling:
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BAM_MERGE_INDEX_SAMTOOLS                            } from '../bam_merge_index_samtools/main'
include { GATK4_HAPLOTYPECALLER                               } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS            as MERGE_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/mergevcfs/main'

// ╔══════════════════════════════════════════════════════════════════╗
// ║  Main Workflow: BAM_VARIANT_CALLING_HAPLOTYPECALLER             ║
// ║  主工作流程: 用 HaplotypeCaller 从 BAM 文件检测变异              ║
// ╚══════════════════════════════════════════════════════════════════╝
// Takes aligned BAM files and outputs variant call files (VCF/gVCF)
// 接收比对后的 BAM 文件，输出变异检测结果 (VCF/gVCF 文件)

workflow BAM_VARIANT_CALLING_HAPLOTYPECALLER {
    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Input Channels (输入数据流)                                  ║
    // ╚══════════════════════════════════════════════════════════════╝
    take:
    cram                         // channel: [mandatory] [ meta, cram, crai, interval.bed ]
                                 // 压缩的比对文件 (CRAM) - 比 BAM 更节省空间，包含 meta 信息和区间文件
    fasta                        // channel: [mandatory]
                                 // 参考基因组 FASTA 文件
    fasta_fai                    // channel: [mandatory]
                                 // FASTA 索引文件 (加速对参考序列的查询)
    dict                         // channel: [mandatory]
                                 // GATK SAM 序列字典文件 (包含染色体信息和长度)
    dbsnp                        // channel: [optional]
                                 // dbSNP VCF 文件 (已知 SNP 数据库，用于 variant recalibration)
    dbsnp_tbi                    // channel: [optional]
                                 // dbSNP VCF 索引文件
    intervals                    // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
                                 // 基因组区间文件 - 用于并行化处理 (每个区间独立分析以加速)
                                 // Interval splitting strategy: 将全基因组分割成若干块，每块独立运行
                                 // 如果 num_intervals=1，整个基因组作为一个区间处理

    main:
    versions = Channel.empty()

    vcf           = Channel.empty()
    realigned_bam = Channel.empty()

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Prepare Data for Parallel Processing (准备并行处理的数据)    ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Strategy: Spread-and-Gather Pattern (分散-聚合策略)          ║
    // ║                                                               ║
    // ║ 基本思想: 将基因组分成多个区间，每个区间独立运行             ║
    // ║ HaplotypeCaller，最后合并结果                                ║
    // ║ This enables parallel processing for faster analysis          ║
    // ║ 这种并行化策略可以大大加速分析                                ║
    // ╚══════════════════════════════════════════════════════════════╝

    // Combine cram and intervals for spread and gather strategy
    // 将 CRAM 文件和基因组区间组合在一起，实现 spread-and-gather 策略
    cram_intervals = cram.combine(intervals)
        // ╔══════════════════════════════════════════════════════════╗
        // ║  Update Metadata Map (更新样本元信息)                    ║
        // ╚══════════════════════════════════════════════════════════╝
        // Move num_intervals to meta map
        // 将"区间数量"(num_intervals) 字段添加到 meta 字典中
        // Add interval_name to allow correct merging with interval files
        // 添加"区间名"(interval_name) 字段，用于最后合并阶段识别和匹配分割的 VCF
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ interval_name:intervals.baseName, num_intervals:num_intervals, variantcaller:'haplotypecaller' ], cram, crai, intervals, [] ] }

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Run HaplotypeCaller (执行 HaplotypeCaller)                  ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Process:                                                     ║
    // ║ 1. Local reassembly: 在每个区间内局部重组装 reads            ║
    // ║ 2. Bayesian inference: 使用贝叶斯模型推断变异位点            ║
    // ║ 3. Genotyping: 对每个样本进行基因型分型                      ║
    // ║                                                               ║
    // ║ Output: gVCF (genomic VCF) with all sites, not just variants ║
    // ║ 输出: gVCF (包含所有位点，不仅仅是变异位点)                  ║
    // ╚══════════════════════════════════════════════════════════════╝
    GATK4_HAPLOTYPECALLER(
        cram_intervals,
        fasta,
        fasta_fai,
        dict,
        dbsnp,
        dbsnp_tbi)

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Prepare gVCF for Joint Genotyping (准备用于联合分型的 gVCF)  ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ gVCF (Genomic VCF): Includes ALL sites, even non-variants    ║
    // ║ gVCF 文件包含所有位点，包括没有变异的位点                    ║
    // ║ Essential for multi-sample joint genotyping                  ║
    // ║ 多样本联合分型时必须使用 gVCF                                ║
    // ║                                                               ║
    // ║ Example gVCF content:                                        ║
    // ║   chr1 1000 . A . . PASS ... GT:DP:GQ 0/0:50:99             ║
    // ║   chr1 2000 . G A . PASS ... GT:DP:GQ 0/1:45:88             ║
    // ║ (first line shows "." = no variant; second shows variant)   ║
    // ╚══════════════════════════════════════════════════════════════╝
    // For joint genotyping
    // 这部分输出用于后续的多样本联合分型步骤
    gvcf_tbi_intervals = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnMismatch: true)
        .join(cram_intervals, failOnMismatch: true)
        .map{ meta, gvcf, tbi, cram, crai, intervals, dragstr_model -> [ meta, gvcf, tbi, intervals ] }

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Handling Interval Splitting (处理区间分割)                   ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Branch VCF outputs into two categories:                      ║
    // ║ 根据区间数量将 VCF 输出分为两类:                              ║
    // ║                                                               ║
    // ║ intervals:    num_intervals > 1   (多个区间 → 需要合并)       ║
    // ║ no_intervals: num_intervals <= 1  (单个区间 → 无需合并)       ║
    // ║                                                               ║
    // ║ This allows conditional processing: if multiple regions,    ║
    // ║ merge them; if single region, use directly                  ║
    // ╚══════════════════════════════════════════════════════════════╝

    // Figuring out if there is one or more vcf(s) from the same sample
    // 判断同一样本是否产生了多个 VCF (对应多个基因组区间) 或仅一个 VCF
    haplotypecaller_vcf = GATK4_HAPLOTYPECALLER.out.vcf.map{
            meta, vcf -> [ meta - meta.subMap('interval_name'), vcf]
        }
        .branch{
        // Use meta.num_intervals to asses number of intervals
        // 根据 num_intervals 字段判断是否产生了多个区间的 VCF
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Figuring out if there is one or more tbi(s) from the same sample
    // 同样处理 TBI (VCF 索引文件)
    haplotypecaller_tbi = GATK4_HAPLOTYPECALLER.out.tbi.map{
            meta, tbi -> [ meta - meta.subMap('interval_name'), tbi]
        }.branch{
        // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Figuring out if there is one or more bam(s) from the same sample
    // HaplotypeCaller 在 local reassembly 过程中也生成中间 BAM 文件
    // 同样对这些 BAM 文件进行分支处理
    haplotypecaller_bam = GATK4_HAPLOTYPECALLER.out.bam.map{
            meta, bam -> [ meta - meta.subMap('interval_name'), bam]
        }.branch{
        // Use meta.num_intervals to asses number of intervals
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Merge VCF Files from Multiple Intervals (合并多区间 VCF)    ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ When num_intervals > 1:                                      ║
    // ║ 当一个样本的变异检测结果分散在多个区间时:                    ║
    // ║                                                               ║
    // ║ 1. groupKey(meta, num_intervals): 创建分组键，按样本分组    ║
    // ║ 2. groupTuple(): 将相同样本的 VCF 文件合并为一个列表         ║
    // ║ 3. MERGE_HAPLOTYPECALLER: 合并多个 VCF 成单一文件            ║
    // ║                                                               ║
    // ║ Example: Sample1 产生了 3 个区间 VCF                          ║
    // ║   Sample1_chr1.vcf.gz, Sample1_chr2.vcf.gz, Sample1_chr3... ║
    // ║ → 合并成 Sample1_final.vcf.gz                                ║
    // ╚══════════════════════════════════════════════════════════════╝

    // Only when using intervals
    // 仅当划分了多个基因组区间时才执行合并步骤
    MERGE_HAPLOTYPECALLER(haplotypecaller_vcf.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ] }.groupTuple(), dict)

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Consolidate Final VCF Output (统一最终 VCF 输出)            ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Mix two sources:                                             ║
    // ║ 1. Merged VCF (from multiple intervals)                      ║
    // ║ 2. Single VCF (from no_intervals)                            ║
    // ║                                                               ║
    // ║ Result: Regardless of interval split, each sample has one   ║
    // ║ final VCF file for downstream analysis                      ║
    // ╚══════════════════════════════════════════════════════════════╝
    haplotypecaller_vcf = Channel.empty().mix(
            MERGE_HAPLOTYPECALLER.out.vcf,
            haplotypecaller_vcf.no_intervals)

    haplotypecaller_tbi = Channel.empty().mix(
            MERGE_HAPLOTYPECALLER.out.tbi,
            haplotypecaller_tbi.no_intervals)

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Merge and Index Realigned BAM Files (合并和索引重组装 BAM)   ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ HaplotypeCaller performs local reassembly, generating        ║
    // ║ realigned BAM files that contain locally re-aligned reads   ║
    // ║ These are optionally merged and indexed for inspection       ║
    // ║                                                               ║
    // ║ HaplotypeCaller 在局部重组装过程中生成中间 BAM 文件           ║
    // ║ 这些 BAM 文件可用于检查比对的质量                             ║
    // ╚══════════════════════════════════════════════════════════════╝
    // BAM output
    BAM_MERGE_INDEX_SAMTOOLS(haplotypecaller_bam.intervals
        .map{ meta, bam -> [ groupKey(meta, meta.num_intervals), bam ] }
        .groupTuple()
        .mix(haplotypecaller_bam.no_intervals))

    realigned_bam = BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Collect Software Versions (记录软件版本)                    ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Track versions for reproducibility and troubleshooting       ║
    // ║ 记录所有使用的工具版本，确保结果可重复                        ║
    // ╚══════════════════════════════════════════════════════════════╝
    versions = versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
    versions = versions.mix(MERGE_HAPLOTYPECALLER.out.versions)

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Clean Up Metadata (清理元数据)                              ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Remove temporary fields used during interval processing      ║
    // ║ 删除在区间处理过程中临时添加的字段                            ║
    // ║ Result: Clean meta map for downstream tools                  ║
    // ║ 结果: 清晰的元数据，可用于下游分析工具                       ║
    // ╚══════════════════════════════════════════════════════════════╝
    // Remove no longer necessary field: num_intervals
    vcf = haplotypecaller_vcf.map{ meta, vcf -> [ meta - meta.subMap('num_intervals'), vcf ] }
    tbi = haplotypecaller_tbi.map{ meta, tbi -> [ meta - meta.subMap('num_intervals'), tbi ] }

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Emit Output Channels (输出数据流)                            ║
    // ╚══════════════════════════════════════════════════════════════╝
    emit:
    gvcf_tbi_intervals // For joint genotyping (用于多样本联合分型)
                       // gVCF 包含所有位点，允许与其他样本一起进行联合分型
    realigned_bam      // Optional (可选)
                       // HaplotypeCaller 局部重组装后的 BAM 文件
    vcf                // vcf (Final variant call file)
                       // 最终的变异检测结果 VCF 文件
    tbi                // tbi (VCF tabix index)
                       // VCF 索引文件，加速查询

    versions           // Software versions for reproducibility
                       // 所有工具的版本信息
}
