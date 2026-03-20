// ╔══════════════════════════════════════════════════════════════════╗
// ║  FASTQ Alignment Subworkflow                                     ║
// ║  FASTQ 序列比对子流程                                             ║
// ╚══════════════════════════════════════════════════════════════════╝
// This subworkflow aligns raw sequencing reads (FASTQ) to a reference genome
// 该子流程将原始测序reads (FASTQ文件) 比对到参考基因组上
// Supported aligners: BWA-MEM, BWA-MEM2, DRAGMAP, Sentieon
// 支持的比对工具: BWA-MEM, BWA-MEM2, DRAGMAP, Sentieon

//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { DRAGMAP_ALIGN          } from '../../../modules/nf-core/dragmap/align/main'
include { SENTIEON_BWAMEM        } from '../../../modules/nf-core/sentieon/bwamem/main'

// ╔══════════════════════════════════════════════════════════════════╗
// ║  Main Workflow: FASTQ_ALIGN                                      ║
// ║  主工作流程: FASTQ_ALIGN                                          ║
// ╚══════════════════════════════════════════════════════════════════╝
// This workflow takes raw FASTQ reads and aligns them to a reference genome
// 该工作流程接收原始 FASTQ reads，并使用指定的比对工具将其对标到参考基因组
// Produces: BAM files containing the aligned reads (输出: 包含比对结果的 BAM 文件)

workflow FASTQ_ALIGN {
    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Input Channels (输入数据流)                                  ║
    // ╚══════════════════════════════════════════════════════════════╝
    take:
    reads // channel: [mandatory] meta, reads
          // 包含样本元信息和测序 reads 文件的数据流
    index // channel: [mandatory] index
          // 参考基因组的预建索引文件 (已通过 bwa index 或 bwamem2 index 等生成)
    sort  // boolean: [mandatory] true -> sort, false -> don't sort
          // 是否对输出的 BAM 文件进行排序 (true=按基因组坐标排序; false=不排序)
    fasta // 参考基因组 FASTA 文件 (核酸序列数据库)
    fasta_fai // FASTA 索引文件 (加速基因组访问)

    main:

    versions = Channel.empty()
    reports = Channel.empty()

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Aligner Selection (选择比对工具)                             ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Only ONE aligner will run based on configuration             ║
    // ║ 根据配置文件只会运行其中一个比对工具                           ║
    // ║                                                               ║
    // ║ BWA-MEM: Seed-and-extend algorithm for short reads            ║
    // ║ 使用种子延伸算法将短 reads 比对到参考基因组                   ║
    // ║ Most commonly used for WGS (全基因组测序最常用)               ║
    // ╚══════════════════════════════════════════════════════════════╝

    // If aligner is bwa-mem (legacy version)
    // 如果选择使用 BWA-MEM (较早版本)
    BWAMEM1_MEM(reads, index, [[id:'no_fasta'], []], sort)

    // If aligner is bwa-mem2 (optimized version for modern CPUs)
    // 如果选择使用 BWA-MEM2 (针对现代 CPU 优化的版本)
    BWAMEM2_MEM(reads, index, [[id:'no_fasta'], []], sort)

    // If aligner is dragmap (fast alignment tool from Illumina)
    // 如果选择使用 DRAGMAP (Illumina 的快速比对工具)
    DRAGMAP_ALIGN(reads, index, [[id:'no_fasta'], []], sort)

    // The sentieon-bwamem-module does sorting as part of the conversion from sam to bam.
    // Sentieon BWAMEM 会在 SAM 转换为 BAM 的过程中自动进行排序
    // If using Sentieon (accelerated variant calling with commercial licensing)
    // 如果选择使用 Sentieon (商业加速变异检测)
    SENTIEON_BWAMEM(reads, index, fasta, fasta_fai)

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Collect Outputs from Aligners (汇总所有比对工具的输出)      ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Get the bam files from the aligner                           ║
    // ║ 从各个比对工具提取 BAM 文件 (只有激活的工具会产生实际输出)   ║
    // ╚══════════════════════════════════════════════════════════════╝
    // Only one aligner is run
    bam = Channel.empty()
    bam = bam.mix(BWAMEM1_MEM.out.bam)
    bam = bam.mix(BWAMEM2_MEM.out.bam)
    bam = bam.mix(DRAGMAP_ALIGN.out.bam)
    bam = bam.mix(SENTIEON_BWAMEM.out.bam_and_bai.map{ meta, bam, bai -> [ meta, bam ] })

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  BAM Index (.bai) Files (BAM 索引文件)                       ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ BAI is a binary index for quick access to specific regions   ║
    // ║ .bai 索引文件加速对 BAM 文件中特定基因组区域的查询           ║
    // ║ Essential for downstream analysis tools                      ║
    // ║ 下游分析工具 (如 IGV 查看器) 必需                             ║
    // ╚══════════════════════════════════════════════════════════════╝
    bai = SENTIEON_BWAMEM.out.bam_and_bai.map{ meta, bam, bai -> [ meta, bai ] }

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Collect QC Reports (汇总质控报告)                            ║
    // ║ Alignment logs and statistics for quality assessment         ║
    // ║ 比对日志和统计信息，用于质量评估                              ║
    // ╚══════════════════════════════════════════════════════════════╝
    // Gather reports of all tools used
    reports = reports.mix(DRAGMAP_ALIGN.out.log)

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Track Software Versions (记录软件版本)                      ║
    // ║ Essential for reproducibility and troubleshooting            ║
    // ║ 对于结果可重复性和问题排查很重要                              ║
    // ╚══════════════════════════════════════════════════════════════╝
    // Gather versions of all tools used
    versions = versions.mix(BWAMEM1_MEM.out.versions)
    versions = versions.mix(BWAMEM2_MEM.out.versions)
    versions = versions.mix(DRAGMAP_ALIGN.out.versions)
    versions = versions.mix(SENTIEON_BWAMEM.out.versions)

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Emit Output Channels (输出数据流)                            ║
    // ╚══════════════════════════════════════════════════════════════╝
    emit:
    bam      // channel: [ [meta], bam ] — Aligned reads in BAM format
             // BAM 格式的比对结果 (包含 read 在基因组上的位置和相关信息)
    bai      // channel: [ [meta], bai ] — BAM index for quick access
             // BAM 索引文件 (加速查询)
    reports  // Quality control logs and statistics
             // 质控日志和统计信息
    versions // channel: [ versions.yml ] — Software version tracking
             // 记录所有使用工具的版本号
}
