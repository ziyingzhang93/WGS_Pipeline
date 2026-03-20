// ╔══════════════════════════════════════════════════════════════════════════╗
// ║                        IMPORT MODULES & SUBWORKFLOWS                     ║
// ║                      导入所有必需的模块和子工作流                         ║
// ╚══════════════════════════════════════════════════════════════════════════╝
//
// BIOLOGICAL PURPOSE:
// This section imports all the bioinformatics tools and workflows needed to
// process raw DNA sequencing reads into analysis-ready BAM/CRAM files.
//
// 生物学目的：
// 本部分导入所有将原始DNA测序数据处理成可用于下游分析的BAM/CRAM文件所需的工具

// Create samplesheets to restart from different steps
// 创建样本列表，允许从不同步骤重新启动流程
include { CHANNEL_ALIGN_CREATE_CSV                          } from '../../../subworkflows/local/channel_align_create_csv/main'
include { CHANNEL_MARKDUPLICATES_CREATE_CSV                 } from '../../../subworkflows/local/channel_markduplicates_create_csv/main'
include { CHANNEL_BASERECALIBRATOR_CREATE_CSV               } from '../../../subworkflows/local/channel_baserecalibrator_create_csv/main'
include { CHANNEL_APPLYBQSR_CREATE_CSV                      } from '../../../subworkflows/local/channel_applybqsr_create_csv/main'

// Convert BAM files to FASTQ files
// 将BAM文件转换回FASTQ格式（用于UMI处理）
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_UMI         } from '../../../subworkflows/local/bam_convert_samtools/main'

// ────────────────────────────────────────────────────────────────────────────
// TRIM/SPLIT FASTQ Files
// Removes adapter sequences and low-quality bases from sequencing reads
// FASTP: 修剪/分割FASTQ文件
// 去除测序reads中的接头序列和低质量的碱基
// ────────────────────────────────────────────────────────────────────────────
include { FASTP                                             } from '../../../modules/nf-core/fastp/main'

// Remove genomic contaminants with bbsplit
// Identifies and removes reads that belong to contaminating genomes
// BBSplit: 使用BBMap检测和去除污染reads
// 识别并去除来自污染物基因组的测序reads（如真菌、细菌污染）
include { BBMAP_BBSPLIT                                     } from '../../../modules/nf-core/bbmap/bbsplit'
//TODO: WHAT ABOUT BBSPLIT RUNS WITH PARABRICKS?

// Create umi consensus bams from fastq
// UMI (Unique Molecular Identifier) allows PCR error correction by grouping reads
// from identical original DNA molecules. This creates a high-quality consensus.
// 从FASTQ创建UMI共识BAM
// UMI（分子条形码）允许通过合并来自同一原始DNA分子的reads来纠正PCR错误，
// 生成高质量的共识序列
include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO                  } from '../../../subworkflows/local/fastq_create_umi_consensus_fgbio/main'

// ────────────────────────────────────────────────────────────────────────────
// Map input reads to reference genome
// Aligns sequencing reads to the reference genome using tools like BWA-MEM
// 将reads比对到参考基因组
// 使用BWA-MEM等工具将测序reads对齐到参考基因组
// ────────────────────────────────────────────────────────────────────────────
include { FASTQ_ALIGN                                       } from '../../../subworkflows/local/fastq_align/main'

// Merge and index BAM files (optional)
// Combines BAM files from the same sample (different sequencing lanes)
// 合并和索引BAM文件（可选）
// 合并来自同一样本的多个BAM文件（来自不同测序lanes）
include { BAM_MERGE_INDEX_SAMTOOLS                          } from '../../../subworkflows/local/bam_merge_index_samtools/main'

// Convert BAM files to CRAM files
// CRAM is a compressed format that saves disk space (30-50% smaller than BAM)
// while preserving all quality information needed for variant calling
// 将BAM文件转换为CRAM文件
// CRAM是一种压缩格式，节省磁盘空间（比BAM小30-50%），
// 同时保留变异检测所需的所有质量信息
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING           } from '../../../modules/nf-core/samtools/convert/main'

// Convert CRAM files (optional)
// Converts CRAM back to BAM if user prefers BAM format for downstream analysis
// 转换CRAM文件（可选）
// 如果用户更喜欢BAM格式用于下游分析，将CRAM转换回BAM
include { SAMTOOLS_CONVERT as CRAM_TO_BAM                   } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_RECAL             } from '../../../modules/nf-core/samtools/convert/main'

// ────────────────────────────────────────────────────────────────────────────
// Copy UMIs from read name to RX tag
// Extracts UMI information embedded in read names and properly tags it
// 从read名称复制UMI到RX标签
// 提取嵌入在read名称中的UMI信息并正确标记
// ────────────────────────────────────────────────────────────────────────────
include { FGBIO_COPYUMIFROMREADNAME                         } from '../../../modules/nf-core/fgbio/copyumifromreadname/main'

// ────────────────────────────────────────────────────────────────────────────
// Mark Duplicates (+QC)
// Identifies PCR duplicate reads that arise from library amplification
// PCR duplicates can create false variant calls and must be marked for exclusion
// 标记重复reads（+质量控制）
// 识别来自文库扩增的PCR重复reads，这些重复可能导致假阳性变异检测，需要标记
// Three options: standard GATK, Spark (for large datasets), or Sentieon (optimized)
// ────────────────────────────────────────────────────────────────────────────
include { BAM_MARKDUPLICATES                                } from '../../../subworkflows/local/bam_markduplicates/main'
include { BAM_MARKDUPLICATES_SPARK                          } from '../../../subworkflows/local/bam_markduplicates_spark/main'
include { BAM_SENTIEON_DEDUP                                } from '../../../subworkflows/local/bam_sentieon_dedup/main'

// QC on CRAM
// Quality control metrics on aligned reads (depth, coverage, contamination)
// CRAM质量控制
// 对比对reads进行质量控制检查（深度、覆盖度、污染情况）
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD        } from '../../../subworkflows/local/cram_qc_mosdepth_samtools/main'

// ────────────────────────────────────────────────────────────────────────────
// Create recalibration tables (BQSR - Base Quality Score Recalibration)
// Builds a statistical model of systematic sequencing errors to improve
// quality score accuracy. Uses known variant databases (dbSNP, indels).
// 创建碱基质量分数重新校准（BQSR）表
// 建立测序错误的统计模型来改进质量分数准确性，使用已知变异数据库
// ────────────────────────────────────────────────────────────────────────────
include { BAM_BASERECALIBRATOR                              } from '../../../subworkflows/local/bam_baserecalibrator/main'
include { BAM_BASERECALIBRATOR_SPARK                        } from '../../../subworkflows/local/bam_baserecalibrator_spark/main'

// ────────────────────────────────────────────────────────────────────────────
// Apply recalibrated quality scores to BAM/CRAM files
// Produces final analysis-ready alignment files with corrected quality scores
// 应用校正后的碱基质量分数
// 生成具有校正后质量分数的最终可用于分析的比对文件
// ────────────────────────────────────────────────────────────────────────────
include { BAM_APPLYBQSR                                     } from '../../../subworkflows/local/bam_applybqsr/main'
include { BAM_APPLYBQSR_SPARK                               } from '../../../subworkflows/local/bam_applybqsr_spark/main'


// ╔══════════════════════════════════════════════════════════════════════════╗
// ║          WORKFLOW: FASTQ_PREPROCESS_GATK                                ║
// ║          工作流：FASTQ预处理（GATK方案）                                ║
// ╚══════════════════════════════════════════════════════════════════════════╝
//
// OVERVIEW:
// This workflow transforms raw sequencing data into variant-calling-ready
// alignment files through 4 major preprocessing steps:
// 1. Quality trimming & contamination removal
// 2. Genome alignment (mapping)
// 3. PCR duplicate marking
// 4. Base quality score recalibration
//
// 概述：
// 该工作流通过4个主要预处理步骤将原始测序数据转换为可用于变异检测的比对文件：
// 1. 质量修剪和污染物去除
// 2. 基因组比对（mapping）
// 3. PCR重复标记
// 4. 碱基质量分数校正

workflow FASTQ_PREPROCESS_GATK {
    take:
        // INPUT CHANNELS
        // All input data channels needed for preprocessing
        input_fastq                    // Raw sequencing reads (FASTQ format)
        input_sample                   // Sample metadata and BAM/CRAM files
        dict                           // GATK sequence dictionary (defines chromosomes)
        fasta                          // Reference genome sequence
        fasta_fai                      // Reference genome index (for fast access)
        index_alignment                // Pre-built aligner index (BWA/Bowtie2)
        intervals_and_num_intervals    // Genomic regions for parallel processing
        intervals_for_preprocessing    // Genomic intervals for BQSR
        known_sites_indels             // VCF file of known indels (for BQSR)
        known_sites_indels_tbi         // Index of known indels file

        // 输入数据通道
        // 预处理所需的所有输入数据通道
        bbsplit_index                  // BBSplit index for contamination detection

    main:

    // ────────────────────────────────────────────────────────────────────────
    // Initialize channels for collecting quality control reports and software versions
    // All QC reports are gathered for MultiQC visualization
    // 初始化用于收集质量控制报告和软件版本的通道
    // 所有质量控制报告都被收集用于MultiQC可视化
    // ────────────────────────────────────────────────────────────────────────
    reports          = Channel.empty()  // QC reports (FastQC, Picard metrics, mosdepth)
    versions         = Channel.empty()  // Version information from each tool

    // PREPROCESSING
    // ────────────────────────────────────────────────────────────────────────

    if (params.step == 'mapping') {

        // ╔══════════════════════════════════════════════════════════════════╗
        // ║  STAGE 0: QUALITY CONTROL & TRIMMING                           ║
        // ║  第0步: 质量控制和修剪                                         ║
        // ╚══════════════════════════════════════════════════════════════════╝
        //
        // BIOLOGICAL PURPOSE:
        // Raw sequencing reads contain sequencing errors and adapter sequences
        // that must be removed before alignment. This stage:
        // - Removes Illumina adapter sequences
        // - Trims low-quality bases (typically last bases of reads)
        // - Optionally generates UMI consensus sequences for error correction
        // - Removes reads from contaminating organisms
        //
        // 生物学目的：
        // 原始测序reads含有测序错误和接头序列，必须在比对前移除。这一步骤：
        // - 去除Illumina接头序列
        // - 修剪低质量的碱基（通常是reads的末端）
        // - 可选地生成UMI共识序列用于误差校正
        // - 去除来自污染生物体的reads
        // ─────────────────────────────────────────────────────────────────

        // UMI CONSENSUS CALLING (optional)
        // UMI共识调用（可选）
        // When UMI (Unique Molecular Identifier) sequences are available:
        // Groups reads with identical UMIs (same original DNA molecule)
        // and creates a consensus sequence. This dramatically improves accuracy
        // by reducing sequencing errors (which appear randomly across reads from
        // the same molecule, and can be averaged out).
        //
        // 当有UMI（分子条形码）序列时：
        // 将具有相同UMI的reads（来自同一原始DNA分子）分组
        // 并创建共识序列。这通过减少测序错误大大改进准确性
        // （测序错误在同一分子的reads中随机出现，可以被平均消除）
        if (params.umi_read_structure) {
            FASTQ_CREATE_UMI_CONSENSUS_FGBIO(
                input_fastq,
                fasta,
                fasta_fai,
                index_alignment,
                params.group_by_umi_strategy)

            // Convert consensus BAM back to FASTQ for downstream processing
            // 将共识BAM转换回FASTQ用于下游处理
            bam_converted_from_fastq = FASTQ_CREATE_UMI_CONSENSUS_FGBIO.out.consensusbam.map{ meta, bam -> [ meta, bam, [] ] }

            // Convert back to fastq for further preprocessing
            // Note: No need for fasta.fai or fasta when converting BAM to FASTQ
            // 转换回fastq以继续预处理
            // 注意：转换BAM到FASTQ时不需要fasta.fai或fasta
            interleave_input = false // Currently don't allow interleaved input
            CONVERT_FASTQ_UMI(
                bam_converted_from_fastq,
                [ [ id:"fasta" ], [] ], // fasta (not needed, empty)
                [ [ id:'null' ], [] ],  // fasta_fai (not needed, empty)
                interleave_input)

            reads_for_fastp = CONVERT_FASTQ_UMI.out.reads

            // Gather used softwares versions for reporting
            // 收集使用过的软件版本用于报告
            versions = versions.mix(CONVERT_FASTQ_UMI.out.versions)
            versions = versions.mix(FASTQ_CREATE_UMI_CONSENSUS_FGBIO.out.versions)
        } else {
            // If no UMI processing, use raw input reads directly
            // 如果没有UMI处理，直接使用原始输入reads
            reads_for_fastp = input_fastq
        }

        // ────────────────────────────────────────────────────────────────
        // TRIMMING AND/OR SPLITTING
        // FASTP: Remove adapters and trim low-quality bases
        //
        // 修剪和/或分割
        // FASTP：去除接头和修剪低质量碱基
        //
        // Why trim?
        // - Illumina sequencing adds adapter sequences at read ends
        // - Quality scores decrease toward the end of reads
        // - Untrimmed adapters can map to wrong locations, causing false variants
        //
        // 为什么要修剪？
        // - Illumina测序在reads末端添加接头序列
        // - 质量分数向reads末端下降
        // - 未修剪的接头可能映射到错误位置，导致假阳性变异
        // ────────────────────────────────────────────────────────────────
        if (params.trim_fastq || params.split_fastq > 0 || params.umi_location) {

            save_trimmed_fail = false  // Don't save reads that fail trimming
            save_merged = false        // Don't save merged paired-end reads
            FASTP(
                reads_for_fastp,
                [], // we are not using any adapter fastas at the moment
                false, // we don't use discard_trimmed_pass at the moment
                save_trimmed_fail,
                save_merged
            )

            // Collect QC reports (JSON and HTML) for MultiQC visualization
            // 收集质量控制报告（JSON和HTML）用于MultiQC可视化
            reports = reports.mix(FASTP.out.json.collect{ _meta, json -> json })
            reports = reports.mix(FASTP.out.html.collect{ _meta, html -> html })

            // If splitting FASTQ files (for parallel processing of large files)
            // Split large FASTQ files into smaller chunks for faster processing
            // 如果分割FASTQ文件（用于大文件的并行处理）
            // 将大的FASTQ文件分割成更小的块以加快处理速度
            if (params.split_fastq) {
                reads_for_bbsplit = FASTP.out.reads.map{ meta, reads ->
                    // Group split files in pairs (R1, R2)
                    // 将分割文件配对分组（R1，R2）
                    def read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                    [ meta + [ n_fastq: read_files.size() ], read_files ]
                }.transpose()
            } else reads_for_bbsplit = FASTP.out.reads

            versions = versions.mix(FASTP.out.versions)

        } else {
            // If trimming is disabled, pass reads through unchanged
            // 如果禁用修剪，直接传递reads
            reads_for_bbsplit = reads_for_fastp
        }

        // ────────────────────────────────────────────────────────────────
        // REMOVE GENOMIC CONTAMINANTS WITH BBSPLIT
        //
        // 使用BBSplit去除基因组污染物
        //
        // Why contamination detection?
        // Many samples can be contaminated with DNA from other organisms:
        // - Bacteria/fungi during sample collection or DNA extraction
        // - Bacterial contamination is common in food samples, skin samples
        // - Contaminating reads will map to bacterial genomes instead
        //   of the human reference, creating incorrect alignments
        //
        // BBSplit maps reads to multiple genomes simultaneously and removes
        // reads that align better to contaminant genomes than human reference
        //
        // 为什么要检测污染物？
        // 许多样本可能被其他生物体的DNA污染：
        // - 样本收集或DNA提取过程中的细菌/真菌
        // - 细菌污染在食品样本、皮肤样本中很常见
        // - 污染reads会映射到细菌基因组而不是人类参考，导致比对错误
        //
        // BBSplit同时将reads映射到多个基因组，并移除
        // 更好地对齐到污染物基因组而不是人类参考的reads
        // ────────────────────────────────────────────────────────────────
        if (params.tools && params.tools.split(',').contains('bbsplit')) {

            reads_for_alignment = BBMAP_BBSPLIT (
                                        reads_for_bbsplit,
                                        bbsplit_index,
                                        [],
                                        [ [], [] ],
                                        false
                                    )
                                    .primary_fastq  // Keep only reads mapping to primary reference

            // Collect contamination statistics
            // 收集污染物统计信息
            reports = reports.mix(BBMAP_BBSPLIT.out.stats.collect{ _meta, stats -> stats })

        } else {
            // If BBSplit is disabled, use all reads
            // 如果禁用BBSplit，使用所有reads
            reads_for_alignment = reads_for_bbsplit
        }


        // ╔══════════════════════════════════════════════════════════════════╗
        // ║  STAGE 1: MAPPING/ALIGNMENT                                    ║
        // ║  第1步: 读取比对/对齐                                          ║
        // ╚══════════════════════════════════════════════════════════════════╝
        //
        // BIOLOGICAL PURPOSE:
        // Sequencing reads must be aligned to a reference genome to:
        // - Determine their genomic location
        // - Identify genetic variants (SNPs, indels)
        // - Calculate coverage at each position
        // - Enable downstream analysis (variant calling, copy number analysis)
        //
        // This stage uses BWA-MEM (or BWA-MEM2), the most accurate aligner
        // for short reads. BWA-MEM uses a seed-and-extend strategy:
        // 1. Seeds: Find short, exact matches (typically 19bp) between read and reference
        // 2. Extend: Extend those seeds with gap-aware alignment
        // 3. Chain: Link multiple extended regions together
        // 4. Realign: Final alignment refinement
        //
        // 生物学目的：
        // 测序reads必须对齐到参考基因组以：
        // - 确定其基因组位置
        // - 识别遗传变异（SNP、插入缺失）
        // - 计算每个位置的测序深度
        // - 启用下游分析（变异检测、拷贝数分析）
        //
        // 本阶段使用BWA-MEM（或BWA-MEM2），最准确的短reads比对工具
        // BWA-MEM使用seed-and-extend策略：
        // 1. Seeds：在read和参考之间查找短的精确匹配（通常19bp）
        // 2. Extend：用间隙感知对齐扩展这些seeds
        // 3. Chain：链接多个扩展区域
        // 4. Realign：最终对齐细化
        // ─────────────────────────────────────────────────────────────────

        // STEP 1a: Calculate number of FASTQ files per sample
        // Some samples have multiple lanes of sequencing (technical replication)
        // We need to group them back together after mapping
        //
        // 第1a步：计算每个样本的FASTQ文件数
        // 某些样本有多个测序lanes（技术重复）
        // 我们需要在比对后将它们重新分组在一起
        reads_for_alignment.map { meta, reads ->
                [ meta.subMap('patient', 'sample', 'sex', 'status'), reads ]  // Extract only relevant metadata
            }
            .groupTuple()  // Group reads from same sample
            .map { meta, reads ->
                meta + [ n_fastq: reads.size() ]  // Count how many FASTQ files for this sample
            }
            .set { reads_grouping_key }

        // Update metadata for single-lane samples
        // 为单lane样本更新元数据
        reads_for_alignment = reads_for_alignment.map{ meta, reads ->
            // If single lane AND single FASTQ file, mark as complete sample
            // 如果单lane且单FASTQ文件，标记为完整样本
            if (meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads ]
            else [ meta, reads ]
        }

        // ────────────────────────────────────────────────────────────────
        // FASTQ_ALIGN: Map reads to reference genome using BWA-MEM
        //
        // FASTQ_ALIGN：使用BWA-MEM将reads映射到参考基因组
        //
        // BWA-MEM algorithm produces a SAM file containing:
        // - Read sequence and quality scores
        // - Mapping position on reference
        // - Mapping quality (confidence in the alignment)
        // - Detailed alignment information (CIGAR string)
        //
        // BWA-MEM算法生成包含以下内容的SAM文件：
        // - Read序列和质量分数
        // - 参考上的映射位置
        // - 映射质量（对齐的信心）
        // - 详细的对齐信息（CIGAR字符串）
        // ────────────────────────────────────────────────────────────────
        sort_bam = true  // Sort BAM by genomic coordinate for faster downstream access
        FASTQ_ALIGN(reads_for_alignment, index_alignment, sort_bam, fasta, fasta_fai)

        aligned_bam = Channel.empty()
        aligned_bai = Channel.empty()

        // STEP 1b: Copy UMI from read names to RX tag (if applicable)
        //
        // UMI (Unique Molecular Identifiers) are DNA barcodes added during library prep
        // Some protocols encode UMIs in the read name (e.g., "@read_name_UMI123")
        // This step copies that UMI to the standard BAM "RX" tag for consistency
        //
        // 第1b步：将UMI从read名称复制到RX标签（如果适用）
        //
        // UMI（分子条形码）是文库准备期间添加的DNA条形码
        // 某些方案在read名称中编码UMI（例如"@read_name_UMI123"）
        // 此步骤将UMI复制到标准BAM"RX"标签以保证一致性
        if (params.umi_in_read_header || params.umi_location) {
            FGBIO_COPYUMIFROMREADNAME(FASTQ_ALIGN.out.bam.map{meta, bam -> [meta, bam, []]})
            aligned_bam = FGBIO_COPYUMIFROMREADNAME.out.bam
            aligned_bai = FGBIO_COPYUMIFROMREADNAME.out.bai
            versions = versions.mix(FGBIO_COPYUMIFROMREADNAME.out.versions)
        } else {
            aligned_bam = FASTQ_ALIGN.out.bam
            aligned_bai = FASTQ_ALIGN.out.bai
        }

        // ────────────────────────────────────────────────────────────────
        // STEP 1c: Group BAM files from same sample
        //
        // 第1c步：将同一样本的BAM文件分组
        //
        // When a sample was sequenced on multiple lanes, we now have multiple BAM files
        // We need to combine them while tracking which ones belong to the same sample
        // groupKey() ensures that related BAMs are processed together
        //
        // 当样本在多个lanes上测序时，我们现在有多个BAM文件
        // 我们需要将它们组合在一起，同时追踪哪些属于同一样本
        // groupKey()确保相关的BAM一起处理
        // ────────────────────────────────────────────────────────────────
        bam_mapped = aligned_bam
            .combine(reads_grouping_key)  // Attach sample grouping information
            .filter { meta1, _bam, meta2 -> meta1.sample == meta2.sample }  // Keep only matching samples
            // Add n_fastq (number of files) to metadata
            .map { meta1, bam, meta2 ->
                [ meta1 + meta2, bam ]
            }
            // Clean up metadata: remove fields no longer needed, mark as BAM
            .map { meta, bam ->
                [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'read_group', 'size', 'sample_lane_id') + [ data_type: 'bam', id: meta.sample ], bam ]
            }
            // Create groupKey: tells workflow to group files and wait for all related files
            .map { meta, bam ->
                [ groupKey( meta, meta.n_fastq), bam ]
            }
            // Perform the grouping
            .groupTuple()

        // Same grouping for BAM index files
        // 对BAM索引文件进行相同的分组
        bai_mapped = aligned_bai
            .combine(reads_grouping_key)
            .filter { meta1, _bai, meta2 -> meta1.sample == meta2.sample }
            .map { meta1, bai, meta2 ->
                [ meta1 + meta2, bai ]
            }
            .map { meta, bai ->
                [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'read_group', 'size', 'sample_lane_id') + [ data_type: 'bai', id: meta.sample ], bai ]
            }
            .map { meta, bai ->
                [ groupKey( meta, meta.n_fastq), bai ]
            }
            .groupTuple()


        // ────────────────────────────────────────────────────────────────
        // OPTIONAL: Merge BAM files and convert to CRAM
        //
        // 可选：合并BAM文件并转换为CRAM
        //
        // Why merge?
        // - GATK MarkDuplicates can process multiple BAM files together
        // - But if user wants to save mapped data, we merge and convert to space-efficient CRAM format
        // 为什么合并？
        // - GATK MarkDuplicates可以同时处理多个BAM文件
        // - 但如果用户想保存比对数据，我们合并并转换为节省空间的CRAM格式
        // ────────────────────────────────────────────────────────────────
        if (
            params.save_mapped ||
            (
                (params.skip_tools && params.skip_tools.split(',').contains('markduplicates')) &&
                !(params.tools && params.tools.split(',').contains('sentieon_dedup'))
            )
        ) {
            // Merge BAM files from same sample, create index, convert to CRAM
            // 合并同一样本的BAM文件，创建索引，转换为CRAM
            BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)

            // Convert BAM to CRAM (30-50% smaller file size, same information content)
            // 将BAM转换为CRAM（文件大小减少30-50%，信息内容相同）
            BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, fasta, fasta_fai)

            // Create CSV metadata file to allow workflow restart from this step
            // 创建CSV元数据文件以允许从此步骤重新启动工作流
            if (params.save_output_as_bam) CHANNEL_ALIGN_CREATE_CSV(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, params.outdir, params.save_output_as_bam)
            else CHANNEL_ALIGN_CREATE_CSV(BAM_TO_CRAM_MAPPING.out.cram.join(BAM_TO_CRAM_MAPPING.out.crai, failOnDuplicate: true, failOnMismatch: true), params.outdir, params.save_output_as_bam)

            // Gather used softwares versions for reporting
            // 收集使用过的软件版本用于报告
            versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
            versions = versions.mix(BAM_TO_CRAM_MAPPING.out.versions)
        }

        // Gather used softwares versions
        // 收集使用过的软件版本
        versions = versions.mix(FASTQ_ALIGN.out.versions)
    }

    if (params.step in ['mapping', 'markduplicates']) {

        // ╔══════════════════════════════════════════════════════════════════╗
        // ║  STAGE 2: MARK DUPLICATES                                      ║
        // ║  第2步: 标记重复序列                                            ║
        // ╚══════════════════════════════════════════════════════════════════╝
        //
        // BIOLOGICAL PURPOSE:
        // During DNA library preparation, PCR amplification creates exact copies of DNA fragments
        // These PCR duplicates are artifacts - they don't represent unique DNA molecules
        //
        // Why is this a problem?
        // - PCR duplicates can create false "supporting evidence" for variants
        // - Example: If a variant appears in 1 original DNA molecule, but PCR amplified it 100x,
        //   the variant appears to have 100x support (making it seem much more confident)
        // - This can lead to calling false variants
        //
        // What this stage does:
        // Identifies which reads are PCR duplicates (reads with identical start position and sequence)
        // and marks them with the SAM "duplicate" flag. These marked reads are then excluded from
        // variant calling in downstream analysis.
        //
        // 生物学目的：
        // 在DNA文库制备期间，PCR扩增会产生DNA片段的精确副本
        // 这些PCR重复是人工制品 - 它们不代表独特的DNA分子
        //
        // 为什么这是个问题？
        // - PCR重复可能为变异创建假的"支持证据"
        // - 示例：如果变异出现在1个原始DNA分子中，但PCR扩增了100倍，
        //   变异似乎有100倍的支持（使其看起来更有信心）
        // - 这可能导致调用假阳性变异
        //
        // 此阶段的作用：
        // 识别哪些reads是PCR重复（具有相同起始位置和序列的reads）
        // 并用SAM"duplicate"标志标记它们。然后在下游分析中的变异检测中排除这些标记的reads
        // ─────────────────────────────────────────────────────────────────

        cram_markduplicates_no_spark = Channel.empty()
        cram_sentieon_dedup          = Channel.empty()
        cram_markduplicates_spark    = Channel.empty()

        // Determine input source for MarkDuplicates
        // If continuing from 'mapping' step: use BAM from FASTQ_ALIGN
        // If starting from 'markduplicates' step: use input BAM/CRAM
        //
        // 确定MarkDuplicates的输入源
        // 如果从"mapping"步骤继续：使用来自FASTQ_ALIGN的BAM
        // 如果从"markduplicates"步骤开始：使用输入BAM/CRAM
        cram_for_markduplicates = params.step == 'mapping' ? bam_mapped : input_sample.map{ meta, input, _index -> [ meta, input ] }

        // Copy UMI if needed (when restarting from markduplicates step)
        // 如果需要复制UMI（从markduplicates步骤重新启动时）
        if(params.step == 'markduplicates' && params.umi_in_read_header) {
            FGBIO_COPYUMIFROMREADNAME(cram_for_markduplicates.map{ meta, bam -> [ meta, bam, [] ] })
            cram_for_markduplicates = FGBIO_COPYUMIFROMREADNAME.out.bam
            versions = versions.mix(FGBIO_COPYUMIFROMREADNAME.out.versions)
        }

        // ────────────────────────────────────────────────────────────────
        // Three options for marking duplicates:
        //
        // 标记重复的三个选项：
        //
        // 1. BAM_MARKDUPLICATES: Standard GATK MarkDuplicates
        //    - Most accurate, uses coordinate-based duplicate detection
        //    - Slower for very large BAM files (>500GB)
        //
        // 2. BAM_MARKDUPLICATES_SPARK: Spark version (for large datasets)
        //    - Distributed processing using Apache Spark
        //    - Faster for large samples (whole genome, deep sequencing)
        //
        // 3. BAM_SENTIEON_DEDUP: Optimized commercial tool
        //    - Fastest option (if Sentieon license available)
        //    - Uses GPU acceleration
        // ────────────────────────────────────────────────────────────────

        cram_skip_markduplicates = Channel.empty()

        if (
            params.skip_tools &&
            params.skip_tools.split(',').contains('markduplicates') &&
            !(params.tools && params.tools.split(',').contains('sentieon_dedup'))
        ) {
            // SKIP MARKING DUPLICATES (not recommended for variant calling)
            // Only convert BAM to CRAM if needed
            //
            // 跳过标记重复（不推荐用于变异检测）
            // 只在需要时将BAM转换为CRAM
            if (params.step == 'mapping') {
                cram_skip_markduplicates = BAM_TO_CRAM_MAPPING.out.cram.join(BAM_TO_CRAM_MAPPING.out.crai, failOnDuplicate: true, failOnMismatch: true)
            } else {
                cram_skip_markduplicates = Channel.empty().mix(input_sample)
            }

            // Run QC on mapped reads (without duplicate marking)
            // 在映射读取上运行质量控制（不标记重复）
            CRAM_QC_NO_MD(cram_skip_markduplicates, fasta, intervals_for_preprocessing)

            // Gather QC reports
            // 收集质量控制报告
            reports = reports.mix(CRAM_QC_NO_MD.out.reports.collect{ _meta, report -> [ report ] })

            // Gather used softwares versions
            // 收集使用过的软件版本
            versions = versions.mix(CRAM_QC_NO_MD.out.versions)

        } else if (params.use_gatk_spark && params.use_gatk_spark.contains('markduplicates')) {
            // OPTION 2: Use Spark version for large datasets
            // 选项2：对大型数据集使用Spark版本
            BAM_MARKDUPLICATES_SPARK(
                cram_for_markduplicates,
                dict,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)
            cram_markduplicates_spark = BAM_MARKDUPLICATES_SPARK.out.cram

            // Gather QC reports (duplicate metrics from Picard)
            // 收集质量控制报告（来自Picard的重复指标）
            reports = reports.mix(BAM_MARKDUPLICATES_SPARK.out.reports.collect{ _meta, report -> [ report ] })

            // Gather used softwares versions
            // 收集使用过的软件版本
            versions = versions.mix(BAM_MARKDUPLICATES_SPARK.out.versions)

        } else if (params.tools && params.tools.split(',').contains('sentieon_dedup')) {
            // OPTION 3: Use Sentieon optimized tool
            // 选项3：使用Sentieon优化工具
            crai_for_markduplicates = params.step == 'mapping'
                ? bai_mapped
                : ( params.umi_in_read_header ? FGBIO_COPYUMIFROMREADNAME.out.bai : input_sample.map{ meta, _input, index -> [ meta, index ] } )
            BAM_SENTIEON_DEDUP(
                cram_for_markduplicates,
                crai_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            cram_sentieon_dedup = BAM_SENTIEON_DEDUP.out.cram

            // Gather QC reports
            // 收集质量控制报告
            reports = reports.mix(BAM_SENTIEON_DEDUP.out.reports.collect{ _meta, report -> [ report ] })

            // Gather used softwares versions
            // 收集使用过的软件版本
            versions = versions.mix(BAM_SENTIEON_DEDUP.out.versions)

        } else {
            // OPTION 1: Standard GATK MarkDuplicates (default)
            // 选项1：标准GATK MarkDuplicates（默认）
            BAM_MARKDUPLICATES(
                cram_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            cram_markduplicates_no_spark = BAM_MARKDUPLICATES.out.cram

            // Gather QC reports (Picard MarkDuplicates metrics)
            // 收集质量控制报告（Picard MarkDuplicates指标）
            reports = reports.mix(BAM_MARKDUPLICATES.out.reports.collect{ _meta, report -> [ report ] })

            // Gather used softwares versions
            // 收集使用过的软件版本
            versions = versions.mix(BAM_MARKDUPLICATES.out.versions)
        }

        // ────────────────────────────────────────────────────────────────
        // Combine outputs from different MarkDuplicates methods
        //
        // 合并来自不同MarkDuplicates方法的输出
        //
        // The workflow can use different deduplication tools depending on:
        // - Dataset size (Spark for large samples)
        // - Available software licenses (Sentieon if licensed)
        // - User preference (standard GATK as default)
        // We combine all outputs into one channel for downstream processing
        //
        // 工作流可以根据以下条件使用不同的去重工具：
        // - 数据集大小（大样本使用Spark）
        // - 可用的软件许可证（如有许可使用Sentieon）
        // - 用户偏好（标准GATK作为默认）
        // 我们将所有输出合并到一个通道中以进行下游处理
        // ────────────────────────────────────────────────────────────────
        ch_md_cram_for_restart = Channel.empty().mix(cram_markduplicates_no_spark, cram_markduplicates_spark, cram_sentieon_dedup)
            // Mark all files as CRAM format
            // 将所有文件标记为CRAM格式
            .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

        // OPTIONAL: Convert CRAM back to BAM if user prefers BAM format
        // 可选：如果用户更喜欢BAM格式，将CRAM转换回BAM
        CRAM_TO_BAM(ch_md_cram_for_restart, fasta, fasta_fai)
        versions = versions.mix(CRAM_TO_BAM.out.versions)

        // ────────────────────────────────────────────────────────────────
        // Create checkpoint file for workflow restart
        // 为工作流重新启动创建检查点文件
        // ────────────────────────────────────────────────────────────────
        csv_subfolder = (params.tools && params.tools.split(',').contains('sentieon_dedup')) ? 'sentieon_dedup' : 'markduplicates'

        if (params.save_output_as_bam) CHANNEL_MARKDUPLICATES_CREATE_CSV(CRAM_TO_BAM.out.bam.join(CRAM_TO_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true), csv_subfolder, params.outdir, params.save_output_as_bam)
        else CHANNEL_MARKDUPLICATES_CREATE_CSV(ch_md_cram_for_restart, csv_subfolder, params.outdir, params.save_output_as_bam)
    }

    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration']) {

        // ╔══════════════════════════════════════════════════════════════════╗
        // ║  STAGE 3: BASE QUALITY SCORE RECALIBRATION (BQSR)              ║
        // ║  第3步: 碱基质量分数重新校准（BQSR）                           ║
        // ╚══════════════════════════════════════════════════════════════════╝
        //
        // BIOLOGICAL PURPOSE:
        // Sequencing instruments (Illuminators, Ion Proton, PacBio, Nanopore, etc.)
        // produce quality scores (Phred scores) that estimate the probability of error
        // for each base in each read.
        //
        // Problem: These quality scores are often systematically inaccurate
        // - May overestimate quality of reads from low-diversity regions
        // - May underestimate quality near read ends
        // - May depend on dinucleotide context (some DNA sequences are harder to sequence)
        // - May depend on machine and flowcell batch
        //
        // Solution: BQSR (Base Quality Score Recalibration)
        // BQSR compares the reported quality scores to actual error patterns observed in
        // the data by looking at known variants (from dbSNP, ClinVar, 1000 Genomes).
        //
        // The algorithm:
        // 1. For each base position, calculate what the actual error rate should be
        //    (by comparing alignments to known variant sites)
        // 2. Create a lookup table: given certain covariates (context, position, etc.),
        //    what should the quality score be?
        // 3. Apply this table to adjust all quality scores
        //
        // Result: Better-calibrated quality scores lead to:
        // - More accurate variant confidence estimates
        // - Better variant call filtering
        // - Improved accuracy in genotype calling
        //
        // 生物学目的：
        // 测序仪（Illumina、Ion Proton、PacBio、Nanopore等）
        // 产生质量分数（Phred分数），估计每个read中每个碱基的误差概率
        //
        // 问题：这些质量分数通常系统性地不准确
        // - 可能高估低多样性区域reads的质量
        // - 可能低估read末端附近的质量
        // - 可能取决于二核苷酸背景（某些DNA序列更难测序）
        // - 可能取决于机器和流动池批次
        //
        // 解决方案：BQSR（碱基质量分数重新校准）
        // BQSR通过查看已知变异（来自dbSNP、ClinVar、1000基因组）
        // 将报告的质量分数与数据中观察到的实际错误模式进行比较
        //
        // 算法：
        // 1. 对于每个碱基位置，计算实际错误率应该是多少
        //    （通过将对齐与已知变异位点进行比较）
        // 2. 创建查找表：给定某些协变量（背景、位置等），
        //    质量分数应该是什么？
        // 3. 应用此表调整所有质量分数
        //
        // 结果：更好校准的质量分数导致：
        // - 更准确的变异置信度估计
        // - 更好的变异调用过滤
        // - 改进的基因型调用准确性
        // ─────────────────────────────────────────────────────────────────

        // Determine input source for BaseRecalibrator
        // If starting from 'prepare_recalibration': use input BAM/CRAM directly
        // Otherwise: use output from MarkDuplicates step
        //
        // 确定BaseRecalibrator的输入源
        // 如果从"prepare_recalibration"开始：直接使用输入BAM/CRAM
        // 否则：使用来自MarkDuplicates步骤的输出
        if (params.step == 'prepare_recalibration') {

            ch_cram_for_bam_baserecalibrator = Channel.empty().mix(input_sample)

            // Set the input samples for restart so we generate a samplesheet that contains the input files together with the recalibration table
            // 设置重新启动的输入样本，以便我们生成一个包含输入文件以及校准表的样本列表
            ch_md_cram_for_restart           = ch_cram_for_bam_baserecalibrator

        } else {

            // ch_cram_for_bam_baserecalibrator contains either:
            // - CRAMs from MarkDuplicates (if step == 'mapping')
            // - CRAMs from MarkDuplicates_Spark
            // - CRAMs converted from BAM (if skipping MarkDuplicates)
            // - Input CRAM files (if step == 'markduplicates')
            //
            // ch_cram_for_bam_baserecalibrator包含：
            // - 来自MarkDuplicates的CRAM（如果step == 'mapping'）
            // - 来自MarkDuplicates_Spark的CRAM
            // - 从BAM转换的CRAM（如果跳过MarkDuplicates）
            // - 输入CRAM文件（如果step == 'markduplicates'）
            ch_cram_for_bam_baserecalibrator = Channel.empty().mix(ch_md_cram_for_restart, cram_skip_markduplicates )
                // Ensure correct data type marking
                // 确保正确的数据类型标记
                .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

        }

        // ────────────────────────────────────────────────────────────────
        // Run BaseRecalibrator if not explicitly skipped
        //
        // 如果未明确跳过，运行BaseRecalibrator
        //
        // This step builds the recalibration model:
        // 1. Reads the input BAM/CRAM file
        // 2. Looks at known variant positions (from known_sites_indels VCF)
        // 3. At each known variant position, calculates:
        //    "Does the reported quality score match the actual error rate?"
        // 4. Builds a statistical table of corrections needed
        // 5. Outputs a .recal_table binary file (used in next step)
        //
        // 此步骤构建校准模型：
        // 1. 读取输入BAM/CRAM文件
        // 2. 查看已知变异位置（来自known_sites_indels VCF）
        // 3. 在每个已知变异位置，计算：
        //    "报告的质量分数是否与实际错误率相匹配？"
        // 4. 构建所需更正的统计表
        // 5. 输出.recal_table二进制文件（用于下一步）
        // ────────────────────────────────────────────────────────────────
        if (!(params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'))) {

            ch_table_bqsr_no_spark = Channel.empty()
            ch_table_bqsr_spark    = Channel.empty()

            if (params.use_gatk_spark && params.use_gatk_spark.contains('baserecalibrator')) {
                // Use Spark version for large datasets
                // 对大型数据集使用Spark版本
                BAM_BASERECALIBRATOR_SPARK(
                    ch_cram_for_bam_baserecalibrator,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals_and_num_intervals,
                    known_sites_indels,
                    known_sites_indels_tbi)

                ch_table_bqsr_spark = BAM_BASERECALIBRATOR_SPARK.out.table_bqsr

                // Gather used softwares versions
                // 收集使用过的软件版本
                versions = versions.mix(BAM_BASERECALIBRATOR_SPARK.out.versions)
            } else {
                // Use standard GATK BaseRecalibrator (default)
                // 使用标准GATK BaseRecalibrator（默认）
                BAM_BASERECALIBRATOR(
                    ch_cram_for_bam_baserecalibrator,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals_and_num_intervals,
                    known_sites_indels,
                    known_sites_indels_tbi)

                ch_table_bqsr_no_spark = BAM_BASERECALIBRATOR.out.table_bqsr

                // Gather used softwares versions
                // 收集使用过的软件版本
                versions = versions.mix(BAM_BASERECALIBRATOR.out.versions)
            }

            // Combine recalibration tables from both methods (Spark or standard)
            // 合并来自两种方法（Spark或标准）的校准表
            ch_table_bqsr = Channel.empty().mix(
                ch_table_bqsr_no_spark,
                ch_table_bqsr_spark)

            // Collect recalibration tables for QC/reporting
            // 收集校准表用于质量控制/报告
            reports = reports.mix(ch_table_bqsr.collect{ _meta, table -> [ table ] })

            // Join CRAM files with their corresponding recalibration tables
            // (needed for ApplyBQSR step)
            // 将CRAM文件与其对应的校准表连接
            // （ApplyBQSR步骤需要）
            cram_applybqsr = ch_cram_for_bam_baserecalibrator.join(ch_table_bqsr, failOnDuplicate: true, failOnMismatch: true)

            // Create checkpoint file for workflow restart
            // 为工作流重新启动创建检查点文件
            CHANNEL_BASERECALIBRATOR_CREATE_CSV(ch_md_cram_for_restart.join(ch_table_bqsr, failOnDuplicate: true), params.tools, params.skip_tools, params.outdir, params.save_output_as_bam)
        }
    }

    // ╔══════════════════════════════════════════════════════════════════╗
    // ║  STAGE 4: APPLY BASE QUALITY SCORE RECALIBRATION               ║
    // ║  第4步: 应用碱基质量分数重新校准                               ║
    // ╚══════════════════════════════════════════════════════════════════╝
    //
    // BIOLOGICAL PURPOSE:
    // This is the final step of preprocessing. Using the recalibration model
    // built in Stage 3, we now:
    // 1. Read the BAM/CRAM file again
    // 2. For each base, look up its corrected quality score in the recalibration table
    // 3. Rewrite the BAM/CRAM file with the updated quality scores
    // 4. Output analysis-ready CRAM (or BAM) files
    //
    // Why is this important?
    // - All downstream analysis (variant calling, filtering, genotyping) relies on quality scores
    // - Accurate quality scores = accurate variant calls
    // - Poor quality score calibration = false positives and false negatives in variant calls
    //
    // The output files are now ready for:
    // - GATK HaplotypeCaller (for SNP/indel calling)
    // - Copy number variation analysis
    // - Structural variant detection
    // - Any other downstream analysis
    //
    // 生物学目的：
    // 这是预处理的最后一步。使用在第3步中构建的校准模型，我们现在：
    // 1. 再次读取BAM/CRAM文件
    // 2. 对于每个碱基，在校准表中查找其校正后的质量分数
    // 3. 用更新的质量分数重写BAM/CRAM文件
    // 4. 输出分析就绪的CRAM（或BAM）文件
    //
    // 为什么这很重要？
    // - 所有下游分析（变异检测、过滤、基因型检测）都依赖于质量分数
    // - 准确的质量分数 = 准确的变异检测
    // - 质量分数校准不良 = 变异检测中的假阳性和假阴性
    //
    // 输出文件现在可用于：
    // - GATK HaplotypeCaller（用于SNP/indel检测）
    // - 拷贝数变异分析
    // - 结构变异检测
    // - 任何其他下游分析
    // ─────────────────────────────────────────────────────────────────

    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate']) {

        // If starting directly from 'recalibrate' step, use input directly
        // 如果直接从"recalibrate"步骤开始，直接使用输入
        if (params.step == 'recalibrate') {

            cram_applybqsr = Channel.empty().mix(input_sample)

        }

        // ────────────────────────────────────────────────────────────────
        // Apply recalibration if BaseRecalibrator was run (not skipped)
        //
        // 如果BaseRecalibrator已运行（未跳过），应用校准
        //
        // This reads the .recal_table file created by BaseRecalibrator and applies
        // the quality score adjustments to the BAM/CRAM file.
        //
        // 这读取BaseRecalibrator创建的.recal_table文件，并将
        // 质量分数调整应用到BAM/CRAM文件
        // ────────────────────────────────────────────────────────────────
        if (!(params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'))) {
            cram_variant_calling_no_spark = Channel.empty()
            cram_variant_calling_spark    = Channel.empty()

            if (params.use_gatk_spark && params.use_gatk_spark.contains('baserecalibrator')) {
                // Use Spark version for large datasets
                // 对大型数据集使用Spark版本
                BAM_APPLYBQSR_SPARK(
                    cram_applybqsr,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals_and_num_intervals)

                cram_variant_calling_spark = BAM_APPLYBQSR_SPARK.out.cram

                // Gather used softwares versions
                // 收集使用过的软件版本
                versions = versions.mix(BAM_APPLYBQSR_SPARK.out.versions)

            } else {
                // Use standard GATK ApplyBQSR (default)
                // 使用标准GATK ApplyBQSR（默认）
                BAM_APPLYBQSR(
                    cram_applybqsr,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals_and_num_intervals)

                cram_variant_calling_no_spark = BAM_APPLYBQSR.out.cram

                // Gather used softwares versions
                // 收集使用过的软件版本
                versions = versions.mix(BAM_APPLYBQSR.out.versions)
            }

            // Combine outputs from both processing methods
            // 合并来自两种处理方法的输出
            cram_variant_calling = Channel.empty().mix(
                cram_variant_calling_no_spark,
                cram_variant_calling_spark)

            // ────────────────────────────────────────────────────────────
            // OPTIONAL: Convert CRAM back to BAM if user prefers
            //
            // 可选：如果用户更喜欢，将CRAM转换回BAM
            //
            // Why CRAM?
            // - CRAM files are 30-50% smaller than BAM files
            // - They store only the differences from the reference
            // - Still contain all quality information needed for analysis
            // - Standard format in major genomic databases
            //
            // Some users still prefer BAM for compatibility with older tools
            //
            // 为什么使用CRAM？
            // - CRAM文件比BAM文件小30-50%
            // - 它们只存储与参考的差异
            // - 仍然包含分析所需的所有质量信息
            // - 主要基因组数据库中的标准格式
            //
            // 某些用户仍然更喜欢BAM以与旧工具兼容
            // ────────────────────────────────────────────────────────────
            CRAM_TO_BAM_RECAL(cram_variant_calling, fasta, fasta_fai)
            versions = versions.mix(CRAM_TO_BAM_RECAL.out.versions)

            // Choose output format based on user preference
            // 根据用户偏好选择输出格式
            csv_recalibration = Channel.empty()
            csv_recalibration = params.save_output_as_bam ? CRAM_TO_BAM_RECAL.out.bam.join(CRAM_TO_BAM_RECAL.out.bai, failOnDuplicate: true, failOnMismatch: true) : cram_variant_calling

            // Create checkpoint file for workflow restart
            // 为工作流重新启动创建检查点文件
            CHANNEL_APPLYBQSR_CREATE_CSV(csv_recalibration, params.outdir, params.save_output_as_bam)

        } else if (params.step == 'recalibrate') {
            // If starting from 'recalibrate' with BQSR skipped:
            // Convert input BAMs/CRAMs to CRAMs if needed
            // 如果从"recalibrate"开始且跳过BQSR：
            // 根据需要将输入BAM/CRAM转换为CRAM
            cram_variant_calling = Channel.empty().mix(input_sample.map{ meta, cram, crai, _table -> [ meta, cram, crai ] })
        } else {
            // If BQSR is skipped but not at 'recalibrate' step:
            // Use CRAMs from duplicate marking step
            // 如果BQSR被跳过但不在"recalibrate"步骤：
            // 使用来自重复标记步骤的CRAM
            cram_variant_calling = Channel.empty().mix(ch_cram_for_bam_baserecalibrator)
        }
    }

    // ────────────────────────────────────────────────────────────────────────
    // WORKFLOW OUTPUTS
    //
    // 工作流输出
    //
    // The final output of this preprocessing workflow is:
    // 1. cram_variant_calling: Analysis-ready BAM/CRAM files with:
    //    - Reads aligned to reference genome
    //    - PCR duplicates marked
    //    - Quality scores calibrated
    //    - Ready for variant calling in next step
    //
    // 2. reports: Quality control metrics for all preprocessing steps
    //    (used by MultiQC for visualization)
    //
    // 3. versions: Software version information for documentation
    //
    // 此预处理工作流的最终输出是：
    // 1. cram_variant_calling：分析就绪的BAM/CRAM文件，包括：
    //    - 比对到参考基因组的reads
    //    - 标记的PCR重复
    //    - 校准的质量分数
    //    - 为下一步变异检测做好准备
    //
    // 2. reports：所有预处理步骤的质量控制指标
    //    （由MultiQC用于可视化）
    //
    // 3. versions：软件版本信息用于文档
    // ────────────────────────────────────────────────────────────────────────

    emit:
    cram_variant_calling  // Analysis-ready alignment files for variant calling
    reports               // QC reports for MultiQC
    versions              // Software version information

}
