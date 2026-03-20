/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ╔══════════════════════════════════════════════════════════════════════════════════╗
// ║  SECTION 1: IMPORT TOOLS & MODULES                                              ║
// ║  导入工具和模块                                                                  ║
// ╚══════════════════════════════════════════════════════════════════════════════════╝
//
// Think of "include" like importing tools from a toolbox. The Sarek pipeline uses
// many specialized bioinformatics tools. Each tool is organized into modules
// (individual tools) or subworkflows (collections of related steps).
//
// 把"include"想象成从工具箱导入工具。Sarek流程使用许多专业的生物信息学工具。
// 每个工具被组织成模块（单个工具）或子流程（相关步骤的集合）。

include { paramsSummaryMap                                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                              } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from 'plugin/nf-core-utils'
include { methodsDescriptionText                            } from '../../subworkflows/local/utils_nfcore_sarek_pipeline'

// ┌─ Input/Output Management ─────────────────────────────────────────────────────┐
// │ 输入/输出管理                                                                    │
// └────────────────────────────────────────────────────────────────────────────────┘
// Create samplesheets to restart from different steps
// 创建样本信息表以从不同的步骤重新开始
include { CHANNEL_VARIANT_CALLING_CREATE_CSV                } from '../../subworkflows/local/channel_variant_calling_create_csv'

// ┌─ Input Format Conversion ─────────────────────────────────────────────────────┐
// │ 输入格式转换                                                                    │
// │ Sometimes input comes as BAM files (aligned reads) instead of raw FASTQ files, │
// │ so we convert them back. Some labs also use SPRING compression for storage.   │
// │                                                                                 │
// │ 有时输入以BAM文件（已比对的序列）的形式而不是原始FASTQ文件出现，                │
// │ 所以我们需要将它们转换回来。一些实验室还使用SPRING压缩进行存储。                 │
// └────────────────────────────────────────────────────────────────────────────────┘
// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT       } from '../../subworkflows/local/bam_convert_samtools'

// Convert fastq.gz.spring files to fastq.gz files
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R1_FQ   } from '../../modules/nf-core/spring/decompress'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R2_FQ   } from '../../modules/nf-core/spring/decompress'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_FQ_PAIR } from '../../modules/nf-core/spring/decompress'

// ┌─ Quality Control (QC) Tools ──────────────────────────────────────────────────┐
// │ 质量控制（QC）工具                                                               │
// │ Before doing expensive processing, we check if the DNA/RNA sequencing data     │
// │ is high quality. This helps catch problems early.                             │
// │                                                                                 │
// │ 在进行昂贵的处理之前，我们检查DNA/RNA测序数据的质量。这有助于提前发现问题。       │
// └────────────────────────────────────────────────────────────────────────────────┘
// Run FASTQC - checks read quality before processing
include { FASTQC                                            } from '../../modules/nf-core/fastqc'

// Quality check on CRAM files (compressed aligned reads)
// CRAM是BAM的压缩版本，用于存储已比对的序列
include { CRAM_SAMPLEQC                                     } from '../../subworkflows/local/cram_sampleqc'

// ┌─ Preprocessing/Alignment ─────────────────────────────────────────────────────┐
// │ 预处理/比对                                                                     │
// │ Raw FASTQ files are compared to the reference genome. This is like matching    │
// │ puzzle pieces - we find where each DNA sequence belongs in the genome.         │
// │                                                                                 │
// │ Two pathways available:                                                        │
// │ 1. GATK (standard, CPU-based) - reliable but slower                           │
// │ 2. Parabricks (GPU-accelerated) - much faster if you have a GPU              │
// │                                                                                 │
// │ 原始FASTQ文件与参考基因组进行比较。这就像匹配拼图——我们找到每个DNA序列在        │
// │ 基因组中的位置。                                                               │
// │                                                                                 │
// │ 两种路径可用：                                                                  │
// │ 1. GATK（标准，基于CPU）-可靠但较慢                                             │
// │ 2. Parabricks（GPU加速）-如果有GPU会快得多                                     │
// └────────────────────────────────────────────────────────────────────────────────┘
include { FASTQ_PREPROCESS_GATK                             } from '../../subworkflows/local/fastq_preprocess_gatk'
include { FASTQ_PREPROCESS_PARABRICKS                       } from '../../subworkflows/local/fastq_preprocess_parabricks'

// CRAM to BAM conversion - sometimes we need the uncompressed format
// CRAM到BAM转换 - 有时我们需要未压缩的格式
include { SAMTOOLS_CONVERT as CRAM_TO_BAM                   } from '../../modules/nf-core/samtools/convert'

// ┌─ Variant Calling Tools ───────────────────────────────────────────────────────┐
// │ 变异检测工具                                                                     │
// │ These tools identify genetic differences (mutations) by comparing a patient's  │
// │ DNA to the reference genome. Different tools are used for different scenarios: │
// │                                                                                 │
// │ GERMLINE (inherited): Changes in normal cells inherited from parents          │
// │ TUMOR-ONLY: Cancer mutations when only tumor sample available               │
// │ SOMATIC PAIRED: Cancer mutations detected by comparing tumor vs normal       │
// │                                                                                 │
// │ 这些工具通过将患者DNA与参考基因组进行比较来识别遗传差异（突变）。                │
// │ 不同的工具用于不同的场景：                                                     │
// │                                                                                 │
// │ 生殖系（遗传）：从父母那里继承的正常细胞的变化                                   │
// │ 肿瘤单独：仅有肿瘤样本可用时的癌症突变                                           │
// │ 体细胞配对：通过比较肿瘤与正常来检测的癌症突变                                   │
// └────────────────────────────────────────────────────────────────────────────────┘
// Variant calling on a single normal sample (for germline mutations)
include { BAM_VARIANT_CALLING_GERMLINE_ALL                  } from '../../subworkflows/local/bam_variant_calling_germline_all'

// Variant calling on a single tumor sample (no matched normal)
include { BAM_VARIANT_CALLING_TUMOR_ONLY_ALL                } from '../../subworkflows/local/bam_variant_calling_tumor_only_all'

// Variant calling on tumor/normal pair (best case - can detect somatic mutations)
include { BAM_VARIANT_CALLING_SOMATIC_ALL                   } from '../../subworkflows/local/bam_variant_calling_somatic_all'

// ┌─ Variant Calling Post-Processing ─────────────────────────────────────────────┐
// │ 变异检测后期处理                                                                 │
// │ After finding variants, we need to merge, filter, and normalize them for      │
// │ analysis and annotation.                                                       │
// │                                                                                 │
// │ 发现变异后，我们需要合并、过滤和标准化它们以进行分析和注释。                     │
// └────────────────────────────────────────────────────────────────────────────────┘
include { POST_VARIANTCALLING                               } from '../../subworkflows/local/post_variantcalling'

// Quality control on variant calls (VCF files)
// 变异调用的质量控制（VCF文件）
include { VCF_QC_BCFTOOLS_VCFTOOLS                          } from '../../subworkflows/local/vcf_qc_bcftools_vcftools'

// ┌─ Annotation Tools ────────────────────────────────────────────────────────────┐
// │ 注释工具                                                                        │
// │ After finding variants, we add biological information (annotations) to each   │
// │ variant. This tells us: What gene is affected? Does this mutation have known  │
// │ effects? Is it likely to be disease-causing?                                 │
// │                                                                                 │
// │ 找到变异后，我们为每个变异添加生物学信息（注释）。这告诉我们：哪个基因受到      │
// │ 影响？这个突变是否有已知的影响？它是否可能导致疾病？                            │
// └────────────────────────────────────────────────────────────────────────────────┘
include { VCF_ANNOTATE_ALL                                  } from '../../subworkflows/local/vcf_annotate_all'

// ┌─ Quality Report Aggregation ──────────────────────────────────────────────────┐
// │ 质量报告汇总                                                                     │
// │ Combines all QC metrics into one comprehensive HTML report for easy viewing    │
// │                                                                                 │
// │ 将所有QC指标合并为一个综合HTML报告，便于查看                                    │
// └────────────────────────────────────────────────────────────────────────────────┘
include { MULTIQC                                           } from '../../modules/nf-core/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAREK {
    take:
    // ╔══════════════════════════════════════════════════════════════════════════════╗
    // ║  INPUT PARAMETERS (inputs from the workflow)                                 ║
    // ║  输入参数（来自工作流的输入）                                                 ║
    // ╚══════════════════════════════════════════════════════════════════════════════╝

    input_sample
    aligner
    skip_tools
    step
    tools

    // ┌─ Reference Database Files ────────────────────────────────────────────────────┐
    // │ 参考数据库文件                                                                 │
    // │                                                                                │
    // │ fasta: The human reference genome - like a standard blueprint that we compare│
    // │        each patient's DNA against. Usually GRCh37 or GRCh38 (human builds)   │
    // │                                                                                │
    // │ dbsnp: Database of known genetic variants. If a variant is already known     │
    // │        and not disease-causing, this database contains it. Helps filter out  │
    // │        common, benign variants.                                             │
    // │                                                                                │
    // │ known_sites_indels/snps: Used for calibrating base quality scores (BQSR).   │
    // │        This improves accuracy of variant calling by using known good variants.│
    // │                                                                                │
    // │ germline_resource: Database of common germline (inherited) variants          │
    // │        Used to filter out non-somatic variants in tumor analysis             │
    // │                                                                                │
    // │ pon (Panel of Normals): Variants found in many normal samples from the lab  │
    // │        These are likely sequencing artifacts, not real mutations             │
    // │                                                                                │
    // │ fasta: 人类参考基因组 - 就像一个标准蓝图，我们将每个患者的DNA与之进行比较。    │
    // │        通常是GRCh37或GRCh38（人类版本）                                      │
    // │                                                                                │
    // │ dbsnp: 已知遗传变异数据库。如果一个变异已为人所知且不引起疾病，                │
    // │        这个数据库包含它。有助于过滤掉常见、良性的变异。                       │
    // │                                                                                │
    // │ known_sites_indels/snps: 用于校准碱基质量分数（BQSR）。这通过使用已知        │
    // │        的良好变异来提高变异检测的准确性。                                     │
    // │                                                                                │
    // │ germline_resource: 常见生殖系（遗传）变异数据库。用于过滤肿瘤分析中的        │
    // │        非体细胞变异                                                           │
    // │                                                                                │
    // │ pon (法线板): 在实验室中许多正常样本中发现的变异。这些可能是测序伪影，      │
    // │        不是真实的突变                                                         │
    // └────────────────────────────────────────────────────────────────────────────────┘
    ascat_alleles
    ascat_loci
    ascat_loci_gc
    ascat_loci_rt
    bbsplit_index
    bcftools_annotations
    bcftools_annotations_tbi
    bcftools_columns
    bcftools_header_lines
    cf_chrom_len
    chr_files
    cnvkit_reference
    dbsnp
    dbsnp_tbi
    dbsnp_vqsr
    dict
    fasta
    fasta_fai
    germline_resource
    germline_resource_tbi
    index_alignment
    // ┌─ Interval Files (Genome Division) ────────────────────────────────────────────┐
    // │ 区间文件（基因组分割）                                                         │
    // │                                                                                │
    // │ The human genome is huge (~3 billion base pairs). To speed up processing,    │
    // │ we divide it into smaller chunks called intervals. Different intervals can   │
    // │ be processed in parallel on different computers - like dividing homework     │
    // │ among students instead of one person doing it all.                          │
    // │                                                                                │
    // │ intervals_bed_combined: All intervals merged into one file for operations   │
    // │ that need the complete set.                                                 │
    // │                                                                                │
    // │ intervals_for_preprocessing: Specific intervals used during alignment and    │
    // │ base quality score recalibration.                                           │
    // │                                                                                │
    // │ 人类基因组是巨大的（约30亿个碱基对）。为了加快处理速度，我们将其分成更小的    │
    // │ 块，称为间隔。不同的间隔可以在不同的计算机上并行处理 - 就像将作业分配给学生   │
    // │ 而不是一个人完成所有工作。                                                     │
    // │                                                                                │
    // │ intervals_bed_combined: 所有间隔合并到一个文件中，供需要完整集合的操作使用。 │
    // │                                                                                │
    // │ intervals_for_preprocessing: 在比对和碱基质量分数重新校准期间使用的特定间隔。 │
    // └────────────────────────────────────────────────────────────────────────────────┘
    intervals_and_num_intervals
    intervals_bed_combined
    intervals_bed_combined_for_variant_calling
    intervals_bed_gz_tbi_and_num_intervals
    intervals_bed_gz_tbi_combined
    intervals_for_preprocessing
    known_indels_vqsr
    known_sites_indels
    known_sites_indels_tbi
    known_sites_snps
    known_sites_snps_tbi
    known_snps_vqsr
    mappability
    msisensor2_models
    msisensorpro_scan
    ngscheckmate_bed
    pon
    pon_tbi
    sentieon_dnascope_model
    varlociraptor_scenario_germline
    varlociraptor_scenario_somatic
    varlociraptor_scenario_tumor_only
    snpeff_cache
    snpeff_db
    vep_cache
    vep_cache_version
    vep_extra_files
    vep_fasta
    vep_genome
    vep_species
    snpsift_db                  // channel: [[databases], [tbis], [vardbs], [fields], [prefixes]]
    versions

    main:
    // To gather all QC reports for MultiQC
    ch_multiqc_files = channel.empty()
    multiqc_publish = channel.empty()
    multiqc_report = channel.empty()
    reports = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 1: INPUT VALIDATION & QUALITY CONTROL
        步骤1：输入验证和质量控制
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // ╔══════════════════════════════════════════════════════════════════════════════╗
    // ║  STEP 1A: INPUT FORMAT DETECTION & CONVERSION                               ║
    // ║  步骤1A：输入格式检测和转换                                                   ║
    // ╚══════════════════════════════════════════════════════════════════════════════╝
    //
    // The pipeline can accept inputs in three formats:
    // 1. FASTQ (.fastq.gz) - raw sequencing data, exactly as it comes from the sequencer
    // 2. BAM - already aligned reads (DNA compared to reference genome)
    // 3. SPRING - compressed FASTQ format used by some laboratories
    //
    // We need to detect which format and convert to FASTQ if necessary for processing.
    //
    // 该流程可以接受三种格式的输入：
    // 1. FASTQ (.fastq.gz) - 原始测序数据，正如它从测序仪中出现的那样
    // 2. BAM - 已比对的序列（DNA与参考基因组比较）
    // 3. SPRING - 一些实验室使用的压缩FASTQ格式
    //
    // 我们需要检测格式并在必要时转换为FASTQ以进行处理。

    // PREPROCESSING
    if (step == 'mapping') {
        // Figure out if input is bam, fastq, or spring
        // 确定输入是bam、fastq还是spring
        input_sample_type = input_sample.branch {
            bam: it[0].data_type == "bam"
            fastq_gz: it[0].data_type == "fastq_gz"
            one_fastq_gz_spring: it[0].data_type == "one_fastq_gz_spring"
            two_fastq_gz_spring: it[0].data_type == "two_fastq_gz_spring"
        }

        // Two fastq.gz-files
        // Standard paired-end sequencing: forward (R1) and reverse (R2) reads
        // 标准配对末端测序：前向（R1）和反向（R2）序列
        fastq_gz = input_sample_type.fastq_gz.map { meta, files -> addReadgroupToMeta(meta, files) }

        // Just one fastq.gz.spring-file with both R1 and R2
        // 一个fastq.gz.spring文件包含R1和R2
        // SPRING decompression extracts the paired reads from one compressed file
        fastq_gz_pair_from_spring = SPRING_DECOMPRESS_TO_FQ_PAIR(input_sample_type.one_fastq_gz_spring, false)

        one_fastq_gz_from_spring = fastq_gz_pair_from_spring.fastq.map { meta, files -> addReadgroupToMeta(meta, files) }

        // Two fastq.gz.spring-files - one for R1 and one for R2
        // 两个fastq.gz.spring文件 - 一个用于R1，一个用于R2
        r1_fastq_gz_from_spring = SPRING_DECOMPRESS_TO_R1_FQ(
            input_sample_type.two_fastq_gz_spring.map { meta, files ->
                [meta, files[0]]
            },
            true,
        )
        r2_fastq_gz_from_spring = SPRING_DECOMPRESS_TO_R2_FQ(
            input_sample_type.two_fastq_gz_spring.map { meta, files ->
                [meta, files[1]]
            },
            true,
        )

        versions = versions.mix(SPRING_DECOMPRESS_TO_R1_FQ.out.versions)
        versions = versions.mix(SPRING_DECOMPRESS_TO_R2_FQ.out.versions)
        versions = versions.mix(SPRING_DECOMPRESS_TO_FQ_PAIR.out.versions)

        two_fastq_gz_from_spring = r1_fastq_gz_from_spring.fastq.join(r2_fastq_gz_from_spring.fastq).map { meta, fastq_1, fastq_2 -> [meta, [fastq_1, fastq_2]] }

        two_fastq_gz_from_spring = two_fastq_gz_from_spring.map { meta, files -> addReadgroupToMeta(meta, files) }

        // ┌─ BAM to FASTQ Conversion ─────────────────────────────────────────────────┐
        // │ BAM转FASTQ转换                                                             │
        // │                                                                            │
        // │ If user provides already-aligned BAM files, we must convert them back to  │
        // │ FASTQ format because we want to realign them (potentially with updated   │
        // │ reference or different alignment parameters). This is like unscrambling   │
        // │ an egg - we extract the original sequences.                              │
        // │                                                                            │
        // │ 如果用户提供已比对的BAM文件，我们必须将它们转换回FASTQ格式，因为我们想   │
        // │ 重新比对它们（可能使用更新的参考或不同的比对参数）。这就像解开一个蛋   │
        // │ - 我们提取原始序列。                                                      │
        // └────────────────────────────────────────────────────────────────────────────┘
        // Convert any bam input to fastq
        // fasta are not needed when converting bam to fastq -> [ id:"fasta" ], []
        // No need for fasta.fai -> []
        // Currently don't allow interleaved input
        interleave_input = false
        CONVERT_FASTQ_INPUT(
            input_sample_type.bam,
            [[id: "fasta"], []],
            [[id: 'null'], []],
            interleave_input,
        )

        versions = versions.mix(CONVERT_FASTQ_INPUT.out.versions)

        // Gather fastq (inputed or converted)
        // 收集fastq（输入或转换）
        // Theorically this could work on mixed input (fastq for one sample and bam for another)
        // But not sure how to handle that with the samplesheet
        // Or if we really want users to be able to do that
        input_fastq = fastq_gz.mix(CONVERT_FASTQ_INPUT.out.reads).mix(one_fastq_gz_from_spring).mix(two_fastq_gz_from_spring)

        // ╔══════════════════════════════════════════════════════════════════════════════╗
        // ║  STEP 1B: FASTQC QUALITY CONTROL                                           ║
        // ║  步骤1B：FASTQC质量控制                                                      ║
        // ╚══════════════════════════════════════════════════════════════════════════════╝
        //
        // FASTQC analyzes raw sequencing reads to check:
        // - Are reads long enough? (usually 50-150 bp)
        // - Is the quality score distribution reasonable?
        // - Are there adapter sequences left over? (should be removed)
        // - Are there GC content biases? (unusual high/low GC can indicate problems)
        // - Are there quality score dropoffs at the ends? (common in old sequencers)
        //
        // Think of it like inspecting a package before processing - if quality is bad,
        // we might want to resequence instead of wasting time analyzing bad data.
        //
        // FASTQC分析原始测序序列以检查：
        // - 序列是否足够长？（通常为50-150 bp）
        // - 质量分数分布是否合理？
        // - 是否有遗留的接头序列？（应该被移除）
        // - 是否有GC含量偏差？（不寻常的高/低GC可能表明问题）
        // - 末尾是否有质量分数下降？（在旧的测序仪中很常见）
        //
        // 可以将其视为在处理前检查包装 - 如果质量差，我们可能想重新测序，
        // 而不是浪费时间分析差的数据。

        // QC
        // `--skip_tools fastqc` to skip fastqc
        if (!(skip_tools.split(',').contains('fastqc'))) {
            FASTQC(input_fastq)

            reports = reports.mix(FASTQC.out.zip.collect { _meta, logs -> logs })
            versions = versions.mix(FASTQC.out.versions)
        }
    }
    else {
        input_fastq = channel.empty().mix(input_sample)
    }

    // ╔══════════════════════════════════════════════════════════════════════════════╗
    // ║  STEP 2: PREPROCESSING (ALIGNMENT & BASE QUALITY SCORE RECALIBRATION)       ║
    // ║  步骤2：预处理（比对和碱基质量分数重新校准）                                 ║
    // ╚══════════════════════════════════════════════════════════════════════════════╝
    //
    // This is the most computationally intensive part of the pipeline. We take raw
    // DNA sequences and align them to a reference genome, then improve the quality
    // scores using machine learning.
    //
    // Key substeps:
    // 1. Trim adapters - remove short plastic sequences added during library prep
    // 2. Align - find where each read belongs in the genome using BWA or other aligners
    // 3. Mark duplicates - PCR amplification creates duplicate copies (same start/end).
    //                     These must be marked because they inflate variant frequencies
    // 4. Base Quality Score Recalibration (BQSR) - improve confidence scores using
    //    known variants. If a position usually has a T but a read has G, either the
    //    read is mutated or the quality score was wrong. BQSR learns from this.
    //
    // Output: CRAM or BAM file - a compressed representation of aligned reads
    //
    // 这是流程中计算量最大的部分。我们取原始DNA序列并将其比对到参考基因组，
    // 然后使用机器学习来改进质量分数。
    //
    // 关键子步骤：
    // 1. 修剪接头 - 移除在文库制备期间添加的短塑料序列
    // 2. 比对 - 使用BWA或其他比对器找到每个序列在基因组中的位置
    // 3. 标记重复序列 - PCR扩增会创建重复副本（相同的开始/结束）。
    //               这些必须被标记，因为它们会增加变异频率
    // 4. 碱基质量分数重新校准（BQSR）- 使用已知变异来改进置信度分数。
    //    如果一个位置通常有T，但一个序列有G，那么要么该序列被突变，
    //    要么质量分数是错误的。BQSR从这个学习。
    //
    // 输出：CRAM或BAM文件 - 已比对序列的压缩表示

    if (step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate']) {

        if (aligner == 'parabricks') {
            // ┌─ PARABRICKS GPU-ACCELERATED PREPROCESSING ─────────────────────────────┐
            // │ PARABRICKS GPU加速预处理                                                 │
            // │                                                                          │
            // │ Parabricks uses NVIDIA GPUs to dramatically speed up alignment and      │
            // │ processing. Instead of hours on a CPU, the same work takes minutes.     │
            // │ This is like using a supercomputer instead of a laptop.                │
            // │                                                                          │
            // │ Parabricks使用NVIDIA GPU显著加快比对和处理速度。                        │
            // │ 与CPU上的数小时相比，相同的工作需要几分钟。                             │
            // │ 这就像使用超级计算机而不是笔记本电脑。                                   │
            // └──────────────────────────────────────────────────────────────────────────┘
            // PREPROCESSING WITH PARABRICKS
            FASTQ_PREPROCESS_PARABRICKS(
                input_fastq,
                fasta,
                index_alignment,
                intervals_bed_combined,
                known_sites_indels,
                channel.value("cram"),
            )

            // Gather preprocessing output
            cram_variant_calling = channel.empty()
            cram_variant_calling = cram_variant_calling.mix(FASTQ_PREPROCESS_PARABRICKS.out.cram)

            // Gather used softwares versions
            reports = reports.mix(FASTQ_PREPROCESS_PARABRICKS.out.reports)
            versions = versions.mix(FASTQ_PREPROCESS_PARABRICKS.out.versions)
        }
        else {
            // ┌─ GATK STANDARD PREPROCESSING ─────────────────────────────────────────────┐
            // │ GATK标准预处理                                                            │
            // │                                                                            │
            // │ The gold standard approach recommended by GATK and BROAD Institute.      │
            // │ Slower than Parabricks but produces publication-quality results and     │
            // │ is widely accepted in clinical genomics.                                │
            // │                                                                            │
            // │ The key advantage of GATK BQSR: it uses reference bases known to be     │
            // │ correct to teach the model "when we see this sequencing chemistry,      │
            // │ quality scores are usually overestimated by 2 points". This training    │
            // │ improves all subsequent variant calling.                                │
            // │                                                                            │
            // │ GATK推荐的金标准方法，由BROAD Institute推荐。                             │
            // │ 比Parabricks慢，但能产生出版质量的结果，                                  │
            // │ 并在临床基因组学中广泛接受。                                             │
            // │                                                                            │
            // │ GATK BQSR的关键优势：它使用已知正确的参考碱基来教授模型                   │
            // │ "当我们看到这个测序化学时，质量分数通常被高估2个点"。                    │
            // │ 这个训练改进了所有后续的变异检测。                                       │
            // └────────────────────────────────────────────────────────────────────────────┘
            // PREPROCESSING
            FASTQ_PREPROCESS_GATK(
                input_fastq,
                input_sample,
                dict,
                fasta,
                fasta_fai,
                index_alignment,
                intervals_and_num_intervals,
                intervals_for_preprocessing,
                known_sites_indels,
                known_sites_indels_tbi,
                bbsplit_index,
            )

            // Gather preprocessing output
            cram_variant_calling = channel.empty()
            cram_variant_calling = cram_variant_calling.mix(FASTQ_PREPROCESS_GATK.out.cram_variant_calling)

            // Gather used softwares versions
            reports = reports.mix(FASTQ_PREPROCESS_GATK.out.reports)
            versions = versions.mix(FASTQ_PREPROCESS_GATK.out.versions)
        }
    }

    if (step == 'variant_calling') {

        cram_variant_calling = channel.empty().mix(input_sample)
    }

    if (step == 'annotate') {

        cram_variant_calling = channel.empty()
    }

    // ╔══════════════════════════════════════════════════════════════════════════════╗
    // ║  STEP 2B: CRAM QUALITY CONTROL (Post-Preprocessing QC)                      ║
    // ║  步骤2B：CRAM质量控制（预处理后QC）                                          ║
    // ╚══════════════════════════════════════════════════════════════════════════════╝
    //
    // After alignment, we check:
    // - Are reads properly aligned to the genome?
    // - Is there contamination from other samples? (NGSCheckmate test)
    // - What's the coverage distribution? (are some regions under-covered?)
    //
    // RUN CRAM QC on the recalibrated CRAM files or when starting from step variant calling.
    // NGSCheckmate should be run also on non-recalibrated CRAM files
    // NGSCheckmate应该也在非重新校准的CRAM文件上运行
    CRAM_SAMPLEQC(
        cram_variant_calling,
        ngscheckmate_bed,
        fasta,
        skip_tools.split(',').contains('baserecalibrator'),
        intervals_for_preprocessing,
    )

    reports = reports.mix(CRAM_SAMPLEQC.out.reports)
    versions = versions.mix(CRAM_SAMPLEQC.out.versions)

    if (tools) {

        bam_variant_calling = channel.empty()

        // ┌─ Format Conversion for Tools Requiring BAM ───────────────────────────────┐
        // │ 为需要BAM的工具进行格式转换                                                │
        // │                                                                            │
        // │ Some specialized tools (CNVkit for copy number variation, MSIsensor2 for  │
        // │ microsatellite instability, MUSE for variant calling) require uncompressed│
        // │ BAM format instead of CRAM. Since CRAM saves ~80% storage but these tools │
        // │ can't use it, we convert back to BAM for these specific analyses.        │
        // │                                                                            │
        // │ 一些专门的工具（用于拷贝数变异的CNVkit、用于微卫星不稳定性的MSIsensor2、  │
        // │ 用于变异检测的MUSE）需要未压缩的BAM格式而不是CRAM。                        │
        // │ 由于CRAM节省了约80%的存储空间，但这些工具无法使用它，                      │
        // │ 我们为这些特定分析转换回BAM。                                             │
        // └────────────────────────────────────────────────────────────────────────────┘
        //  For cnvkit, msisensor2 and muse we need to use bam input and not cram
        if (tools.split(',').contains('cnvkit') || tools.split(',').contains('msisensor2') || tools.split(',').contains('muse')) {

            // Differentiate between bam and cram files
            cram_variant_calling_status_tmp = cram_variant_calling.branch { meta, file, index ->
                bam: file.toString().endsWith('.bam')
                cram: file.toString().endsWith('.cram')
            }

            // convert cram files
            CRAM_TO_BAM(cram_variant_calling_status_tmp.cram, fasta, fasta_fai)

            // gather all bam files
            bam_variant_calling = CRAM_TO_BAM.out.bam
                .join(CRAM_TO_BAM.out.bai, by: [0])
                .mix(cram_variant_calling_status_tmp.bam)
                .map { meta, bam, bai ->
                    [meta + [data_type: 'bam'], bam, bai]
                }

            versions = versions.mix(CRAM_TO_BAM.out.versions)
        }

        // ╔══════════════════════════════════════════════════════════════════════════════╗
        // ║  STEP 3: SAMPLE CLASSIFICATION                                             ║
        // ║  步骤3：样本分类                                                             ║
        // ╚══════════════════════════════════════════════════════════════════════════════╝
        //
        // We need to handle samples differently based on their type:
        //
        // status == 0: Normal/germline samples
        //             These are from healthy individuals or represent normal tissue
        //             Used for germline variant detection (inherited mutations)
        //
        // status == 1: Tumor samples
        //             Can be analyzed alone (tumor-only) or with matched normal
        //             Tumor-only requires filtering normal variants from other people
        //             Paired analysis is best: compares tumor vs patient's own normal
        //
        // The logic here is like sorting mail by destination - each type goes to
        // different downstream analyses.
        //
        // 我们需要根据样本类型以不同的方式处理：
        //
        // status == 0: 正常/生殖系样本
        //             来自健康个体或代表正常组织
        //             用于生殖系变异检测（遗传突变）
        //
        // status == 1: 肿瘤样本
        //             可以单独分析（仅肿瘤）或与匹配的正常样本
        //             肿瘤单独需要过滤来自其他人的正常变异
        //             配对分析最好：比较肿瘤与患者自己的正常
        //
        // 这里的逻辑就像按目的地分类邮件 - 每种类型都会进行
        // 不同的下游分析。

        // Logic to separate germline samples, tumor samples with no matched normal, and combine tumor-normal pairs
        cram_variant_calling_status = cram_variant_calling.branch { meta, file, index ->
            normal: meta.status == 0
            tumor: meta.status == 1
        }

        // Follow the same logic with bam as we have with cram
        bam_variant_calling_status = bam_variant_calling.branch { meta, file, index ->
            normal: meta.status == 0
            tumor: meta.status == 1
        }

        // All Germline samples
        cram_variant_calling_normal_to_cross = cram_variant_calling_status.normal.map { meta, cram, crai -> [meta.patient, meta, cram, crai] }
        bam_variant_calling_normal_to_cross = bam_variant_calling_status.normal.map { meta, bam, bai -> [meta.patient, meta, bam, bai] }

        // All tumor samples
        cram_variant_calling_pair_to_cross = cram_variant_calling_status.tumor.map { meta, cram, crai -> [meta.patient, meta, cram, crai] }
        bam_variant_calling_pair_to_cross = bam_variant_calling_status.tumor.map { meta, bam, bai -> [meta.patient, meta, bam, bai] }

        // ┌─ Tumor-Only Sample Logic ─────────────────────────────────────────────────┐
        // │ 肿瘤单独样本逻辑                                                            │
        // │                                                                            │
        // │ Tumor-only analysis is used when:                                         │
        // │ - Matched normal tissue wasn't collected                                  │
        // │ - Patient has no germline DNA available                                   │
        // │ - Resources don't allow paired analysis                                   │
        // │                                                                            │
        // │ Challenge: How do we know if a variant is somatic (cancer-specific)       │
        // │ vs. germline (inherited by the patient since birth)?                      │
        // │                                                                            │
        // │ Solution: We use the germline_resource (gnomAD, ClinVar) to filter out   │
        // │ common germline variants. Any remaining variant is likely somatic.        │
        // │                                                                            │
        // │ 在以下情况下使用肿瘤单独分析：                                             │
        // │ - 未收集匹配的正常组织                                                    │
        // │ - 患者没有可用的生殖系DNA                                                  │
        // │ - 资源不允许进行配对分析                                                  │
        // │                                                                            │
        // │ 挑战：我们如何知道一个变异是体细胞（癌症特异性）                           │
        // │ 还是生殖系（患者从出生起就继承的）？                                      │
        // │                                                                            │
        // │ 解决方案：我们使用germline_resource（gnomAD, ClinVar）过滤掉                │
        // │ 常见的生殖系变异。任何剩余的变异都可能是体细胞的。                         │
        // └────────────────────────────────────────────────────────────────────────────┘

        // Tumor only samples
        // 1. Group together all tumor samples by patient ID [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ] ]

        // Downside: this only works by waiting for all tumor samples to finish preprocessing, since no group size is provided
        cram_variant_calling_tumor_grouped = cram_variant_calling_pair_to_cross.groupTuple()
        bam_variant_calling_tumor_grouped = bam_variant_calling_pair_to_cross.groupTuple()

        // 2. Join with normal samples, in each channel there is one key per patient now. Patients without matched normal end up with: [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ], null ]
        cram_variant_calling_tumor_joined = cram_variant_calling_tumor_grouped.join(cram_variant_calling_normal_to_cross, failOnDuplicate: true, remainder: true)
        bam_variant_calling_tumor_joined = bam_variant_calling_tumor_grouped.join(bam_variant_calling_normal_to_cross, failOnDuplicate: true, remainder: true)

        // 3. Filter out entries with last entry null
        // This keeps tumor-only samples (without matched normal)
        cram_variant_calling_tumor_filtered = cram_variant_calling_tumor_joined.filter { it -> !(it.last()) }
        bam_variant_calling_tumor_filtered = bam_variant_calling_tumor_joined.filter { it -> !(it.last()) }

        // 4. Transpose [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ] ] back to [ patient1, meta1, [ cram1, crai1 ], null ] [ patient1, meta2, [ cram2, crai2 ], null ]
        // and remove patient ID field & null value for further processing [ meta1, [ cram1, crai1 ] ] [ meta2, [ cram2, crai2 ] ]
        cram_variant_calling_tumor_only = cram_variant_calling_tumor_filtered.transpose().map { it -> [it[1], it[2], it[3]] }
        bam_variant_calling_tumor_only = bam_variant_calling_tumor_filtered.transpose().map { it -> [it[1], it[2], it[3]] }

        if (params.only_paired_variant_calling) {
            // ┌─ Paired-Only Mode: Filter Normal Samples ──────────────────────────────┐
            // │ 仅配对模式：过滤正常样本                                                 │
            // │                                                                        │
            // │ When only_paired_variant_calling is enabled, we only analyze samples  │
            // │ that have matched pairs (normal + tumor). Normal samples without a    │
            // │ corresponding tumor are discarded.                                    │
            // │                                                                        │
            // │ 启用only_paired_variant_calling时，我们仅分析具有匹配对                │
            // │ （正常+肿瘤）的样本。没有对应肿瘤的正常样本被丢弃。                    │
            // └────────────────────────────────────────────────────────────────────────┘
            // Normal only samples

            // 1. Join with tumor samples, in each channel there is one key per patient now. Patients without matched tumor end up with: [ patient1, [ meta1 ], [ cram1, crai1 ], null ] as there is only one matched normal possible
            cram_variant_calling_normal_joined = cram_variant_calling_normal_to_cross.join(cram_variant_calling_tumor_grouped, failOnDuplicate: true, remainder: true)
            bam_variant_calling_normal_joined = bam_variant_calling_normal_to_cross.join(bam_variant_calling_tumor_grouped, failOnDuplicate: true, remainder: true)

            // 2. Filter out entries with last entry null
            cram_variant_calling_normal_filtered = cram_variant_calling_normal_joined.filter { it -> !(it.last()) }
            bam_variant_calling_normal_filtered = bam_variant_calling_normal_joined.filter { it -> !(it.last()) }

            // 3. Remove patient ID field & null value for further processing [ meta1, [ cram1, crai1 ] ] [ meta2, [ cram2, crai2 ] ] (no transposing needed since only one normal per patient ID)
            cram_variant_calling_status_normal = cram_variant_calling_normal_filtered.map { it -> [it[1], it[2], it[3]] }
            bam_variant_calling_status_normal = bam_variant_calling_normal_filtered.map { it -> [it[1], it[2], it[3]] }
        }
        else {
            cram_variant_calling_status_normal = cram_variant_calling_status.normal
            bam_variant_calling_status_normal = bam_variant_calling_status.normal
        }

        // ┌─ Tumor-Normal Pairing ────────────────────────────────────────────────────┐
        // │ 肿瘤-正常配对                                                              │
        // │                                                                            │
        // │ The CROSS operation combines each tumor with each normal from the same   │
        // │ patient. This is the BEST scenario for somatic variant calling because:  │
        // │                                                                            │
        // │ - Any variant in the tumor that's NOT in normal = likely somatic         │
        // │ - We can measure contamination (normal DNA in tumor sample)              │
        // │ - We can detect subclonal populations (some cells have variant, others   │
        // │   don't) - these are much harder to call from tumor-only data           │
        // │                                                                            │
        // │ Special case: Multi-tumor per patient (recurrent tumors, metastases)     │
        // │ - One normal can be paired with multiple tumor samples                   │
        // │ - This is why we use CROSS instead of JOIN                              │
        // │                                                                            │
        // │ CROSS操作将每个肿瘤与同一患者的每个正常组织结合。                          │
        // │ 这是体细胞变异检测的最佳情景，因为：                                      │
        // │                                                                            │
        // │ - 肿瘤中任何不在正常的变异=可能是体细胞                                   │
        // │ - 我们可以测量污染（肿瘤样本中的正常DNA）                                │
        // │ - 我们可以检测亚克隆群体（某些细胞有变异，其他则没有）                    │
        // │   - 这些从仅肿瘤数据中更难调用                                            │
        // │                                                                            │
        // │ 特殊情况：每个患者的多个肿瘤（复发性肿瘤、转移）                          │
        // │ - 一个正常可以与多个肿瘤样本配对                                          │
        // │ - 这就是为什么我们使用CROSS而不是JOIN                                    │
        // └────────────────────────────────────────────────────────────────────────────┘

        // Tumor - normal pairs
        // Use cross to combine normal with all tumor samples, i.e. multi tumor samples from recurrences
        cram_variant_calling_pair = cram_variant_calling_normal_to_cross
            .cross(cram_variant_calling_pair_to_cross)
            .map { normal, tumor ->
                def meta = [:]

                meta.id = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id = normal[1].sample
                meta.patient = normal[0]
                meta.sex = normal[1].sex
                meta.tumor_id = tumor[1].sample
                meta.contamination = tumor[1].contamination

                [meta, normal[2], normal[3], tumor[2], tumor[3]]
            }
        bam_variant_calling_pair = bam_variant_calling_normal_to_cross
            .cross(bam_variant_calling_pair_to_cross)
            .map { normal, tumor ->
                def meta = [:]

                meta.id = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id = normal[1].sample
                meta.patient = normal[0]
                meta.sex = normal[1].sex
                meta.tumor_id = tumor[1].sample

                [meta, normal[2], normal[3], tumor[2], tumor[3]]
            }

        // ╔══════════════════════════════════════════════════════════════════════════════╗
        // ║  STEP 4: VARIANT CALLING (THREE PATHWAYS)                                   ║
        // ║  步骤4：变异检测（三条路径）                                                 ║
        // ╚══════════════════════════════════════════════════════════════════════════════╝
        //
        // Based on sample type, we use different variant calling strategies:
        // 根据样本类型，我们使用不同的变异检测策略：

        // ┌─ PATHWAY 1: GERMLINE VARIANT CALLING ──────────────────────────────────┐
        // │ 路径1：生殖系变异检测                                                    │
        // │                                                                        │
        // │ For normal/healthy samples, we detect inherited genetic variants.    │
        // │ These tools find SNVs (single base changes) and indels.              │
        // │                                                                        │
        // │ Tools used:                                                           │
        // │ - HaplotypeCaller: GATK's standard germline variant caller           │
        // │ - DeepVariant: AI-based variant caller (excellent accuracy)          │
        // │ - FreeBayes: Conservative caller, good for low-coverage regions      │
        // │ - Strelka2: Fast caller, originally for tumor but good for germline  │
        // │ - TIDDIT: Detects structural variants (large deletions/duplications) │
        // │                                                                        │
        // │ 对于正常/健康的样本，我们检测遗传的遗传变异。                         │
        // │ 这些工具找到SNVs（单碱基变化）和插入缺失。                             │
        // │                                                                        │
        // │ 使用的工具：                                                          │
        // │ - HaplotypeCaller: GATK的标准生殖系变异调用器                        │
        // │ - DeepVariant: 基于AI的变异检测器（准确性优异）                       │
        // │ - FreeBayes: 保守的调用器，对低覆盖范围区域很好                      │
        // │ - Strelka2: 快速调用器，最初用于肿瘤，但对生殖系很好                 │
        // │ - TIDDIT: 检测结构变异（大删除/重复）                                 │
        // └────────────────────────────────────────────────────────────────────────┘
        // GERMLINE VARIANT CALLING
        //   No bwa index for TIDDIT
        //   intervals handling
        //     intervals_bed_combined: [] if no_intervals, else interval_bed_combined
        //     intervals_bed_gz_tbi_combined, [] if no_intervals, else interval_bed_combined_gz_tbi
        //     intervals_bed_combined_for_variant_calling, no_intervals.bed if no intervals, else interval_bed_combined.bed
        //   skip_tools.split(',').contains('haplotypecaller_filter') is true if filtering should be skipped
        BAM_VARIANT_CALLING_GERMLINE_ALL(
            tools,
            skip_tools,
            bam_variant_calling_status_normal,
            cram_variant_calling_status_normal,
            [[id: 'bwa'], []],
            cnvkit_reference,
            dbsnp,
            dbsnp_tbi,
            dbsnp_vqsr,
            dict,
            fasta,
            fasta_fai,
            intervals_and_num_intervals,
            intervals_bed_combined,
            intervals_bed_gz_tbi_combined,
            intervals_bed_combined_for_variant_calling,
            intervals_bed_gz_tbi_and_num_intervals,
            known_indels_vqsr,
            known_sites_indels,
            known_sites_indels_tbi,
            known_sites_snps,
            known_sites_snps_tbi,
            known_snps_vqsr,
            params.joint_germline,
            skip_tools.split(',').contains('haplotypecaller_filter'),
            params.sentieon_haplotyper_emit_mode,
            params.sentieon_dnascope_emit_mode,
            params.sentieon_dnascope_pcr_indel_model,
            sentieon_dnascope_model,
        )

        // ┌─ PATHWAY 2: TUMOR-ONLY VARIANT CALLING ───────────────────────────────┐
        // │ 路径2：肿瘤单独变异检测                                                │
        // │                                                                        │
        // │ When only tumor tissue is available, we must be very careful about   │
        // │ what we call "cancer mutations". Two main approaches:                │
        // │                                                                        │
        // │ 1. SNV/INDEL callers (Mutect2, Strelka2):                           │
        // │    - Use Panel of Normals (PoN) - variants from healthy people      │
        // │    - Any variant in PoN is likely sequencing artifact, not mutation │
        // │    - Uses germline_resource to filter common inherited variants     │
        // │                                                                        │
        // │ 2. Structural variant callers (Manta):                              │
        // │    - Detect large deletions, duplications, inversions               │
        // │    - Common in cancer but rare in normal genomes                    │
        // │                                                                        │
        // │ 3. Copy number analysis (ControlFREEC, CNVkit):                      │
        // │    - Detect regions where tumor has extra (gain) or missing (loss) │
        // │      copies of DNA - very common in cancer                         │
        // │                                                                        │
        // │ 当仅有肿瘤组织可用时，我们必须非常小心地定义"癌症突变"。               │
        // │ 两个主要方法：                                                        │
        // │                                                                        │
        // │ 1. SNV/INDEL调用器（Mutect2, Strelka2）：                           │
        // │    - 使用法线板（PoN）- 来自健康人的变异                             │
        // │    - PoN中的任何变异都可能是测序伪影，而不是突变                     │
        // │    - 使用germline_resource过滤常见的遗传变异                         │
        // │                                                                        │
        // │ 2. 结构变异调用器（Manta）：                                        │
        // │    - 检测大删除、重复、反向                                          │
        // │    - 在癌症中常见，但在正常基因组中罕见                              │
        // │                                                                        │
        // │ 3. 拷贝数分析（ControlFREEC, CNVkit）：                             │
        // │    - 检测肿瘤具有额外（增益）或缺少（丢失）DNA拷贝的区域              │
        // │      - 在癌症中非常常见                                              │
        // └────────────────────────────────────────────────────────────────────────┘
        // TUMOR ONLY VARIANT CALLING
        //   No bwa index for TIDDIT
        //   intervals handling
        //     intervals_bed_combined: [] if no_intervals, else interval_bed_combined
        //     intervals_bed_gz_tbi_combined, [] if no_intervals, else interval_bed_combined_gz_tbi
        BAM_VARIANT_CALLING_TUMOR_ONLY_ALL(
            tools,
            bam_variant_calling_tumor_only,
            cram_variant_calling_tumor_only,
            [[id: 'bwa'], []],
            cf_chrom_len,
            chr_files,
            cnvkit_reference,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            germline_resource,
            germline_resource_tbi,
            intervals_and_num_intervals,
            intervals_bed_gz_tbi_and_num_intervals,
            intervals_bed_combined,
            intervals_bed_gz_tbi_combined,
            mappability,
            msisensor2_models,
            pon,
            pon_tbi,
            params.joint_mutect2,
            params.wes,
        )

        // ┌─ PATHWAY 3: SOMATIC PAIRED VARIANT CALLING ────────────────────────────┐
        // │ 路径3：体细胞配对变异检测                                               │
        // │                                                                         │
        // │ BEST CASE: Both tumor and normal tissue available from same patient    │
        // │                                                                         │
        // │ Advantages:                                                            │
        // │ - Directly compare tumor vs this patient's own normal tissue          │
        // │ - Any variant in tumor but not in normal = definitely somatic         │
        // │ - Can measure tumor purity (% of tumor cells vs infiltrating normal)  │
        // │ - Can detect LOH (Loss of Heterozygosity) - cancer hiding one copy    │
        // │ - Can use ASCAT for allele-specific copy number (which copy has loss?)│
        // │                                                                         │
        // │ Tools:                                                                 │
        // │ - Mutect2: Specialized for detecting somatic SNVs/indels             │
        // │ - Strelka2: Fast, accurate for paired samples                        │
        // │ - Manta: Structural variants in tumor                                │
        // │ - ASCAT: Allele-specific copy number analysis                        │
        // │ - ControlFREEC: Another CNV caller, doesn't need matched normal      │
        // │ - MSIsensor/Pro: Microsatellite instability (sign of mismatch repair │
        // │                  deficiency - common in cancer)                       │
        // │                                                                         │
        // │ 最佳情况：同一患者的肿瘤和正常组织都可用                               │
        // │                                                                         │
        // │ 优点：                                                                  │
        // │ - 直接比较肿瘤与这个患者自己的正常组织                                 │
        // │ - 肿瘤中的任何变异但不在正常的 = 绝对是体细胞                          │
        // │ - 可以测量肿瘤纯度（% 肿瘤细胞与浸润的正常）                            │
        // │ - 可以检测LOH（杂合子丧失）- 癌症隐藏一个拷贝                          │
        // │ - 可以使用ASCAT进行等位基因特异性拷贝数（哪个拷贝有丢失？）            │
        // │                                                                         │
        // │ 工具：                                                                  │
        // │ - Mutect2: 专门用于检测体细胞SNVs/插入缺失                            │
        // │ - Strelka2: 快速、准确的配对样本                                      │
        // │ - Manta: 肿瘤中的结构变异                                             │
        // │ - ASCAT: 等位基因特异性拷贝数分析                                      │
        // │ - ControlFREEC: 另一个CNV调用器，不需要匹配的正常                    │
        // │ - MSIsensor/Pro: 微卫星不稳定性（错配修复缺陷的标志 - 癌症中常见）  │
        // └─────────────────────────────────────────────────────────────────────────┘
        // PAIR VARIANT CALLING
        //   No bwa index for TIDDIT
        //   intervals handling
        //     intervals_bed_combined: [] if no_intervals, else interval_bed_combined
        //     intervals_bed_gz_tbi_combined, [] if no_intervals, else interval_bed_combined_gz_tbi
        BAM_VARIANT_CALLING_SOMATIC_ALL(
            tools,
            bam_variant_calling_pair,
            cram_variant_calling_pair,
            [[id: 'bwa'], []],
            cf_chrom_len,
            chr_files,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            germline_resource,
            germline_resource_tbi,
            intervals_and_num_intervals,
            intervals_bed_gz_tbi_and_num_intervals,
            intervals_bed_combined,
            intervals_bed_gz_tbi_combined,
            mappability,
            msisensorpro_scan,
            pon,
            pon_tbi,
            ascat_alleles,
            ascat_loci,
            ascat_loci_gc,
            ascat_loci_rt,
            params.joint_mutect2,
            params.wes,
        )

        // ╔══════════════════════════════════════════════════════════════════════════════╗
        // ║  STEP 5: VARIANT CALLING QUALITY CONTROL                                    ║
        // ║  步骤5：变异检测质量控制                                                      ║
        // ╚══════════════════════════════════════════════════════════════════════════════╝
        //
        // After finding variants, we assess their quality:
        // - Transition/Transversion ratio (Ts/Tv): Natural bias in genetics
        //   Real mutations are ~2:1 Ts/Tv. Low ratio suggests sequencing errors
        // - Filter summary: How many variants passed quality thresholds?
        // - Coverage statistics: Are variants in well-covered regions?
        //
        // 找到变异后，我们评估其质量：
        // - 转变/颠换比（Ts/Tv）：遗传学中的自然偏差
        //   真正的突变约为2:1 Ts/Tv。低比率表明测序错误
        // - 过滤摘要：有多少变异通过了质量阈值？
        // - 覆盖统计：变异是否在覆盖良好的区域中？

        // QC on raw variant calls
        VCF_QC_BCFTOOLS_VCFTOOLS(
            BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_all.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_ALL.out.vcf_all).mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.vcf_all),
            intervals_bed_combined,
        )

        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.bcftools_stats.collect { _meta, stats -> [stats] })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_counts.collect { _meta, counts -> [counts] })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_qual.collect { _meta, qual -> [qual] })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_filter_summary.collect { _meta, summary -> [summary] })
        reports = reports.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.out_indexcov.collect { _meta, indexcov -> indexcov.flatten() })
        reports = reports.mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.out_indexcov.collect { _meta, indexcov -> indexcov.flatten() })

        // ╔══════════════════════════════════════════════════════════════════════════════╗
        // ║  STEP 6: POST-VARIANT CALLING PROCESSING                                    ║
        // ║  步骤6：变异检测后处理                                                        ║
        // ╚══════════════════════════════════════════════════════════════════════════════╝
        //
        // VCF files from different callers may:
        // - Have different formats (need normalization)
        // - Include low-confidence variants (need filtering)
        // - Need merging for joint analysis
        // - Need left-alignment (move indels to canonical positions)
        //
        // This step also supports consensus calling - if multiple callers agree on
        // a variant, confidence goes up dramatically.
        //
        // VCF文件来自不同的调用器可能：
        // - 具有不同的格式（需要标准化）
        // - 包括低置信度变异（需要过滤）
        // - 需要合并以进行联合分析
        // - 需要左对齐（将插入缺失移动到规范位置）
        //
        // 该步骤还支持共识调用 - 如果多个调用器同意一个变异，
        // 置信度会急剧上升。

        // POST VARIANTCALLING
        POST_VARIANTCALLING(
            tools,
            cram_variant_calling_status_normal,
            BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_all,
            BAM_VARIANT_CALLING_GERMLINE_ALL.out.tbi_all,
            cram_variant_calling_tumor_only,
            BAM_VARIANT_CALLING_TUMOR_ONLY_ALL.out.vcf_all,
            BAM_VARIANT_CALLING_TUMOR_ONLY_ALL.out.tbi_all,
            cram_variant_calling_pair,
            BAM_VARIANT_CALLING_SOMATIC_ALL.out.vcf_all,
            BAM_VARIANT_CALLING_SOMATIC_ALL.out.tbi_all,
            fasta,
            fasta_fai,
            params.concatenate_vcfs,
            params.filter_vcfs,
            params.snv_consensus_calling,
            params.normalize_vcfs,
            params.varlociraptor_chunk_size,
            varlociraptor_scenario_germline,
            varlociraptor_scenario_somatic,
            varlociraptor_scenario_tumor_only,
        )

        // Gather vcf files for annotation and QC
        // POST_VARIANTCALLING always outputs VCFs - either processed or pass-through originals
        vcf_to_annotate = POST_VARIANTCALLING.out.vcfs

        CHANNEL_VARIANT_CALLING_CREATE_CSV(vcf_to_annotate, params.outdir)

        // Gather used variant calling softwares versions
        versions = versions.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.versions)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.versions)
        versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_ALL.out.versions)
        versions = versions.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.versions)
        versions = versions.mix(POST_VARIANTCALLING.out.versions)

        // ╔══════════════════════════════════════════════════════════════════════════════╗
        // ║  STEP 7: VARIANT ANNOTATION                                                ║
        // ║  步骤7：变异注释                                                              ║
        // ╚══════════════════════════════════════════════════════════════════════════════╝
        //
        // At this point, we have lists of genomic coordinates with variants, but
        // we don't know WHAT THOSE VARIANTS MEAN biologically.
        //
        // This step adds critical information:
        //
        // 1. GENE ASSIGNMENT: Which gene is affected?
        //    - UCSC: maps coordinates to genes using transcripts
        //
        // 2. CONSEQUENCE PREDICTION: What does this mutation do?
        //    - Synonymous (silent): doesn't change amino acid - usually harmless
        //    - Missense: changes amino acid - may affect protein function
        //    - Frameshift: changes reading frame - usually very damaging
        //    - Splice site: affects RNA splicing - often harmful
        //    - Stop codon: truncates protein - usually very harmful
        //
        // 3. PATHOGENICITY PREDICTION:
        //    - SIFT: "Is this amino acid change tolerated?"
        //    - PolyPhen: "Is this probably damaging?"
        //    - These use machine learning on protein structure
        //
        // 4. POPULATION FREQUENCY:
        //    - gnomAD: Variant frequency in 140,000+ people
        //    - ClinVar: Clinical significance from NCBI
        //    - dbSNP: Known variants (most common variants known)
        //
        // 5. DISEASE ASSOCIATIONS:
        //    - COSMIC: Variants found in cancer samples
        //    - ClinVar: Pathogenic variants associated with diseases
        //
        // Tools used:
        // - VEP (Variant Effect Predictor): Comprehensive annotation
        // - SnpEff: Fast, reliable annotations
        // - BCFtools: Adds custom annotations from other databases
        // - SnpSift: Adds clinical significance annotations
        //
        // 此时，我们有带有变异的基因组坐标列表，但
        // 我们不知道这些变异在生物学上意味着什么。
        //
        // 此步骤添加关键信息：
        //
        // 1. 基因分配：哪个基因受到影响？
        //    - UCSC：使用转录本将坐标映射到基因
        //
        // 2. 后果预测：这个突变做什么？
        //    - 同义（沉默）：不改变氨基酸 - 通常无害
        //    - 错义：改变氨基酸 - 可能影响蛋白质功能
        //    - 移码：改变阅读框 - 通常非常有害
        //    - 剪接位点：影响RNA剪接 - 通常有害
        //    - 停止密码子：截断蛋白质 - 通常非常有害
        //
        // 3. 致病性预测：
        //    - SIFT："这个氨基酸改变被容许吗？"
        //    - PolyPhen："这可能有害吗？"
        //    - 这些使用蛋白质结构上的机器学习
        //
        // 4. 群体频率：
        //    - gnomAD：140,000+人群中的变异频率
        //    - ClinVar：来自NCBI的临床意义
        //    - dbSNP：已知变异（最常见的已知变异）
        //
        // 5. 疾病关联：
        //    - COSMIC：在癌症样本中发现的变异
        //    - ClinVar：与疾病相关的致病性变异
        //
        // 使用的工具：
        // - VEP（变异效果预测器）：全面的注释
        // - SnpEff：快速、可靠的注释
        // - BCFtools：从其他数据库添加自定义注释
        // - SnpSift：添加临床意义注释

        // ANNOTATE
        if (step == 'annotate') {
            vcf_to_annotate = input_sample
        }

        if (tools.split(',').contains('merge') || tools.split(',').contains('snpeff') || tools.split(',').contains('vep') || tools.split(',').contains('bcfann') || tools.split(',').contains('snpsift')) {

            vep_fasta = params.vep_include_fasta ? fasta : [[id: 'null'], []]

            VCF_ANNOTATE_ALL(
                vcf_to_annotate.map { meta, vcf -> [meta + [file_name: vcf.baseName], vcf] },
                vep_fasta,
                tools,
                snpeff_db,
                snpeff_cache,
                vep_genome,
                vep_species,
                vep_cache_version,
                vep_cache,
                vep_extra_files,
                bcftools_annotations,
                bcftools_annotations_tbi,
                bcftools_columns,
                bcftools_header_lines,
                snpsift_db,
            )

            // Gather used softwares versions
            versions = versions.mix(VCF_ANNOTATE_ALL.out.versions)
            reports = reports.mix(VCF_ANNOTATE_ALL.out.reports)
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 8: SOFTWARE VERSION TRACKING & MULTIQC REPORT
        步骤8：软件版本跟踪和MultiQC报告
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //
    // Collate and save software versions
    // 整理和保存软件版本
    //
    // For reproducibility, document exactly which version of every tool was used
    // This is critical for:
    // - Reproducing results in another lab
    // - Reporting in publications
    // - Regulatory compliance (FDA, EMA for clinical use)
    //
    // 为了可重现性，记录使用的每个工具的确切版本
    // 这对以下方面至关重要：
    // - 在另一个实验室中重现结果
    // - 在出版物中报告
    // - 法规遵从性（FDA、EMA用于临床使用）
    //
    def version_yaml = channel.empty()
    if (!(skip_tools.split(',').contains('versions'))) {
        version_yaml = softwareVersionsToYAML(
            softwareVersions: versions.mix(channel.topic("versions")),
            nextflowVersion: workflow.nextflow.version,
        ).collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'sarek_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
    }

    // ╔══════════════════════════════════════════════════════════════════════════════╗
    // ║  STEP 9: MULTIQC - COMPREHENSIVE QUALITY REPORT                            ║
    // ║  步骤9：MultiQC - 综合质量报告                                               ║
    // ╚══════════════════════════════════════════════════════════════════════════════╝
    //
    // MultiQC aggregates all QC reports into ONE beautiful interactive HTML report.
    // Instead of looking at 20 different output files, you see everything together:
    //
    // - FASTQC: Input read quality
    // - Samtools: Alignment statistics
    // - BCFtools: Variant calling statistics
    // - Picard: Duplicate rates, coverage
    // - And many more...
    //
    // Think of it like a dashboard - one screen showing all important metrics.
    // You can quickly see if something went wrong with preprocessing, alignment,
    // or variant calling.
    //
    // MultiQC将所有QC报告聚合为一个漂亮的交互式HTML报告。
    // 而不是查看20个不同的输出文件，您看到所有内容一起：
    //
    // - FASTQC：输入读取质量
    // - Samtools：比对统计
    // - BCFtools：变异调用统计
    // - Picard：重复率、覆盖
    // - 以及许多更多...
    //
    // 可以将其视为仪表板 - 一个屏幕显示所有重要指标。
    // 您可以快速看到预处理、比对或变异检测是否出现问题。
    //
    if (!(skip_tools.split(',').contains('multiqc'))) {
        ch_multiqc_config = channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
        ch_multiqc_logo = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()

        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(version_yaml)
        ch_multiqc_files = ch_multiqc_files.mix(reports)
        ch_multiqc_files = ch_multiqc_files.mix(channel.topic("multiqc_files").map { _meta, _process, _tool, report -> report })

        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

        MULTIQC(
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            [],
        )
        multiqc_publish = MULTIQC.out.data.mix(MULTIQC.out.plots, MULTIQC.out.report)
        multiqc_report = MULTIQC.out.report.toList()
    }

    emit:
    multiqc_report // channel: /path/to/multiqc_report.html
    multiqc_publish
    versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  HELPER FUNCTION: Add Read Group Information                               ║
// ║  辅助函数：添加读取组信息                                                     ║
// ╚══════════════════════════════════════════════════════════════════════════════╝
//
// Read group information is metadata that travels with each DNA read through
// the pipeline. It tells downstream tools:
// - Which sample the read came from
// - Which sequencing run/flowcell
// - Which sequencing lane
// - Sequencing center/platform
//
// Why this matters:
// - Allows filtering of low-quality lanes
// - Tracks contamination sources
// - Required for GATK best practices
// - Essential for matching reads back to lanes for troubleshooting
//
// Format: @RG\tID:...\tSM:...\tLB:...\tPL:...
// Each \t is a tab-separated field
//
// 读取组信息是随着每个DNA读取通过流程传输的元数据。
// 它告诉下游工具：
// - 读取来自哪个样本
// - 哪个测序运行/流动池
// - 哪个测序车道
// - 测序中心/平台
//
// 为什么这很重要：
// - 允许过滤低质量车道
// - 跟踪污染源
// - GATK最佳实践需要
// - 对于匹配读取回车道以进行故障排除至关重要
//
// 格式：@RG\tID:...\tSM:...\tLB:...\tPL:...
// 每个\t都是制表符分隔的字段
//
// Add readgroup to meta and remove lane
def addReadgroupToMeta(meta, files) {
    def CN = params.seq_center ? "CN:${params.seq_center}\\t" : ''
    def flowcell = flowcellLaneFromFastq(files[0])

    // Check if flowcell ID matches
    // If R1 and R2 come from different flowcells, something is wrong with the input
    if (flowcell && flowcell != flowcellLaneFromFastq(files[1])) {
        error("Flowcell ID does not match for paired reads of sample ${meta.id} - ${files}")
    }

    // If we cannot read the flowcell ID from the fastq file, then we don't use it
    def sample_lane_id = flowcell ? "${flowcell}.${meta.sample}.${meta.lane}" : "${meta.sample}.${meta.lane}"

    // Don't use a random element for ID, it breaks resuming
    // UMI (Unique Molecular Identifiers) = short barcodes to track duplicates better
    def read_group = params.umi_read_structure
        ? "\"@RG\\tID:${meta.sample}\\t${CN}PU:consensus\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
        : "\"@RG\\tID:${sample_lane_id}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

    meta = meta - meta.subMap('lane') + [read_group: read_group.toString(), sample_lane_id: sample_lane_id.toString()]
    return [meta, files]
}

// ╔══════════════════════════════════════════════════════════════════════════════╗
// ║  HELPER FUNCTION: Extract Flowcell & Lane from FASTQ Header               ║
// ║  辅助函数：从FASTQ标头中提取流动池和车道                                      ║
// ╚══════════════════════════════════════════════════════════════════════════════╝
//
// Every FASTQ file has a header line that looks like:
// @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x-pos>:<y-pos> <description>
//
// Example (Illumina NextSeq):
// @NS500564:123:HCNCTALXX:1:11001:1000:1000 1:N:0:
//
// We parse this to extract:
// - Flowcell ID: identifies which physical flowcell (chip) the reads came from
// - Lane number: which lane on the flowcell
//
// Why extract this?
// - Quality can differ between lanes (one lane might be dusty)
// - Can identify which specific chip had a problem
// - Important for multiplexed samples (multiple samples on one flowcell)
//
// 每个FASTQ文件都有一个看起来像的标头行：
// @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x-pos>:<y-pos> <description>
//
// 例子（Illumina NextSeq）：
// @NS500564:123:HCNCTALXX:1:11001:1000:1000 1:N:0:
//
// 我们解析这个以提取：
// - 流动池ID：标识读取来自哪个物理流动池（芯片）
// - 车道号：流动池上的哪个车道
//
// 为什么提取这个？
// - 质量在车道之间可能不同（一个车道可能是多尘的）
// - 可以识别哪个特定芯片有问题
// - 对于多路复用样本（一个流动池上的多个样本）很重要
//
// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // First line of FASTQ file contains sequence identifier plus optional description
    def firstLine = readFirstLineOfFastq(path)
    def flowcell_id = null

    // Expected format from ILLUMINA
    // cf https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers
    // Five fields:
    // @<instrument>:<lane>:<tile>:<x-pos>:<y-pos>...
    // Seven fields or more (from CASAVA 1.8+):
    // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>..."

    def fields = firstLine ? firstLine.split(':') : []
    if (fields.size() == 5) {
        // Get the instrument name as flowcell ID
        flowcell_id = fields[0].substring(1)
    }
    else if (fields.size() >= 7) {
        // Get the actual flowcell ID
        flowcell_id = fields[2]
    }
    else if (fields.size() != 0) {
        log.warn("FASTQ file(${path}): Cannot extract flowcell ID from ${firstLine}")
    }
    return flowcell_id
}

// Get first line of a FASTQ file
// 获取FASTQ文件的第一行
//
// We need to read gzip-compressed FASTQ files line by line.
// This function opens the file, decompresses it, and reads just the first line.
//
// 我们需要逐行读取gzip压缩的FASTQ文件。
// 此函数打开文件、解压缩并仅读取第一行。
//
def readFirstLineOfFastq(path) {
    def line = null
    try {
        path.withInputStream {
            def InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
            def Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
            def BufferedReader buffered = new BufferedReader(decoder)
            line = buffered.readLine()
            assert line.startsWith('@')
        }
    }
    catch (Exception e) {
        log.warn("FASTQ file(${path}): Error streaming")
        log.warn("${e.message}")
    }
    return line
}
