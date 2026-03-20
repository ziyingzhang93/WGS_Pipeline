// ╔══════════════════════════════════════════════════════════════════╗
// ║  GATK4 HaplotypeCaller Core Module                              ║
// ║  GATK4 HaplotypeCaller 核心模块                                   ║
// ╚══════════════════════════════════════════════════════════════════╝
// This is the actual GATK HaplotypeCaller execution module
// 这是执行 GATK HaplotypeCaller 变异检测的核心模块
// Performs: Local read reassembly + Bayesian inference for variant genotyping
// 执行: 局部 reads 重组装 + 贝叶斯推断进行基因型分型
// Output: gVCF (genomic VCF with all sites) and optionally realigned BAM

// ╔══════════════════════════════════════════════════════════════════╗
// ║  Process Definition: GATK4_HAPLOTYPECALLER                      ║
// ║  进程定义: GATK4_HAPLOTYPECALLER                                  ║
// ╚══════════════════════════════════════════════════════════════════╝

process GATK4_HAPLOTYPECALLER {
    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Process Metadata and Resource Configuration                 ║
    // ║  进程元数据和资源配置                                         ║
    // ╚══════════════════════════════════════════════════════════════╝
    tag "$meta.id"
    // tag: Label for this process execution in the workflow logs
    // 为流程日志中的这个进程执行添加标签，通常是样本 ID
    // Makes it easier to track which sample is being processed
    // 便于追踪正在处理哪个样本

    label 'process_low'
    // label: Specifies compute resource requirements
    // 指定计算资源需求 (在 nextflow.config 中定义)
    // 'process_low' = low CPU, low memory allocation
    // 'process_low' = 低 CPU，低内存分配 (通常用于小任务)

    conda "${moduleDir}/environment.yml"
    // Specify conda environment with required dependencies
    // 指定包含所需依赖的 conda 环境
    // HaplotypeCaller requires: gatk4, samtools, etc.
    // HaplotypeCaller 需要: gatk4, samtools 等

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data':
        'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd' }"
    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Container Image Configuration (容器镜像配置)                ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ If using Singularity: Use local container image URL          ║
    // ║ 如果使用 Singularity: 使用本地容器镜像 URL                    ║
    // ║                                                               ║
    // ║ If using Docker (default): Use Wave/Seqera cloud image       ║
    // ║ 如果使用 Docker (默认): 使用 Wave/Seqera 云镜像              ║
    // ║                                                               ║
    // ║ These are prebuilt container images with GATK4, samtools,   ║
    // ║ and other bioinformatics tools installed                    ║
    // ║ 这些是预构建的容器镜像，已安装 GATK4、samtools 等工具       ║
    // ╚══════════════════════════════════════════════════════════════╝

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Input: Data Channels and Files (输入: 数据流和文件)         ║
    // ╚══════════════════════════════════════════════════════════════╝
    input:
    tuple val(meta),  path(input), path(input_index), path(intervals), path(dragstr_model)
    // meta: Sample metadata map (样本元数据)
    //   - meta.id: Sample identifier (样本 ID)
    //   - meta.interval_name: Name of this genomic interval (基因组区间名)
    //   - meta.num_intervals: Total number of intervals for this sample
    // input: BAM/CRAM file (比对文件)
    // input_index: BAI/CRAI index file (比对索引)
    // intervals: BED file specifying genomic region (基因组区间 BED 文件)
    // dragstr_model: Optional DragStr model file for HMM parameterization

    tuple val(meta2), path(fasta)
    // meta2: Metadata for reference genome (参考基因组元数据)
    // fasta: Reference genome sequence file (参考基因组 FASTA)

    tuple val(meta3), path(fai)
    // fai: FASTA index file (FASTA 索引文件)
    // Enables fast random access to genomic sequences
    // 加快对基因组序列的随机访问

    tuple val(meta4), path(dict)
    // dict: SAM sequence dictionary file (SAM 序列字典)
    // Contains chromosome names, lengths, and MD5 checksums
    // 包含染色体名称、长度和 MD5 校验和
    // Required by GATK for coordinate validation
    // GATK 需要这个文件进行坐标验证

    tuple val(meta5), path(dbsnp)
    // dbsnp: Known variants from dbSNP database (dbSNP 已知变异数据库)
    // Used to annotate known variants and for recalibration
    // 用于标注已知变异和重新校正
    // Optional but recommended
    // 可选但推荐使用

    tuple val(meta6), path(dbsnp_tbi)
    // dbsnp_tbi: Tabix index for dbSNP VCF file (dbSNP VCF 索引)
    // Enables fast lookup of known variants
    // 加速已知变异的查询

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Output: Generated Files (输出: 生成的文件)                  ║
    // ╚══════════════════════════════════════════════════════════════╝
    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf
    // vcf: gVCF file (genomic VCF)
    // Naming: Sample_interval_name.vcf.gz
    // Contains all sites (variant and non-variant)
    // 包含所有位点 (变异和非变异)
    // gVCF 格式允许后续多样本联合分型

    tuple val(meta), path("*.tbi")          , optional:true, emit: tbi
    // tbi: Tabix index for the VCF file
    // Enables quick lookup of specific genomic regions
    // 加速特定基因组区域的查询
    // Generated automatically by GATK4 (bcftools)
    // 由 GATK4 自动生成

    tuple val(meta), path("*.realigned.bam"), optional:true, emit: bam
    // bam: Realigned BAM file (optional)
    // 重新比对的 BAM 文件 (可选)
    // Contains reads that were locally reassembled by HaplotypeCaller
    // 包含被 HaplotypeCaller 局部重组装过的 reads
    // Only generated if --bam-writer-type is specified
    // 仅当指定 --bam-writer-type 时生成
    // Useful for quality assessment and visualization in IGV
    // 用于质量评估和在 IGV 中可视化

    path "versions.yml"                     , emit: versions
    // versions.yml: Version tracking file
    // 版本跟踪文件，记录所有使用的工具版本
    // Used for reproducibility and troubleshooting
    // 用于结果重复性和问题排查

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Execution Control (执行控制)                                ║
    // ╚══════════════════════════════════════════════════════════════╝
    when:
    task.ext.when == null || task.ext.when
    // Conditional execution: Only run if task.ext.when is true
    // 条件执行: 仅当 task.ext.when 为 true 时运行
    // Allows configuration-based enable/disable of this process
    // 允许通过配置启用/禁用此进程

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Main Script Block (主脚本块)                                ║
    // ╚══════════════════════════════════════════════════════════════╝
    script:
    // Groovy script block: defines the actual command to execute
    // Groovy 脚本块: 定义要执行的实际命令

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Variable Setup (变量设置)                                   ║
    // ╚══════════════════════════════════════════════════════════════╝
    def args = task.ext.args ?: ''
    // Collect extra command-line arguments from configuration
    // 从配置中收集额外的命令行参数
    // Example: --standard-min-confidence-threshold-for-calling 10

    def prefix = task.ext.prefix ?: "${meta.id}"
    // Output file prefix (output filename prefix)
    // 输出文件前缀 (默认为样本 ID)
    // Output files: {prefix}.vcf.gz, {prefix}.vcf.gz.tbi

    def dbsnp_command = dbsnp ? "--dbsnp $dbsnp" : ""
    // Conditionally add dbsnp argument if file provided
    // 如果提供了 dbsnp 文件，添加 --dbsnp 参数
    // dbSNP helps annotate known variants

    def interval_command = intervals ? "--intervals $intervals" : ""
    // Conditionally add intervals argument if BED file provided
    // 如果提供了区间 BED 文件，添加 --intervals 参数
    // Limits analysis to specified genomic regions
    // 限制分析到指定的基因组区域

    def dragstr_command = dragstr_model ? "--dragstr-params-path $dragstr_model" : ""
    // Conditionally add DragStr model if provided
    // 如果提供了 DragStr 模型，添加参数
    // DragStr: Dynamic read pair graph model for short tandem repeats
    // DragStr: 用于短串联重复序列 (STR) 的动态读对图模型
    // Optional but improves accuracy in repeat regions
    // 可选但提高重复区域的准确性

    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""
    // Conditionally generate realigned BAM output
    // 如果用户请求，生成重新比对的 BAM 文件
    // Only if --bam-writer-type is in arguments
    // 仅当参数中包含 --bam-writer-type 时

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Memory Management (内存管理)                                ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Calculate available Java heap memory for GATK                ║
    // ║ 计算 GATK 可用的 Java 堆内存                                  ║
    // ║ Default: 3GB if not specified                               ║
    // ║ 默认: 如果未指定则为 3GB                                      ║
    // ╚══════════════════════════════════════════════════════════════╝
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
        // Use 80% of allocated memory (safety margin for JVM overhead)
        // 使用分配内存的 80% (为 JVM 开销留出安全边际)
    }

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Execute HaplotypeCaller (执行 HaplotypeCaller)             ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ Command structure:                                           ║
    // ║ 1. gatk: GATK wrapper script                                ║
    // ║ 2. HaplotypeCaller: Main tool name                          ║
    // ║ 3. Input flags: --input, --reference, etc.                  ║
    // ║ 4. Optional flags: --dbsnp, --intervals, etc.               ║
    // ║ 5. Resource flags: --native-pair-hmm-threads                ║
    // ╚══════════════════════════════════════════════════════════════╝
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        HaplotypeCaller \\
        --input $input \\
        // Input BAM/CRAM file (输入 BAM/CRAM 文件)
        --output ${prefix}.vcf.gz \\
        // Output gVCF filename (输出 gVCF 文件名)
        --reference $fasta \\
        // Reference genome sequence (参考基因组序列)
        --native-pair-hmm-threads ${task.cpus} \\
        // ╔══════════════════════════════════════════════════════════╗
        // ║  Multi-threading for Pair HMM Computation              ║
        // ║  配对 HMM 计算的多线程                                   ║
        // ╠══════════════════════════════════════════════════════════╣
        // ║ Pair HMM: Uses hidden Markov model to compute            ║
        // ║ likelihood of alignments vs. variant hypotheses          ║
        // ║                                                           ║
        // ║ 配对 HMM: 计算比对相对于变异假设的似然度                ║
        // ║ This is the most computationally expensive step          ║
        // ║ 这是计算上最昂贵的步骤                                    ║
        // ║ Using multiple threads significantly speeds up analysis  ║
        // ║ 使用多线程可显著加快分析速度                              ║
        // ╚══════════════════════════════════════════════════════════╝
        $dbsnp_command \\
        // Optional: Known variants for annotation
        // 可选: 用于标注的已知变异
        $interval_command \\
        // Optional: Limit analysis to specific genomic regions
        // 可选: 限制分析到特定基因组区域
        $dragstr_command \\
        // Optional: Model for short tandem repeat calling
        // 可选: 用于短串联重复序列的模型
        $bamout_command \\
        // Optional: Output realigned BAM file
        // 可选: 输出重新比对的 BAM 文件
        --tmp-dir . \\
        // Use current directory for temporary files
        // 使用当前目录存放临时文件 (避免填满系统 /tmp)
        $args
        // Additional arguments from configuration
        // 来自配置的额外参数

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Record Software Versions (记录软件版本)                    ║
    // ╚══════════════════════════════════════════════════════════════╝
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        // Extract GATK4 version number from version output
        // 从版本输出中提取 GATK4 版本号
        // Used for reproducibility tracking
        // 用于可重复性追踪
    END_VERSIONS
    """

    // ╔══════════════════════════════════════════════════════════════╗
    // ║  Stub Block: Test Mode (存根块: 测试模式)                   ║
    // ╠══════════════════════════════════════════════════════════════╣
    // ║ When using -stub option, creates dummy output files          ║
    // ║ 使用 -stub 选项时，创建虚拟输出文件                          ║
    // ║ Allows testing workflow without running actual analysis      ║
    // ║ 允许在不运行实际分析的情况下测试工作流                        ║
    // ║ Useful for validating workflow structure and parameters     ║
    // ║ 用于验证工作流结构和参数                                      ║
    // ╚══════════════════════════════════════════════════════════════╝
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""

    def stub_realigned_bam = bamout_command ? "touch ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""
    // Create empty realigned BAM file if requested
    // 如果请求，创建空的重新比对 BAM 文件

    """
    touch ${prefix}.vcf.gz
    // Create empty gVCF file (输出空的 gVCF 文件)

    touch ${prefix}.vcf.gz.tbi
    // Create empty index file (输出空的索引文件)

    ${stub_realigned_bam}
    // Create empty realigned BAM if requested (如果请求，创建空 BAM)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
