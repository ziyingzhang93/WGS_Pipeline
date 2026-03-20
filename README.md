# WGS_Pipeline

Annotated nf-core/sarek WGS/WES pipeline guide with bilingual (EN/CN) comments for biomedical researchers.

## Online Guide

**[View the Interactive Pipeline Guide](https://ziyingzhang93.github.io/WGS_Pipeline/README_Guide.html)**

The guide includes a language toggle button — click "Show Chinese" to display Chinese translations alongside English content.

## Quick Start

```bash
# Germline WGS analysis
nextflow run nf-core/sarek \
  --input samplesheet.csv \
  --genome GATK.GRCh38 \
  --tools haplotypecaller \
  -profile docker

# WES analysis (add --wes and --intervals)
nextflow run nf-core/sarek \
  --input samplesheet.csv \
  --genome GATK.GRCh38 \
  --wes --intervals targets.bed \
  --tools haplotypecaller \
  -profile docker
```

## References

- [nf-core/sarek GitHub](https://github.com/nf-core/sarek)
- [sarek Documentation](https://nf-co.re/sarek)
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us)
