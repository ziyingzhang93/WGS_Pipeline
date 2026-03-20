# WGS_Pipeline

Annotated nf-core/sarek WGS pipeline guide with bilingual (EN/CN) comments for biomedical researchers.

## Online Guide

**[View the Interactive Pipeline Guide](https://ziyingzhang93.github.io/WGS_Pipeline/README_Guide.html)**

The guide includes a language toggle button — click "Show Chinese" to display Chinese translations alongside English content.

## Contents

- **README_Guide.html** — Interactive bilingual guide covering:
-   1. What is sarek?
    2.   2. WGS Pipeline Overview (9-step flowchart)
         3.   3. Code Structure Guide
              4.   4. Key Files Quick Reference
                   5.   5. Core Concepts (Channel, Process, Subworkflow, Meta Map)
                        6.   6. Biology Background Table
                             7.   7. Common Parameters
                                  8.   8. Three Variant Calling Modes (Germline / Somatic / Tumor-only)
                                       9.   9. Example Commands (WGS, WES, Somatic, Test data)
                                            10.   10. WGS vs WES — Key Differences
                                               
                                                  11. - **Annotated source code** — Key sarek pipeline files with bilingual comments:
                                                      -   - `workflows/sarek/main.nf` — Main pipeline orchestration
                                                          -   - `subworkflows/local/` — Preprocessing, alignment, variant calling
                                                              -   - `modules/nf-core/gatk4/haplotypecaller/` — GATK HaplotypeCaller module
                                                               
                                                                  - ## Quick Start
                                                               
                                                                  - ```bash
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
                                                                    - - [sarek Documentation](https://nf-co.re/sarek)
                                                                      - - [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us)
