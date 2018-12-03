# xms-rnaseq-nf
A kallisto quantification pipeline based on nextflow

## Requirements

- Unix-like operating system (Linux, macOS, etc)
- Java 8

## Quickstart

1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
```

2. clone `xms-rnaseq-nf` to loacl

```bash
git clone https://github.com/bioShaun/xms-rnaseq-nf.git
```

3. Launch the pipeline execution:

```bash

./nextflow run xms-rnaseq-nf \
    --reads fq_reads_pattern \
    --gene2tr gene_transcript_id_map_file \
    --outdir /path/to/outdir \
    --fasta /path/to/reference/fasta

# optional params
# --kallisto_index: using prepared kallisto index
# -profile slurm : launch jobs using slurm, default is on local machine

```

## Components

xms-rnaseq-nf uses the following software components and tools:

- [kallisto](https://pachterlab.github.io/kallisto/) 0.43.1