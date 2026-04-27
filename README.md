# winnowmap
Winnowmap is a long-read mapping algorithm optimized for mapping ONT and PacBio reads to repetitive reference sequences. Winnowmap development began on top of minimap2 codebase, and since then we have incorporated the following two ideas to improve mapping accuracy within repeats.
# winnowmap-ont-pipeline

A production-ready Bash pipeline for aligning Oxford Nanopore Technology (ONT) reads using **Winnowmap2**, with full **methylation tag (MM/ML) preservation** and support for multi-chunk uBAM input from P2 Solo and other ONT sequencers.

---

## Features

- Automated **meryl k-mer counting** and repetitive k-mer extraction
- **Single uBAM** or **directory of uBAM chunks** as input (auto-merge)
- Alignment with **Winnowmap2** in SV-aware ONT mode
- **MM/ML methylation tag transplant** from uBAM into aligned BAM via pysam
- Per-step timing and colored log output
- Optional intermediate file cleanup
- Flagstat QC and chromosome coverage summary saved automatically

---

## Dependencies

| Tool | Version tested | Install |
|---|---|---|
| [Winnowmap2](https://github.com/marbl/Winnowmap) | 2.03 | `conda install -c bioconda winnowmap` |
| [meryl](https://github.com/marbl/meryl) | r977 | `conda install -c bioconda meryl` |
| [samtools](https://www.htslib.org/) | ≥ 1.17 | `conda install -c bioconda samtools` |
| [pysam](https://pysam.readthedocs.io/) | ≥ 0.22 | `pip install pysam` |
| Python | ≥ 3.8 | — |

---

## Installation

```bash
git clone https://github.com/<your-username>/winnowmap-ont-pipeline.git
cd winnowmap-ont-pipeline
chmod +x winnowmap_align.sh
```

---

## Usage

```
bash winnowmap_align.sh -i <ubam_or_dir> -r <reference.fa> -o <output_dir> [options]
```

### Required arguments

| Flag | Description |
|---|---|
| `-i` | Input uBAM file **or** directory of uBAM chunks |
| `-r` | Reference FASTA (e.g. `hg38.fa`, must be `samtools faidx` indexed) |
| `-o` | Output directory (created if it does not exist) |

### Optional arguments

| Flag | Description | Default |
|---|---|---|
| `-s` | Sample name | Derived from uBAM filename |
| `-t` | Threads | `16` |
| `-k` | k-mer size for Winnowmap2 and meryl | `15` |
| `-w` | Meryl database directory | `<ref_dir>/<ref_base>.meryl` |
| `-f` | High-frequency k-mer file | `<ref_dir>/<ref_base>_repetitive_k<k>.txt` |
| `-m` | Skip meryl steps if k-mer file already exists | `false` |
| `-c` | Remove intermediate files after completion | `false` |
| `-h` | Show help message | — |

---

## Examples

### Single uBAM file

```bash
bash winnowmap_align.sh \
  -i /data/ubam/sample.bam \
  -r /data/reference/hg38.fa \
  -o /data/output \
  -s MySample \
  -t 16 \
  -c
```

### Directory of uBAM chunks (typical ONT run output)

```bash
bash winnowmap_align.sh \
  -i /data/ubam_chunks/ \
  -r /data/reference/hg38.fa \
  -o /data/output \
  -s MySample \
  -t 16 \
  -c
```

### Reuse existing k-mer file (skip meryl on repeated runs)

```bash
bash winnowmap_align.sh \
  -i /data/ubam/sample.bam \
  -r /data/reference/hg38.fa \
  -o /data/output \
  -f /data/reference/hg38_repetitive_k15.txt \
  -m \
  -c
```

---

## Pipeline steps

```
Input uBAM / chunks
       │
       ▼
[Auto-merge chunks]          ← if directory input with multiple BAMs
       │
       ▼
[meryl k-mer count]          ← skipped if k-mer file already exists (-m)
       │
       ▼
[meryl extract repetitive k-mers]
       │
       ▼
[samtools fastq -T MM,ML]    ← preserves methylation tags in FASTQ
       │
       ▼
[Winnowmap2 -ax map-ont]     ← SV-aware ONT alignment
       │
       ▼
[samtools sort + index]
       │
       ▼
[pysam MM/ML tag transplant] ← injects methylation tags back from uBAM
       │
       ▼
[samtools sort + index]
       │
       ▼
Final: <sample>_aligned_meth_sorted.bam
       + <sample>_flagstat.txt
```

---

## Output files

| File | Description |
|---|---|
| `<sample>_aligned_meth_sorted.bam` | Final coordinate-sorted BAM with MM/ML tags |
| `<sample>_aligned_meth_sorted.bam.bai` | BAM index |
| `<sample>_flagstat.txt` | samtools flagstat QC summary |

When `-c` is used, all intermediate files (FASTQ, SAM, unsorted BAMs) are removed automatically.

---

## Why methylation tags need transplanting

Winnowmap2 aligns reads from FASTQ input but does not propagate BAM tags such as `MM` (base modification) and `ML` (modification likelihood) into its SAM output. These tags, written by Dorado during basecalling, are present in the uBAM but lost during the FASTQ conversion step. This pipeline uses pysam to look up each read by name in the original uBAM and inject the tags back into the aligned BAM, ensuring all downstream methylation tools (Modkit, DSS, bismark) receive complete data.

---

## Notes on meryl k-mer files

The meryl database and repetitive k-mer file only need to be generated **once per reference genome**. For all subsequent samples aligned to the same reference, use `-m -f` to skip this step and save significant time.

| Reference | k=15 k-mer count | Build time (32 threads) |
|---|---|---|
| hg38 (unmasked) | ~66,000 | ~20 min |
| hg38 (masked) | ~10,000–30,000 | ~15 min |

---

## Downstream analysis

The output BAM is compatible with:

| Tool | Purpose |
|---|---|
| [Clair3](https://github.com/HKU-BAL/Clair3) | SNV and indel calling |
| [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) | Structural variant calling |
| [Severus](https://github.com/KolmogorovLab/Severus) | SV calling with phasing |
| [Modkit](https://github.com/nanoporetech/modkit) | CpG methylation pileup |
| [DSS](https://bioconductor.org/packages/DSS/) | Differential methylation |
| [WhatsHap](https://whatshap.readthedocs.io/) | Read phasing |

---

## Citation

If you use Winnowmap2 in your work, please cite:

> Jain et al. (2022). Long-read mapping to repetitive regions with Winnowmap2. *Nature Methods*, 19, 705–710. https://doi.org/10.1038/s41592-022-01457-8

---

## License

MIT
