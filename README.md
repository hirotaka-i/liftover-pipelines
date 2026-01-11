# GWAS Summary Statistics Liftover

A minimal Docker-based tool to convert GWAS summary statistics between genome builds (e.g., hg19 ↔ hg38) with proper effect allele and allele frequency adjustments.

## Features

- Convert summary statistics between genome builds using bcftools +liftover
- Automatically flip effect sizes and allele frequencies when alleles swap
- Standalone Docker/Singularity container for easy distribution
- Simple direct liftover without pre-normalization
- Dockerized for easy setup and reproducibility

## Quick Start
You need to moount the following folders into the container:
* data folder: the folder containing your summary statistics file  
* refs folder: the folder containing reference files (FASTA, chain files)


### Option 1: Docker directly (No download needed!)

```bash
# Pull the pre-built image (works on Mac with Rosetta 2)
docker pull --platform linux/amd64 ghcr.io/hirotaka-i/liftover-pipelines:latest

# Run with your data
docker run --rm --platform linux/amd64 \
  -v /path/to/your/data:/data:rw \
  -v /path/to/references:/refs:ro \
  ghcr.io/hirotaka-i/liftover-pipelines:latest \
  --input /data/sumstats_hg19.txt \
  --output /data/sumstats_hg38.txt \
  --unmatched /data/unmatched.txt \
  --chr-col CHR \
  --pos-col POS \
  --ea-col A1 \
  --ref-col A2 \
  --effect-col Z \
  --eaf-col A1_FREQ \
  --add-chr-prefix \
  --source-fasta /refs/hg19.fa.gz \
  --target-fasta /refs/hg38.fa.gz \
  --chain-file /refs/hg19ToHg38.over.chain.gz
```


### Option 2: Singularity/Apptainer (Manual method for HPC)

If you prefer to use Singularity directly instead of the wrapper script:

```bash
# Convert to Singularity SIF
singularity build liftover-sumstats.sif docker://ghcr.io/hirotaka-i/liftover-pipelines:latest

# Run with Singularity
singularity exec \
  -B /path/to/data:/data:rw \
  -B /path/to/refs:/refs:ro \
  liftover-sumstats.sif \
  liftover \
  --input /data/sumstats_hg19.txt \
  --output /data/sumstats_hg38.txt \
  --unmatched /data/unmatched.txt \
  --chr-col CHR \
  --pos-col POS \
  --ea-col A1 \
  --ref-col A2 \
  --effect-col Z \
  --eaf-col A1_FREQ \
  --add-chr-prefix \
  --source-fasta /refs/hg19.fa.gz \
  --target-fasta /refs/hg38.fa.gz \
  --chain-file /refs/hg19ToHg38.over.chain.gz
```

**Important:** When using Docker directly, all file paths inside the container must match your volume mount paths!

## Required Arguments

- `--input`: Input summary statistics file (txt or txt.gz)
- `--output`: Output lifted summary statistics file
- `--unmatched`: Output file for variants that couldn't be lifted
- `--chr-col`: Chromosome column name in your input file
- `--pos-col`: Position column name
- `--ea-col`: Effect allele column name
- `--ref-col`: Reference allele column name
- `--source-fasta`: Source genome reference FASTA file
- `--target-fasta`: Target genome reference FASTA file
- `--chain-file`: Chain file for liftover

## Optional Arguments

- `--effect-col`: Effect column name to flip when alleles swap (e.g., BETA, Z). Can specify multiple times
- `--eaf-col`: Effect allele frequency column to flip (e.g., EAF). Can specify multiple times
- `--add-chr-prefix`: Add 'chr' prefix to chromosomes if your source FASTA uses chr1, chr2, etc.
- `--keep-temp`: Keep temporary files for debugging

## Getting Reference Files

You'll need to download these reference files once and reuse them:

### Chain Files
```bash
# hg19 to hg38
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# hg38 to hg19
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
```

### Reference Genomes
```bash
# hg19 (GRCh37) - can keep gzipped
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# hg38 (GRCh38) - can keep gzipped
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

## Input File Format

Your summary statistics file should be a tab or space-delimited text file with headers. Example:

```
CHR  POS       EA  REF  BETA      SE        P         EAF
1    12345     A   G    0.123     0.045     0.001     0.35
1    67890     C   T    -0.087    0.038     0.023     0.62
...
```

Gzipped files (.txt.gz) are also supported.

## Output

The tool produces:
1. **Lifted file**: Summary statistics with updated coordinates in the target genome build
2. **Unmatched file**: Variants that couldn't be lifted (no match in chain file or different alleles)

## Real-World Example

```bash
# Using the wrapper script with real paths
./run_liftover.sh \
  --input $HOME/gwas_data/my_gwas_hg19.txt.gz \
  --output $HOME/gwas_data/my_gwas_hg38.txt.gz \
  --unmatched $HOME/gwas_data/my_gwas_unmatched.txt \
  --chr-col CHR \
  --pos-col BP \
  --ea-col A1 \
  --ref-col A2 \
  --effect-col BETA \
  --eaf-col FRQ \
  --source-fasta $HOME/references/hg19.fa.gz \
  --target-fasta $HOME/references/hg38.fa.gz \
  --chain-file $HOME/references/hg19ToHg38.over.chain.gz
```

## Troubleshooting

### Chromosome naming mismatch
If your source FASTA uses "chr1, chr2, ..." but your summary statistics use "1, 2, ...", add the `--add-chr-prefix` flag.

### Many unmatched variants
- Check that your input coordinates match the source genome build
- Verify that chromosome naming is consistent
- Ensure your chain file matches your source→target conversion

### Docker permission errors
Run with your user ID:
```bash
docker run --rm --platform linux/amd64 --user $(id -u):$(id -g) \
  -v /path/to/data:/data:rw \
  ghcr.io/hirotaka-i/liftover-pipelines:latest ...
```

### File not found errors
Make sure:
1. All file paths inside the container match your volume mounts
2. You're using absolute paths for volume mounts
3. Files have read permissions and output directories have write permissions

### Slow on Mac (Apple Silicon)
The image runs via Rosetta 2 emulation on Apple Silicon Macs, which should be fine for most use cases. The `--platform linux/amd64` flag is required for compatibility.

### Out of memory
For very large files, increase Docker's memory limit in Docker Desktop settings.

## Citation

If you use this tool, please cite:
- BCFtools: Danecek P, et al. (2021) Twelve years of SAMtools and BCFtools. GigaScience, 10(2).
- Liftover plugin: Giulio Genovese (freeseek/score repository)

## License

This tool uses:
- BCFtools (MIT/Expat License)
- BCFtools/liftover pluging[MIT License]: 
  - Genovese G., McCarroll S. et al. BCFtools/liftover: an accurate and comprehensive tool to convert genetic variants across genome assemblies. Bioinformatics 40, Issue 2 (2024). [PMID: 38261650] [DOI: 10.1093/bioinformatics/btae038]
- Python pandas, numpy (BSD License)

---

## Important Note on Strand Assumptions

For the input, the assumption is that they are all forward stranded. The pipeline can flip strand while lifting over but it only happens the flipping between different genome builds. If the input file is not aligned to the forward strand of the source genome build, the output would contain a lot of misaligned variants.


```
# How normalization works if C T is the truth >
REF_IN	ALT_IN	OPERATION	REF_OUT	ALT_OUT	Notes
C	T	No change	C	T	REF matches FASTA
C	A	No change	C	A	REF matches FASTA
C	G	No change	C	G	REF matches FASTA
T	C	Allele swap	C	T	ALT matches FASTA, swap to fix
A	C	Allele swap	C	A	ALT matches FASTA, swap to fix
G	C	Allele swap	C	G	ALT matches FASTA, swap to fix
A	T	Ref added	C	T	REF_IN replaced with C, ALT unchanged
A	G	Ref added	C	G	REF_IN replaced with C, ALT unchanged
T	A	Ref added	C	A	REF_IN replaced with C, ALT unchanged
T	G	Ref added	C	G	REF_IN replaced with C, ALT unchanged
G	A	Ref added	C	A	REF_IN replaced with C, ALT unchanged
G	T	Ref added	C	T	REF_IN replaced with C, ALT unchanged

# How bcftools +liftover works
# REF_IN	ALT_IN	OPERATION	REF_OUT	ALT_OUT	Notes
C	T	No change	C	T	REF matches FASTA
C	A	No change	C	A	REF matches FASTA
C	G	No change	C	G	REF matches FASTA
T	C	Allele swap	C	T	ALT matches FASTA, swap to fix
A	C	Allele swap	C	A	ALT matches FASTA, swap to fix
G	C	Allele swap	C	G	ALT matches FASTA, swap to fix
A	T	Ref added	C	A,T	REF_IN added, A,T to ALT
A	G	Ref added	C	G,A	REF_IN replaced with C, ALT unchanged
T	A	Ref added	C	T,A	REF_IN added, A,T to ALT
T	G	Ref added	C	G,T	REF_IN replaced with C, ALT unchanged
G	A	Ref added	C	A,G	REF_IN replaced with C, ALT unchanged
G	T	Ref added	C	T,G	REF_IN replaced with C, ALT unchanged
```