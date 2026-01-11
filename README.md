# GWAS Summary Statistics Liftover

A minimal Docker-based tool to convert GWAS summary statistics between genome builds (e.g., hg19 ↔ hg38) with proper effect allele and allele frequency adjustments.

## Features

- Convert summary statistics between genome builds using bcftools +liftover
- Automatically flip effect sizes and allele frequencies when alleles swap
- Two versions available:
  - **Simple**: Direct liftover without pre-normalization
  - **Standard**: With REF allele normalization before liftover
- Dockerized for easy setup and reproducibility

## Quick Start

For the fastest start, see [QUICKSTART.md](QUICKSTART.md).

### Prerequisites

- Docker installed on your system
- Reference genome files (FASTA) for source and target builds
- Chain file for coordinate conversion (e.g., hg19ToHg38.over.chain.gz)

### 1. Build the Docker Image

```bash
docker build -t liftover-sumstats .
```

### 2. Prepare Your Data

Create directories for your data and reference files:

```bash
mkdir -p data references
```

### 3. Run Liftover

#### Using Docker directly:

```bash
docker run --rm \
  -v $(pwd)/data:/workspace/data \
  -v $(pwd)/references:/references \
  liftover-sumstats \
  python3 /usr/local/bin/liftover_sumstats_simple.py \
    --input /workspace/data/input.txt \
    --output /workspace/data/output.hg38.txt \
    --unmatched /workspace/data/output.unmatched.txt \
    --chr-col CHR \
    --pos-col POS \
    --ea-col EA \
    --ref-col REF \
    --source-fasta /references/hg19.fa \
    --target-fasta /references/hg38.fa \
    --chain-file /references/hg19ToHg38.over.chain.gz \
    --effect-col BETA \
    --eaf-col EAF
```

#### Using the wrapper script (recommended):

```bash
./liftover_sumstats_wrapper_simple.sh \
  -i data/input.txt \
  -o data/output.hg38.txt \
  -u data/output.unmatched.txt \
  --chr-col CHR \
  --pos-col POS \
  --ea-col EA \
  --ref-col REF \
  --source-fasta references/hg19.fa \
  --target-fasta references/hg38.fa \
  --chain-file references/hg19ToHg38.over.chain.gz \
  --effect-col BETA \
  --eaf-col EAF
```

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

### Reference Genomes

Download from UCSC or other sources:

```bash
# hg19 (GRCh37)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa

# hg38 (GRCh38)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

### Chain Files

Download from UCSC:

```bash
# hg19 to hg38
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# hg38 to hg19
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
```

## Input File Format

Your summary statistics file should be a tab or space-delimited text file with headers. Example:

```
CHR  POS       EA  REF  BETA      SE        P         EAF
1    12345     A   G    0.123     0.045     0.001     0.35
1    67890     C   T    -0.087    0.038     0.023     0.62
```

## Output

The tool produces:
1. **Lifted file**: Summary statistics with updated coordinates
2. **Unmatched file**: Variants that couldn't be lifted (no match in chain file or different alleles)

## Project Structure

```
.
├── Dockerfile                              # Minimal Docker image for liftover
├── docker-compose.yml                      # Docker compose configuration
├── README.md                               # This file
├── liftover_sumstats.py                   # Standard version with normalization
├── liftover_sumstats_simple.py            # Simple version without normalization
├── liftover_sumstats_wrapper.sh           # Wrapper for standard version
├── liftover_sumstats_wrapper_simple.sh    # Wrapper for simple version
├── data/                                   # Your data files (create this)
├── references/                             # Reference genomes (create this)
├── examples/                               # Example/test files
└── docker/                                 # Full pipeline Dockerfile
    └── Dockerfile.ubuntu22
```

## Troubleshooting

### Issue: Chromosome naming mismatch

If your source FASTA uses "chr1, chr2, ..." but your summary statistics use "1, 2, ...", add the `--add-chr-prefix` flag.

### Issue: Many unmatched variants

- Check that your input coordinates match the source genome build
- Verify that chromosome naming is consistent
- Ensure your chain file matches your source→target conversion

### Issue: Docker permission errors

Make sure your data and references directories have proper permissions:

```bash
chmod -R 755 data references
```

## Examples

See the `examples/` directory for sample input files and test cases.

## License

This tool uses:
- BCFtools (MIT/Expat License)
- UCSC liftOver tool (free for academic use)
- Python pandas, numpy (BSD License)

## Citation

If you use this tool, please cite:
- BCFtools: Danecek P, et al. (2021) Twelve years of SAMtools and BCFtools. GigaScience, 10(2).
- Liftover plugin: Giulio Genovese (freeseek/score repository)

# MEMO
```
Flip mechanism

The assumption is that they are all forward strand. 
They can flip strand while lifting over but
1. first align to hg19 (no strand flip)
2. if hg38 forward is the flip of hg19, then flip

So flipped strand as input will be never fixed and create problematic output. 


How normalization work
if C T is the truth >
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

How bcftools liftover work
REF_IN	ALT_IN	OPERATION	REF_OUT	ALT_OUT	Notes
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