#!/bin/bash
#
# Simple setup script to help users get started
#

echo "=== GWAS Summary Statistics Liftover Tool Setup ==="
echo ""

# Check for Docker or Singularity
if command -v docker &> /dev/null && docker info &> /dev/null 2>&1; then
    echo "✓ Docker detected"
    RUNTIME="docker"
elif command -v singularity &> /dev/null; then
    echo "✓ Singularity detected"
    RUNTIME="singularity"
elif command -v apptainer &> /dev/null; then
    echo "✓ Apptainer detected"
    RUNTIME="apptainer"
else
    echo "✗ No container runtime found. Please install Docker or Singularity."
    exit 1
fi

echo ""
echo "To use this tool, you need:"
echo "1. Reference genome files (hg19.fa.gz, hg38.fa.gz)"
echo "2. Chain file (hg19ToHg38.over.chain.gz or hg38ToHg19.over.chain.gz)"
echo ""
echo "Download reference files:"
echo "  mkdir -p references"
echo "  wget -P references http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
echo "  wget -P references http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
echo "  wget -P references http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
echo ""

if [[ "$RUNTIME" == "docker" ]]; then
    echo "Getting the Docker image..."
    
    # Check if we should pull from GitHub or build locally
    if [[ -f "Dockerfile" ]]; then
        echo "Dockerfile found. Building locally..."
        docker build -t liftover-sumstats:latest .
    else
        echo "Pulling from GitHub Container Registry..."
        docker pull ghcr.io/hirotaka-i/liftover_sumstats:latest
    fi
else
    echo "For Singularity/Apptainer, you can:"
    echo "1. Pull from GitHub: singularity build liftover-sumstats.sif docker://ghcr.io/hirotaka-i/liftover_sumstats:latest"
    echo "2. Or build from local Docker: singularity build liftover-sumstats.sif docker-daemon://liftover-sumstats:latest"
fi

echo ""
echo "Example usage:"
echo "./run_liftover.sh \\"
echo "  --input example/sumstats_hg19.txt \\"
echo "  --output sumstats_hg38.txt \\"
echo "  --unmatched unmatched.txt \\"
echo "  --chr-col CHR \\"
echo "  --pos-col POS \\"
echo "  --ea-col A1 \\"
echo "  --ref-col A2 \\"
echo "  --effect-col Z \\"
echo "  --eaf-col A1_FREQ \\"
echo "  --source-fasta references/hg19.fa.gz \\"
echo "  --target-fasta references/hg38.fa.gz \\"
echo "  --chain-file references/hg19ToHg38.over.chain.gz"
echo ""
echo "See README.md for more details!"
