#!/bin/bash
#
# Example usage with the wrapper script
# This shows how users would run the tool
#

# Set your reference directory (adjust this path!)
REFERENCE_DIR="${HOME}/references"

# Example 1: Using the wrapper script (easiest - auto-detects Docker/Singularity)
echo "=== Example 1: Using wrapper script ==="
./run_liftover.sh \
  --input example/sumstats_hg19.txt \
  --output output/sumstats_hg38.txt \
  --unmatched output/unmatched.txt \
  --chr-col CHR \
  --pos-col POS \
  --ea-col A1 \
  --ref-col A2 \
  --effect-col Z \
  --eaf-col A1_FREQ \
  --add-chr-prefix \
  --source-fasta ${REFERENCE_DIR}/hg19.fa.gz \
  --target-fasta ${REFERENCE_DIR}/hg38.fa.gz \
  --chain-file ${REFERENCE_DIR}/hg19ToHg38.over.chain.gz

echo -e "\n=== Example 2: Direct Docker usage ==="
# Example 2: Direct Docker usage
# Note: All paths must be mapped through volume mounts!
docker run --rm \
  -v $(pwd)/example:/example:rw \
  -v $(pwd)/output:/output:rw \
  -v ${REFERENCE_DIR}:${REFERENCE_DIR}:ro \
  ghcr.io/hirotaka-i/liftover-pipelines:latest \
  --input /example/sumstats_hg19.txt \
  --output /output/sumstats_hg38.txt \
  --unmatched /output/unmatched.txt \
  --chr-col CHR \
  --pos-col POS \
  --ea-col A1 \
  --ref-col A2 \
  --effect-col Z \
  --eaf-col A1_FREQ \
  --add-chr-prefix \
  --source-fasta ${REFERENCE_DIR}/hg19.fa.gz \
  --target-fasta ${REFERENCE_DIR}/hg38.fa.gz \
  --chain-file ${REFERENCE_DIR}/hg19ToHg38.over.chain.gz

echo -e "\n=== Example 3: Singularity usage (if available) ==="
# Example 3: Singularity usage (for HPC systems)
if command -v singularity &> /dev/null; then
    singularity exec \
      -B $(pwd)/example:/example:rw \
      -B $(pwd)/output:/output:rw \
      -B ${REFERENCE_DIR}:${REFERENCE_DIR}:ro \
      liftover-sumstats.sif \
      liftover \
      --input /example/sumstats_hg19.txt \
      --output /output/sumstats_hg38.txt \
      --unmatched /output/unmatched.txt \
      --chr-col CHR \
      --pos-col POS \
      --ea-col A1 \
      --ref-col A2 \
      --effect-col Z \
      --eaf-col A1_FREQ \
      --add-chr-prefix \
      --source-fasta ${REFERENCE_DIR}/hg19.fa.gz \
      --target-fasta ${REFERENCE_DIR}/hg38.fa.gz \
      --chain-file ${REFERENCE_DIR}/hg19ToHg38.over.chain.gz
else
    echo "Singularity not available, skipping..."
fi
