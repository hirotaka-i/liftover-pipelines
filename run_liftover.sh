#!/bin/bash
#
# Standalone wrapper script for running liftover-sumstats with Docker or Singularity
# Automatically detects which container runtime is available
#

set -e

# Default container image
CONTAINER_IMAGE="ghcr.io/hirotaka-i/liftover-pipelines:latest"
SINGULARITY_IMAGE=""

# Detect container runtime
if command -v docker &> /dev/null && docker info &> /dev/null 2>&1; then
    RUNTIME="docker"
elif command -v singularity &> /dev/null; then
    RUNTIME="singularity"
elif command -v apptainer &> /dev/null; then
    RUNTIME="apptainer"
else
    echo "Error: Neither Docker nor Singularity/Apptainer found" >&2
    exit 1
fi

# Get the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Parse arguments to find file paths for volume mounts
INPUT_FILE=""
OUTPUT_FILE=""
UNMATCHED_FILE=""
SOURCE_FASTA=""
TARGET_FASTA=""
CHAIN_FILE=""
MOUNT_DIRS=()

# Function to get absolute path
get_abs_path() {
    local path="$1"
    if [[ "$path" = /* ]]; then
        echo "$path"
    else
        echo "$(cd "$(dirname "$path")" 2>/dev/null && pwd)/$(basename "$path")"
    fi
}

# Function to add unique mount directory
add_mount_dir() {
    local dir="$1"
    # Check if already in array
    for existing in "${MOUNT_DIRS[@]}"; do
        if [[ "$existing" == "$dir" ]]; then
            return
        fi
    done
    MOUNT_DIRS+=("$dir")
}

# Parse arguments to extract file paths
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input|-i)
            INPUT_FILE="$(get_abs_path "$2")"
            add_mount_dir "$(dirname "$INPUT_FILE")"
            shift 2
            ;;
        --output|-o)
            OUTPUT_FILE="$(get_abs_path "$2")"
            add_mount_dir "$(dirname "$OUTPUT_FILE")"
            shift 2
            ;;
        --unmatched|-u)
            UNMATCHED_FILE="$(get_abs_path "$2")"
            add_mount_dir "$(dirname "$UNMATCHED_FILE")"
            shift 2
            ;;
        --source-fasta)
            SOURCE_FASTA="$(get_abs_path "$2")"
            add_mount_dir "$(dirname "$SOURCE_FASTA")"
            shift 2
            ;;
        --target-fasta)
            TARGET_FASTA="$(get_abs_path "$2")"
            add_mount_dir "$(dirname "$TARGET_FASTA")"
            shift 2
            ;;
        --chain-file)
            CHAIN_FILE="$(get_abs_path "$2")"
            add_mount_dir "$(dirname "$CHAIN_FILE")"
            shift 2
            ;;
        --singularity-image)
            SINGULARITY_IMAGE="$2"
            shift 2
            ;;
        --container-image)
            CONTAINER_IMAGE="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

# Build the container command based on runtime
if [[ "$RUNTIME" == "docker" ]]; then
    # Build Docker command with volume mounts
    DOCKER_CMD="docker run --rm --platform linux/amd64"
    
    for mount_dir in "${MOUNT_DIRS[@]}"; do
        DOCKER_CMD="$DOCKER_CMD -v ${mount_dir}:${mount_dir}:rw"
    done
    
    DOCKER_CMD="$DOCKER_CMD $CONTAINER_IMAGE liftover"
    
    # Add back the original arguments
    DOCKER_CMD="$DOCKER_CMD $@"
    
    echo "Running with Docker..." >&2
    echo "Command: $DOCKER_CMD" >&2
    eval "$DOCKER_CMD"
    
else
    # Singularity/Apptainer command
    if [[ -z "$SINGULARITY_IMAGE" ]]; then
        # Try to find a .sif file in the current directory
        if [[ -f "liftover-sumstats.sif" ]]; then
            SINGULARITY_IMAGE="liftover-sumstats.sif"
        elif [[ -f "${SCRIPT_DIR}/liftover-sumstats.sif" ]]; then
            SINGULARITY_IMAGE="${SCRIPT_DIR}/liftover-sumstats.sif"
        else
            SINGULARITY_IMAGE="docker://${CONTAINER_IMAGE}"
        fi
    fi
    
    # Build Singularity bind arguments
    BIND_ARGS=""
    for mount_dir in "${MOUNT_DIRS[@]}"; do
        BIND_ARGS="$BIND_ARGS -B ${mount_dir}:${mount_dir}:rw"
    done
    
    echo "Running with ${RUNTIME^}..." >&2
    echo "Image: $SINGULARITY_IMAGE" >&2
    
    $RUNTIME exec $BIND_ARGS "$SINGULARITY_IMAGE" liftover "$@"
fi
