#!/usr/bin/env python3
"""
Liftover summary statistics between genome builds using bcftools +liftover
Simple version without pre-normalization
"""

import argparse
import gzip
import sys
import subprocess
import tempfile
import os
import pandas as pd
from pathlib import Path
import gc
import logging


def open_file(filepath):
    """Open file, handling gzipped files"""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def write_file(filepath):
    """Write file, handling gzipped files"""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'wt')
    return open(filepath, 'w')


def sumstats_to_vcf(sumstats_file, vcf_file, chr_col, pos_col, ea_col, non_ea_col, 
                    add_chr_prefix=False, source_fasta=None):
    """
    Convert summary statistics to VCF format
    
    Args:
        sumstats_file: Input summary statistics file
        vcf_file: Output VCF file path
        chr_col: Chromosome column name
        pos_col: Position column name
        ea_col: Effect allele column name
        non_ea_col: Non-effect allele column name
        add_chr_prefix: Whether to add 'chr' prefix to chromosome names
        source_fasta: Source reference fasta file (for contig definitions only)
    """
    logging.info(f"Reading summary statistics from {sumstats_file}...")
    
    # Read summary stats
    df = pd.read_csv(sumstats_file, sep='\t', low_memory=False)
    
    # Check required columns
    required_cols = [chr_col, pos_col, ea_col, non_ea_col]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    logging.info(f"Loaded {len(df)} variants")
    
    # Clean chromosome names - handle numeric chromosomes (may be float like 1.0)
    # Convert to string first to avoid dtype issues
    df[chr_col] = df[chr_col].astype(str).str.strip()
    
    # Try to convert numeric-looking chromosomes (remove .0)
    try:
        numeric_df = pd.to_numeric(df[chr_col], errors='coerce')
        numeric_mask = numeric_df.notna()
        # Convert to int to remove .0, then back to string
        df.loc[numeric_mask, chr_col] = numeric_df[numeric_mask].astype(int).astype(str)
    except:
        pass
    
    # Add chr prefix if requested (to match reference fasta)
    if add_chr_prefix:
        # Only add to numeric chromosomes and X, Y, MT, M
        simple_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT', 'M']
        mask = df[chr_col].isin(simple_chroms)
        df.loc[mask, chr_col] = 'chr' + df.loc[mask, chr_col]
        logging.info(f"Added 'chr' prefix to {mask.sum()} chromosomes")
    
    # Debug: show unique chromosome values
    unique_chrs = df[chr_col].unique()[:20]  # First 20 unique values
    logging.debug(f"Sample chromosome values after cleaning: {unique_chrs}")
    
    # Filter to valid chromosomes (accept both with and without chr prefix)
    valid_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT', 'M']
    valid_chroms_with_chr = ['chr' + c for c in valid_chroms]
    all_valid = valid_chroms + valid_chroms_with_chr
    df = df[df[chr_col].isin(all_valid)].copy()
    logging.info(f"Kept {len(df)} variants on valid chromosomes")
    
    # Sort by chromosome and position (AFTER adding chr prefix)
    df[pos_col] = pd.to_numeric(df[pos_col], errors='coerce')
    df = df.dropna(subset=[pos_col])
    df[pos_col] = df[pos_col].astype(int)
    
    # Create sort key for chromosomes
    # Handle both with and without chr prefix
    def chr_sort_key(chr_str):
        chr_clean = chr_str.replace('chr', '')
        if chr_clean.isdigit():
            return int(chr_clean)
        elif chr_clean == 'X':
            return 23
        elif chr_clean == 'Y':
            return 24
        elif chr_clean in ['MT', 'M']:
            return 25
        else:
            return 99
    
    df['_chr_sort'] = df[chr_col].apply(chr_sort_key)
    df = df.sort_values(['_chr_sort', pos_col]).drop('_chr_sort', axis=1)
    logging.info(f"Sorted {len(df)} variants by chromosome and position")
    
    # Write VCF - write uncompressed first, then compress with bcftools
    logging.info(f"Writing VCF to {vcf_file}...")
    
    # Write to uncompressed file first
    temp_vcf = vcf_file.replace('.vcf.gz', '.vcf')
    with open(temp_vcf, 'w') as f:
        # Write header
        f.write("##fileformat=VCFv4.2\n")
        
        # Add contig definitions ONLY for chromosomes present in data (memory efficient)
        if source_fasta:
            logging.info(f"Adding contig definitions for present chromosomes...")
            try:
                # Get unique chromosomes in our data
                unique_chroms = set(df[chr_col].unique())
                
                # Read fasta index to get contig lengths ONLY for present chromosomes
                fai_file = source_fasta + '.fai'
                if os.path.exists(fai_file):
                    with open(fai_file, 'r') as fai:
                        for line in fai:
                            parts = line.strip().split('\t')
                            if len(parts) >= 2:
                                contig = parts[0]
                                # Only write contigs for chromosomes in our data
                                if contig in unique_chroms:
                                    length = parts[1]
                                    f.write(f"##contig=<ID={contig},length={length}>\n")
            except Exception as e:
                logging.warning(f"Could not read contig definitions: {e}")
        
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # Build VCF columns using vectorized pandas operations (much faster than iterrows)
        vcf_df = pd.DataFrame()
        vcf_df['CHROM'] = df[chr_col]
        vcf_df['POS'] = df[pos_col]
        
        # Create variant IDs in chr:pos:ref:alt format
        vcf_df['ID'] = (df[chr_col].astype(str) + ':' + df[pos_col].astype(str) + ':' + 
                       df[non_ea_col].astype(str).str.upper() + ':' + df[ea_col].astype(str).str.upper())
        
        vcf_df['REF'] = df[non_ea_col].astype(str).str.upper()
        vcf_df['ALT'] = df[ea_col].astype(str).str.upper()
        vcf_df['QUAL'] = '.'
        vcf_df['FILTER'] = '.'
        vcf_df['INFO'] = '.'
        
        logging.info(f"Writing {len(vcf_df)} variants...")
        
        # Write all at once (much faster than iterrows)
        vcf_df.to_csv(f, sep='\t', index=False, header=False)
    
    logging.info(f"Wrote {len(df)} variants to VCF")
    
    # Compress with bcftools (creates proper BGZF compression)
    logging.info(f"Compressing with bcftools...")
    subprocess.run(['bcftools', 'view', '-Oz', '-o', vcf_file, temp_vcf], check=True)
    
    # Remove uncompressed file
    os.remove(temp_vcf)
    
    # Index the VCF file
    logging.info(f"Indexing {vcf_file}...")
    subprocess.run(['bcftools', 'index', '-t', vcf_file], check=True)
    
    return len(df)


def run_liftover(input_vcf, output_vcf, chain_file, source_fasta, target_fasta):
    """
    Run bcftools +liftover
    
    Args:
        input_vcf: Input VCF file
        output_vcf: Output lifted VCF file
        chain_file: Chain file for liftover
        source_fasta: Reference fasta for source build
        target_fasta: Reference fasta for target build
    """
    logging.info(f"Running bcftools +liftover...")
    
    # Create rejected file path
    rejected_vcf = output_vcf.replace('.vcf', '.rejected.vcf')
    
    cmd = [
        'bcftools', '+liftover',
        '--no-version',
        '-Oz',
        '-o', output_vcf,
        input_vcf,
        '--',
        '-s', source_fasta,
        '-f', target_fasta,
        '-c', chain_file,
        '--reject', rejected_vcf,
        '-O', 'z'
    ]
    
    logging.info(f"Command: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Print stdout and stderr for debugging
    if result.stdout:
        logging.debug(f"bcftools stdout: {result.stdout}")
    if result.stderr:
        logging.info(f"bcftools stderr: {result.stderr}")
    
    if result.returncode != 0:
        logging.error(f"Error running bcftools +liftover:")
        raise RuntimeError(f"bcftools +liftover failed with exit code {result.returncode}")
    
    # Check if output file was created
    if not os.path.exists(output_vcf):
        raise RuntimeError(f"bcftools +liftover did not create output file: {output_vcf}")
    
    if os.path.getsize(output_vcf) == 0:
        logging.warning(f"Output VCF file is empty: {output_vcf}")
    
    # Sort the lifted VCF (liftover can produce unsorted output)
    logging.info(f"Sorting {output_vcf}...")
    sorted_vcf = output_vcf + '.sorting.tmp'
    subprocess.run(['bcftools', 'sort', '-Oz', '-o', sorted_vcf, output_vcf], check=True)
    # Atomic rename - overwrites output_vcf safely
    os.rename(sorted_vcf, output_vcf)
    
    # Index the output
    logging.info(f"Indexing {output_vcf}...")
    subprocess.run(['bcftools', 'index', '-t', output_vcf], check=True)
    
    if os.path.exists(rejected_vcf):
        logging.info(f"Indexing {rejected_vcf}...")
        subprocess.run(['bcftools', 'index', '-t', rejected_vcf], check=True)
    
    logging.info("Liftover complete")
    return rejected_vcf


def parse_lifted_vcf(lifted_vcf, rejected_vcf):
    """
    Parse lifted VCF and extract liftover information
    
    Returns:
        dict: Mapping from original variant ID to lifted information
        list: List of rejected variant IDs
    """
    logging.info(f"Parsing lifted VCF...")
    
    lifted_variants = {}
    
    # Parse lifted variants - get SWAP and FLIP tags
    cmd = ['bcftools', 'query', '-f', '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/SWAP\t%INFO/FLIP\n', lifted_vcf]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) >= 5:
            variant_id = parts[0]
            chrom = parts[1]
            pos = parts[2]
            ref = parts[3]
            alt = parts[4]
            swap = parts[5] if len(parts) > 5 and parts[5] != '.' else '0'
            flip = parts[6] if len(parts) > 6 and parts[6] != '.' else '0'
            
            # Effect flipped if alleles were swapped by liftover
            effect_flipped = (swap == '1')
            
            lifted_variants[variant_id] = {
                'chr': chrom,
                'pos': int(pos),
                'ref': ref,
                'alt': alt,
                'swap': swap == '1',
                'flip': flip == '1',
                'effect_flipped': effect_flipped,  # True if effect direction needs flipping
                'status': 'lifted'
            }
    
    logging.info(f"Successfully lifted {len(lifted_variants)} variants")
    
    # Parse rejected variants
    rejected_variants = []
    if os.path.exists(rejected_vcf):
        cmd = ['bcftools', 'query', '-f', '%ID\n', rejected_vcf]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        rejected_variants = [line.strip() for line in result.stdout.strip().split('\n') if line.strip()]
        logging.info(f"Failed to lift {len(rejected_variants)} variants")
    
    return lifted_variants, rejected_variants


def update_sumstats(input_file, output_file, unmatched_file, lifted_info, rejected_ids,
                    chr_col, pos_col, ea_col, non_ea_col, effect_cols, eaf_cols=None, add_chr_prefix=False):
    """
    Update summary statistics with lifted coordinates and flip effects if needed
    
    Args:
        input_file: Original summary statistics file
        output_file: Output updated summary statistics file
        unmatched_file: Output file for unmatched variants
        lifted_info: Dictionary of lifted variant information
        rejected_ids: List of rejected variant IDs
        chr_col: Chromosome column name
        pos_col: Position column name
        ea_col: Effect allele column name
        non_ea_col: Non-effect allele column name
        effect_cols: List of effect columns to flip (e.g., Z, BETA)
        eaf_cols: List of effect allele frequency columns to flip (e.g., EAF, EAF_UKB)
        add_chr_prefix: Whether chr prefix was added during conversion
    """
    logging.info(f"Updating summary statistics...")
    
    # Force garbage collection before loading large dataset
    gc.collect()
    
    # Read original data with explicit dtypes to save memory
    logging.info(f"Reading {input_file}...")
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    
    # Clean chromosome names to match VCF format
    df['_chr_clean'] = df[chr_col].astype(str).str.strip()
    
    # Try to convert numeric-looking chromosomes (remove .0)
    try:
        numeric_df = pd.to_numeric(df['_chr_clean'], errors='coerce')
        numeric_mask = numeric_df.notna()
        df.loc[numeric_mask, '_chr_clean'] = numeric_df[numeric_mask].astype(int).astype(str)
    except:
        pass
    
    # Add chr prefix if it was added during conversion
    if add_chr_prefix:
        simple_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT', 'M']
        mask = df['_chr_clean'].isin(simple_chroms)
        df.loc[mask, '_chr_clean'] = 'chr' + df.loc[mask, '_chr_clean']
    
    # Create variant IDs for matching using chr:pos:ref:alt format (with uppercase alleles)
    df['_variant_id'] = (df['_chr_clean'] + ':' + df[pos_col].astype(str) + ':' + 
                        df[non_ea_col].astype(str).str.upper() + ':' + df[ea_col].astype(str).str.upper())
    
    # Add status columns
    df['LIFTOVER_STATUS'] = 'unknown'
    df['ALLELE_SWAP'] = False
    df['STRAND_FLIP'] = False
    df['LIFTED_CHR'] = df[chr_col].astype(str)
    # Convert position to numeric, allowing NaN values (will be replaced for lifted variants)
    df['LIFTED_POS'] = pd.to_numeric(df[pos_col], errors='coerce')
    # Add lifted allele columns (keep original alleles unchanged)
    ea_lifted_col = f'{ea_col}_lifted'
    non_ea_lifted_col = f'{non_ea_col}_lifted'
    df[ea_lifted_col] = ''
    df[non_ea_lifted_col] = ''
    df['AMBIGUOUS'] = False
    
    # Vectorized update using map operations
    # Create mask for lifted variants
    lifted_mask = df['_variant_id'].isin(lifted_info.keys())
    rejected_mask = df['_variant_id'].isin(rejected_ids)
    
    # Update lifted variants
    df.loc[lifted_mask, 'LIFTED_CHR'] = df.loc[lifted_mask, '_variant_id'].map(lambda x: lifted_info[x]['chr'])
    df.loc[lifted_mask, 'LIFTED_POS'] = df.loc[lifted_mask, '_variant_id'].map(lambda x: lifted_info[x]['pos']).astype('Int64')
    df.loc[lifted_mask, 'LIFTOVER_STATUS'] = 'lifted'
    df.loc[lifted_mask, 'ALLELE_SWAP'] = df.loc[lifted_mask, '_variant_id'].map(lambda x: lifted_info[x]['swap']).astype(bool)
    df.loc[lifted_mask, 'STRAND_FLIP'] = df.loc[lifted_mask, '_variant_id'].map(lambda x: lifted_info[x]['flip']).astype(bool)
    
    # Store lifted alleles in new columns (keep original alleles unchanged)
    df.loc[lifted_mask, non_ea_lifted_col] = df.loc[lifted_mask, '_variant_id'].map(lambda x: lifted_info[x]['ref'])
    df.loc[lifted_mask, ea_lifted_col] = df.loc[lifted_mask, '_variant_id'].map(lambda x: lifted_info[x]['alt'])
    
    # Get effect_flipped status for each variant
    df['_effect_flipped'] = False
    df.loc[lifted_mask, '_effect_flipped'] = df.loc[lifted_mask, '_variant_id'].map(lambda x: lifted_info[x]['effect_flipped']).astype(bool)
    
    # Determine AMBIGUOUS variants
    # AMBIGUOUS if: 1) alleles changed (A1:A2 != A1_lifted:A2_lifted), or 2) palindromic (A/T or C/G)
    def is_palindromic(a1, a2):
        a1_upper = str(a1).upper()
        a2_upper = str(a2).upper()
        return (a1_upper == 'A' and a2_upper == 'T') or (a1_upper == 'T' and a2_upper == 'A') or \
               (a1_upper == 'C' and a2_upper == 'G') or (a1_upper == 'G' and a2_upper == 'C')
    
    def alleles_changed(orig_a1, orig_a2, lift_a1, lift_a2):
        # Check if alleles changed (ignoring order)
        orig_set = {str(orig_a1).upper(), str(orig_a2).upper()}
        lift_set = {str(lift_a1).upper(), str(lift_a2).upper()}
        return orig_set != lift_set
    
    # Apply AMBIGUOUS flag for lifted variants
    for idx in df[lifted_mask].index:
        orig_a1 = df.at[idx, ea_col]
        orig_a2 = df.at[idx, non_ea_col]
        lift_a1 = df.at[idx, ea_lifted_col]
        lift_a2 = df.at[idx, non_ea_lifted_col]
        
        is_palindrome = is_palindromic(orig_a1, orig_a2)
        changed = alleles_changed(orig_a1, orig_a2, lift_a1, lift_a2)
        
        df.at[idx, 'AMBIGUOUS'] = is_palindrome or changed
    
    # Fix effect and frequency for effect_flipped variants (when alleles were swapped)
    flip_mask = lifted_mask & df['_effect_flipped']
    
    if flip_mask.any():
        # Swap A1 and A2 columns to keep A1 as effect allele
        temp_ea = df.loc[flip_mask, ea_col].copy()
        df.loc[flip_mask, ea_col] = df.loc[flip_mask, non_ea_col]
        df.loc[flip_mask, non_ea_col] = temp_ea
        
        # Flip effect sizes (change sign)
        if effect_cols:
            for col in effect_cols:
                if col in df.columns:
                    df.loc[flip_mask, col] = -df.loc[flip_mask, col]
        
        # Flip effect allele frequency (EAF -> 1-EAF)
        if eaf_cols:
            for col in eaf_cols:
                if col in df.columns:
                    df.loc[flip_mask, col] = 1 - df.loc[flip_mask, col]
    
    # Update rejected variants
    df.loc[rejected_mask, 'LIFTOVER_STATUS'] = 'rejected'
    
    # Calculate summary statistics
    matched_count = lifted_mask.sum()
    strand_flipped_count = (df['STRAND_FLIP'] == True).sum()
    flipped_count = flip_mask.sum()
    
    # Drop temporary columns
    df = df.drop(['_effect_flipped', '_chr_clean'], axis=1)
    
    # Update chromosome and position with lifted values
    df[chr_col] = df['LIFTED_CHR']
    df[pos_col] = df['LIFTED_POS']
    
    # Calculate summary statistics before splitting
    matched_count = lifted_mask.sum()
    strand_flipped_count = (df['STRAND_FLIP'] == True).sum()
    flipped_count = flip_mask.sum()
    ambiguous_count = (df['AMBIGUOUS'] == True).sum()
    total_count = len(df)
    
    # Write outputs in chunks to reduce memory usage
    logging.info(f"Writing lifted variants to {output_file}...")
    temp_output = output_file + '.tmp'
    # Determine compression based on output filename
    compression = 'gzip' if output_file.endswith('.gz') else None
    first_chunk = True
    for chunk in [df[df['LIFTOVER_STATUS'] == 'lifted']]:
        chunk_clean = chunk.drop(['_variant_id', 'LIFTED_CHR', 'LIFTED_POS'], axis=1)
        chunk_clean.to_csv(temp_output, sep='\t', index=False, mode='w' if first_chunk else 'a', header=first_chunk, compression=compression)
        first_chunk = False
        del chunk_clean
        gc.collect()
    # Only rename to final name after successful write
    if os.path.exists(output_file):
        os.remove(output_file)
    os.rename(temp_output, output_file)
    
    logging.info(f"Writing unmatched variants to {unmatched_file}...")
    temp_unmatched = unmatched_file + '.tmp'
    # Determine compression based on output filename
    compression_unmatched = 'gzip' if unmatched_file.endswith('.gz') else None
    first_chunk = True
    for chunk in [df[df['LIFTOVER_STATUS'] != 'lifted']]:
        chunk_clean = chunk.drop(['_variant_id', 'LIFTED_CHR', 'LIFTED_POS'], axis=1)
        chunk_clean.to_csv(temp_unmatched, sep='\t', index=False, mode='w' if first_chunk else 'a', header=first_chunk, compression=compression_unmatched)
        first_chunk = False
        del chunk_clean
        gc.collect()
    # Only rename to final name after successful write
    if os.path.exists(unmatched_file):
        os.remove(unmatched_file)
    os.rename(temp_unmatched, unmatched_file)
    
    logging.info(f"\nSummary:")
    logging.info(f"  Total variants: {total_count}")
    logging.info(f"  Successfully lifted: {matched_count}")
    logging.info(f"  Strand flipped: {strand_flipped_count}")
    logging.info(f"  Allele swapped: {flipped_count}")
    logging.info(f"  Ambiguous (palindromic or alleles changed): {ambiguous_count}")
    logging.info(f"  Failed to lift: {total_count - matched_count}")


def main():
    parser = argparse.ArgumentParser(
        description='Liftover summary statistics between genome builds (simple version without pre-normalization)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example:
  %(prog)s \\
    --input sumstats.txt.gz \\
    --output sumstats.hg38.txt.gz \\
    --chr-col CHR \\
    --pos-col POS \\
    --ea-col A1 \\
    --non-ea-col A2 \\
    --effect-col Z \\
    --eaf-col EAF_UKB \\
    --source-fasta hg19.fa.gz \\
    --target-fasta hg38.fa.gz \\
    --chain-file hg19ToHg38.over.chain.gz

Note: This version does NOT normalize REF alleles before liftover.
      Your input data should have correct REF alleles matching the source reference genome.
      Unmatched variants will be saved to <output>.unmatched.txt (or .txt.gz if output is compressed).
        '''
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input summary statistics file (txt or txt.gz)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output lifted summary statistics file')
    
    parser.add_argument('--chr-col', required=True,
                        help='Chromosome column name')
    parser.add_argument('--pos-col', required=True,
                        help='Position column name')
    parser.add_argument('--ea-col', required=True,
                        help='Effect allele column name')
    parser.add_argument('--non-ea-col', required=True,
                        help='Non-effect allele column name')
    parser.add_argument('--effect-col', action='append',
                        help='Effect column name(s) to flip (e.g., Z, BETA). Can specify multiple times.')
    parser.add_argument('--eaf-col', action='append',
                        help='Effect allele frequency column name(s) to flip when alleles swap (e.g., EAF, EAF_UKB). Can specify multiple times.')
    parser.add_argument('--source-fasta', required=True,
                        help='Source reference fasta file (for contig definitions)')
    parser.add_argument('--target-fasta', required=True,
                        help='Target reference fasta file')
    parser.add_argument('--chain-file', required=True,
                        help='Chain file for liftover')
    
    parser.add_argument('--add-chr-prefix', action='store_true',
                        help='Add "chr" prefix to chromosome names (use if source fasta has chr1, chr2, etc.)')
    
    parser.add_argument('--temp-dir', 
                        help='Temporary directory for intermediate files')
    parser.add_argument('--keep-temp', action='store_true',
                        help='Keep temporary files')
    
    args = parser.parse_args()
    
    # Set up logging
    log_file = args.output.replace('.txt.gz', '.log').replace('.txt', '.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stderr)
        ]
    )
    
    logging.info(f"Starting liftover process (simple version - no pre-normalization)")
    logging.info(f"Log file: {log_file}")
    
    # Create temp directory
    if args.temp_dir:
        os.makedirs(args.temp_dir, exist_ok=True)
        temp_dir = args.temp_dir
    else:
        temp_dir = tempfile.mkdtemp(prefix='liftover_sumstats_simple_')
    
    logging.info(f"Using temporary directory: {temp_dir}")
    
    try:
        # Step 1: Convert to VCF
        input_vcf = os.path.join(temp_dir, 'input.vcf.gz')
        if os.path.exists(input_vcf) and os.path.exists(input_vcf + '.tbi'):
            logging.info(f"Found existing VCF: {input_vcf}, skipping conversion...")
        else:
            sumstats_to_vcf(
                args.input, input_vcf,
                args.chr_col, args.pos_col, args.ea_col, args.non_ea_col,
                args.add_chr_prefix, args.source_fasta
            )
        
        # Step 2: Run liftover
        lifted_vcf = os.path.join(temp_dir, 'lifted.vcf.gz')
        rejected_vcf_path = os.path.join(temp_dir, 'lifted.rejected.vcf.gz')
        if os.path.exists(lifted_vcf) and os.path.exists(lifted_vcf + '.tbi'):
            logging.info(f"Found existing lifted VCF: {lifted_vcf}, skipping liftover...")
            rejected_vcf = rejected_vcf_path
        else:
            rejected_vcf = run_liftover(input_vcf, lifted_vcf, args.chain_file, args.source_fasta, args.target_fasta)
        
        # Step 3: Parse lifted VCF
        lifted_info, rejected_ids = parse_lifted_vcf(lifted_vcf, rejected_vcf)
        
        # Clear memory before final step
        gc.collect()
        
        # Step 4: Update summary statistics
        # Auto-generate unmatched filename from output filename
        if args.output.endswith('.txt.gz'):
            unmatched_file = args.output.replace('.txt.gz', '.unmatched.txt.gz')
        elif args.output.endswith('.gz'):
            unmatched_file = args.output.replace('.gz', '.unmatched.txt.gz')
        elif args.output.endswith('.txt'):
            unmatched_file = args.output.replace('.txt', '.unmatched.txt')
        else:
            unmatched_file = args.output + '.unmatched.txt'
        
        update_sumstats(
            args.input, args.output, unmatched_file,
            lifted_info, rejected_ids,
            args.chr_col, args.pos_col, args.ea_col, args.non_ea_col,
            args.effect_col or [],
            args.eaf_col or [],
            args.add_chr_prefix
        )
        
        logging.info("\nLiftover complete!")
        
    finally:
        # Clean up temp directory
        if not args.keep_temp and not args.temp_dir:
            import shutil
            shutil.rmtree(temp_dir)
            logging.info(f"Cleaned up temporary directory: {temp_dir}")


if __name__ == '__main__':
    main()
