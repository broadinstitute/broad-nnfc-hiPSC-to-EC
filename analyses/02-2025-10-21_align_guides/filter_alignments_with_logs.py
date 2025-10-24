#!/usr/bin/env python3
"""
sam_to_bed_and_mismatch_hist.py

Streams a SAM/BAM (name-sorted OK) and produces:
  • BED (only alignments that PASS filters)
  • TSV of DISCARDED mapped alignments (fail filters)
  • TSV of UNMAPPED alignments
  • CSV of mismatch positions (kept + discarded), normalized by PAM end
  • Optional Altair HTML chart of mismatches (positions as factors), faceted by reason

Filters (enable with flags):
  --require-pam-clean   Keep only alignments whose last 2 query bases are aligned and NOT mismatches.

Optional BED trimming:
  --emit-protospacer    BED spans protospacer only (drop last 3 read-end bases).
  --drop-leading-g      With --emit-protospacer, also drop 1 leading base (e.g., U6 'G').

Outputs (paths configurable via flags):
  --out-bed BED                (kept only)
  --out-discarded TSV          (mapped but failed filter)
  --out-unmapped TSV           (unmapped)
  --out-mismatch-csv CSV       (kept + discarded, for plotting)
  --plot-html HTML             (Altair factor chart; requires pandas+altair)

Notes:
  - Mismatch placement uses MD + CIGAR; no FASTA required.
  - Query length comes from CIGAR (works even when SEQ='*').
  - BED coordinates are 0-based, half-open.
"""

import sys
import re
import csv
import argparse
import logging
from typing import List, Tuple, Optional, Dict, Set, Union
from dataclasses import dataclass
from pathlib import Path
import pysam  # pip install pysam

# Optional (only if --plot-html is used)
try:
    import pandas as pd
    import altair as alt
except ImportError:
    pd = None
    alt = None

# Constants
PAM_LENGTH = 3
MIN_PAM_CLEAN_BASES = 2
MISSING_SEQUENCE_PLACEHOLDER = '*'

# DNA complement mapping
COMPLEMENT_MAP = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')

# CIGAR operations that consume query sequence
QUERY_CONSUMING_OPS = frozenset(['M', 'I', 'S', '=', 'X'])
ALIGNMENT_OPS = frozenset(['M', '=', 'X'])
INSERTION_SOFT_CLIP_OPS = frozenset(['I', 'S'])
DELETION_SKIP_OPS = frozenset(['D', 'N'])

# Output file headers
DISCARDED_HEADER = "read_name\tchromosome\tpos1\tflag\tmapq\tstrand\tcigar\tNM\tAS\tMD\treason\n"
UNMAPPED_HEADER = "read_name\tflag\n"
MISMATCH_FIELDNAMES = [
    "read_name", "chromosome", "strand", "read_length",
    "query_index", "distance_from_pam", "nm", "as", "reason"
]

# Regex patterns
CIGAR_PATTERN = re.compile(r'(\d+)([MIDNSHP=X])')
MD_PATTERN = re.compile(r'(\d+|\^[A-Za-z]+|[A-Za-z])')

@dataclass
class AlignmentStats:
    """Track alignment processing statistics."""
    total_seen: int = 0
    total_kept: int = 0
    total_discarded: int = 0
    total_unmapped: int = 0

@dataclass
class MismatchRecord:
    """Structured mismatch data for CSV output."""
    read_name: str
    chromosome: str
    strand: str
    read_length: int
    query_index: int
    distance_from_pam: int
    nm: int
    as_score: int
    reason: str

def setup_logging() -> None:
    """Configure logging for the application."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s: %(message)s',
        stream=sys.stderr
    )

def parse_cigar(cigar_string: str) -> List[Tuple[int, str]]:
    """Parse CIGAR string into list of (length, operation) tuples."""
    if not cigar_string:
        return []
    return [(int(length), op) for length, op in CIGAR_PATTERN.findall(cigar_string)]

def parse_md_tag(md_tag: str) -> List[str]:
    """Parse MD tag into list of tokens (matches, mismatches, deletions)."""
    return MD_PATTERN.findall(md_tag) if md_tag else []

def calculate_query_length_from_cigar(cigar_ops: List[Tuple[int, str]]) -> int:
    """Calculate query sequence length from CIGAR operations."""
    return sum(length for length, op in cigar_ops if op in QUERY_CONSUMING_OPS)

def build_query_to_ref_map(ref_start_1based: int, cigar_string: str) -> Tuple[List[Optional[int]], int]:
    """
    Build mapping from query indices to reference coordinates.
    
    Args:
        ref_start_1based: Reference start position (1-based from SAM)
        cigar_string: CIGAR string from alignment
        
    Returns:
        Tuple of (query_to_ref_map, query_length)
        - query_to_ref_map: List where index=query_pos, value=ref_pos (0-based) or None
        - query_length: Total query sequence length
    """
    cigar_ops = parse_cigar(cigar_string)
    query_len = calculate_query_length_from_cigar(cigar_ops)
    
    if query_len == 0:
        return [], 0
        
    query_to_ref_map = [None] * query_len
    query_pos = 0
    ref_pos = ref_start_1based - 1  # Convert to 0-based

    for length, op in cigar_ops:
        if op in ALIGNMENT_OPS:
            # Match/mismatch: both query and reference advance
            for _ in range(length):
                query_to_ref_map[query_pos] = ref_pos
                query_pos += 1
                ref_pos += 1
        elif op in INSERTION_SOFT_CLIP_OPS:
            # Insertion/soft-clip: only query advances
            query_pos += length
        elif op in DELETION_SKIP_OPS:
            # Deletion/skip: only reference advances
            ref_pos += length
        # H/P operations consume neither query nor reference

    return query_to_ref_map, query_len

def find_mismatch_query_indices(md_tag: str, cigar_string: str) -> Set[int]:
    """
    Identify query positions that are mismatches using MD tag and CIGAR.
    
    Args:
        md_tag: MD tag from SAM record
        cigar_string: CIGAR string from SAM record
        
    Returns:
        Set of query indices (0-based) that are mismatches
    """
    mismatches = set()
    if not md_tag:
        return mismatches
        
    cigar_ops = parse_cigar(cigar_string)
    md_tokens = parse_md_tag(md_tag)
    
    token_index = 0
    query_position = 0
    
    for length, op in cigar_ops:
        if op in ALIGNMENT_OPS:
            remaining = length
            while remaining > 0 and token_index < len(md_tokens):
                token = md_tokens[token_index]
                
                if token.isdigit():
                    # Match run
                    match_count = int(token)
                    consumed = min(match_count, remaining)
                    query_position += consumed
                    remaining -= consumed
                    
                    if match_count > consumed:
                        md_tokens[token_index] = str(match_count - consumed)
                    else:
                        token_index += 1
                        
                elif token.startswith('^'):
                    # Deletion in reference (skip)
                    token_index += 1
                else:
                    # Mismatch
                    mismatches.add(query_position)
                    query_position += 1
                    remaining -= 1
                    token_index += 1
                    
            # Consume any remaining positions as matches
            query_position += remaining
            
        elif op in INSERTION_SOFT_CLIP_OPS:
            query_position += length

    return mismatches

def get_aligned_span_on_reference(query_to_ref_map: List[Optional[int]]) -> Optional[Tuple[int, int]]:
    """
    Get the reference span covered by aligned query positions.
    
    Returns:
        Tuple of (start, end) in 0-based half-open coordinates, or None if no alignment
    """
    aligned_positions = [pos for pos in query_to_ref_map if pos is not None]
    if not aligned_positions:
        return None
    return min(aligned_positions), max(aligned_positions) + 1

def map_query_subspan_to_reference(
    query_to_ref_map: List[Optional[int]], 
    query_start: int, 
    query_end: int
) -> Optional[Tuple[int, int]]:
    """
    Map a query subsequence to reference coordinates.
    
    Args:
        query_to_ref_map: Mapping from query to reference positions
        query_start: Start of query region (inclusive)
        query_end: End of query region (exclusive)
        
    Returns:
        Tuple of (ref_start, ref_end) in 0-based half-open coords, or None if invalid
    """
    if not query_to_ref_map or query_start >= query_end:
        return None
        
    # Find first aligned position in range
    ref_start = None
    for i in range(query_start, min(query_end, len(query_to_ref_map))):
        if query_to_ref_map[i] is not None:
            ref_start = query_to_ref_map[i]
            break
    
    # Find last aligned position in range
    ref_end = None
    for i in range(min(query_end - 1, len(query_to_ref_map) - 1), query_start - 1, -1):
        if query_to_ref_map[i] is not None:
            ref_end = query_to_ref_map[i] + 1
            break
    
    if ref_start is not None and ref_end is not None and ref_end > ref_start:
        return ref_start, ref_end
    return None

def is_pam_region_clean(
    query_to_ref_map: List[Optional[int]], 
    query_length: int, 
    mismatch_indices: Set[int],
    is_reverse_strand: bool = False
) -> bool:
    """
    Check if PAM region is aligned without mismatches.
    
    Args:
        query_to_ref_map: Query to reference position mapping
        query_length: Total query sequence length
        mismatch_indices: Set of query positions that are mismatches
        is_reverse_strand: True if alignment is on minus strand
        
    Returns:
        True if PAM region is clean (aligned and no mismatches)
        
    Note:
        - Plus strand: PAM is at the END (last MIN_PAM_CLEAN_BASES)
        - Minus strand: PAM is at the BEGINNING (first MIN_PAM_CLEAN_BASES)
    """
    if query_length < MIN_PAM_CLEAN_BASES:
        return False
    
    if is_reverse_strand:
        # Minus strand: PAM at beginning of query sequence
        pam_indices = range(0, MIN_PAM_CLEAN_BASES)
    else:
        # Plus strand: PAM at end of query sequence  
        pam_indices = range(query_length - MIN_PAM_CLEAN_BASES, query_length)
    
    # Check alignment
    if not all(i < len(query_to_ref_map) and query_to_ref_map[i] is not None 
               for i in pam_indices):
        return False
    
    # Check for mismatches
    return not any(i in mismatch_indices for i in pam_indices)

def determine_discard_reason(
    query_to_ref_map: List[Optional[int]], 
    query_length: int, 
    mismatch_indices: Set[int],
    is_reverse_strand: bool = False
) -> str:
    """
    Determine specific reason for discarding alignment based on PAM region.
    
    Args:
        query_to_ref_map: Query to reference position mapping
        query_length: Total query sequence length
        mismatch_indices: Set of query positions that are mismatches
        is_reverse_strand: True if alignment is on minus strand
    """
    if query_length < MIN_PAM_CLEAN_BASES:
        return "discarded_tail_unaligned"
    
    if is_reverse_strand:
        # Minus strand: PAM at beginning
        pam_indices = range(0, MIN_PAM_CLEAN_BASES)
    else:
        # Plus strand: PAM at end
        pam_indices = range(query_length - MIN_PAM_CLEAN_BASES, query_length)
    
    pam_aligned = all(i < len(query_to_ref_map) and query_to_ref_map[i] is not None 
                     for i in pam_indices)
    
    if not pam_aligned:
        return "discarded_tail_unaligned"
    else:
        return "discarded_tail_mismatch"


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    
    Args:
        sequence: DNA sequence string (supports ATCGN, case insensitive)
        
    Returns:
        Reverse complement of the input sequence
        
    Note:
        - Handles both uppercase and lowercase bases
        - 'N' bases are preserved as 'N'
        - Empty or invalid sequences return empty string
    """
    if not sequence or sequence == MISSING_SEQUENCE_PLACEHOLDER:
        return ''
    
    return sequence.translate(COMPLEMENT_MAP)[::-1]

class SequenceCache:
    """Cache for sequence data with size limit to prevent memory issues."""
    
    def __init__(self, max_size: int = 10000):
        self.cache = {}
        self.max_size = max_size
    
    def get(self, read_name: str, default: str = '') -> str:
        """Get sequence for read name."""
        return self.cache.get(read_name, default)
    
    def set(self, read_name: str, sequence: str) -> None:
        """Store sequence with cache size management."""
        if len(self.cache) >= self.max_size:
            # Remove oldest entries (simple FIFO)
            to_remove = list(self.cache.keys())[:len(self.cache) // 2]
            for key in to_remove:
                del self.cache[key]
        
        self.cache[read_name] = sequence

def create_output_files(args) -> Tuple:
    """Create and return output file handles with proper headers."""
    bed_out = open(args.out_bed, "w", buffering=1)
    disc_out = open(args.out_discarded, "w", buffering=1)
    unmap_out = open(args.out_unmapped, "w", buffering=1)
    
    # Write headers
    disc_out.write(DISCARDED_HEADER)
    unmap_out.write(UNMAPPED_HEADER)
    
    return bed_out, disc_out, unmap_out

def has_indels_in_region(cigar_string: str, query_start: int, query_end: int) -> bool:
    """
    Check if there are insertions or deletions in a specific query region.
    
    Args:
        cigar_string: CIGAR string from alignment
        query_start: Start of region to check (inclusive)
        query_end: End of region to check (exclusive)
        
    Returns:
        True if there are insertions or deletions in the specified region
    """
    cigar_ops = parse_cigar(cigar_string)
    query_pos = 0
    
    for length, op in cigar_ops:
        if op in QUERY_CONSUMING_OPS:
            # This operation consumes query sequence
            region_start = max(query_pos, query_start)
            region_end = min(query_pos + length, query_end)
            
            # Check if this operation overlaps with our region of interest
            if region_start < region_end and op in ['I', 'D']:
                return True
                
            query_pos += length
        elif op in DELETION_SKIP_OPS:
            # Deletion - check if we're in the region of interest
            # For deletions, we need to check if the current query position is in our region
            if query_start <= query_pos < query_end:
                return True
    
    return False

def process_alignment_record(
    aln, 
    args, 
    seq_cache: SequenceCache,
    mismatch_rows: List[Dict],
    stats: AlignmentStats
) -> Optional[Tuple]:
    """
    Process a single alignment record and return BED data if kept.
    
    Returns:
        Tuple of (bed_line_data) if kept, None if discarded
    """
    read_name = aln.query_name
    chrom = aln.reference_name
    ref_start_1b = aln.reference_start + 1
    mapq = aln.mapping_quality
    strand = '-' if aln.is_reverse else '+'
    is_reverse_strand = aln.is_reverse
    cigar = aln.cigarstring or '*'
    tags = dict(aln.tags)
    nm = int(tags.get('NM', -1))
    as_score = int(tags.get('AS', 0))
    md = tags.get('MD', None)
    

    # Build mappings and find mismatches
    query_to_ref_map, query_length = build_query_to_ref_map(ref_start_1b, cigar)
    mismatch_indices = find_mismatch_query_indices(md, cigar)

    # Determine if alignment passes filters
    reason = "kept"
    if args.require_pam_clean and not is_pam_region_clean(query_to_ref_map, query_length, mismatch_indices, is_reverse_strand):
        reason = determine_discard_reason(query_to_ref_map, query_length, mismatch_indices, is_reverse_strand)

    # Record mismatches for analysis
    for query_index in mismatch_indices:
        if 0 <= query_index < query_length and query_to_ref_map[query_index] is not None:
            # Calculate distance from PAM based on strand
            if is_reverse_strand:
                # Minus strand: distance from beginning (PAM at start)
                distance_from_pam = query_index
            else:
                # Plus strand: distance from end (PAM at end)
                distance_from_pam = (query_length - 1) - query_index
                
            mismatch_rows.append({
                "read_name": read_name,
                "chromosome": chrom,
                "strand": strand,
                "read_length": query_length,
                "query_index": query_index,
                "distance_from_pam": distance_from_pam,
                "nm": nm,
                "as": as_score,
                "reason": reason,
            })

    # Update statistics
    if reason != "kept":
        stats.total_discarded += 1
        return None, (read_name, chrom, ref_start_1b, aln.flag, mapq, strand, cigar, nm, as_score, md or '', reason)

    # Process kept alignment for BED output
    span = get_aligned_span_on_reference(query_to_ref_map)
    if span is None:
        stats.total_discarded += 1
        return None, (read_name, chrom, ref_start_1b, aln.flag, mapq, strand, cigar, nm, as_score, md or '', "discarded_tail_unaligned")

    bed_start, bed_end = span

    
    # Determine query subsequence to extract based on whether we're emitting protospacer
    if args.emit_protospacer:
        if is_reverse_strand:
            # Minus strand: PAM at beginning
            # Skip PAM_LENGTH bases at beginning, optionally skip leading G.
            skip_leading_g = 1 if args.drop_leading_g else 0
            query_start = PAM_LENGTH
            query_end = query_length - skip_leading_g
        else:
            # Plus strand: PAM at end, protospacer at beginning  
            # Optionally skip leading G, then exclude PAM_LENGTH bases at end
            query_start = 1 if args.drop_leading_g else 0
            query_end = max(0, query_length - PAM_LENGTH)
        
        # Check for indels in the protospacer region
        if has_indels_in_region(cigar, query_start, query_end):
            stats.total_discarded += 1
            return None, (read_name, chrom, ref_start_1b, aln.flag, mapq, strand, cigar, nm, as_score, md or '', "discarded_protospacer_indel")
        
        # Map protospacer region to reference coordinates
        protospacer_span = map_query_subspan_to_reference(query_to_ref_map, query_start, query_end)
        if protospacer_span is not None:
            bed_start, bed_end = protospacer_span
        else:
            stats.total_discarded += 1
            return None, (read_name, chrom, ref_start_1b, aln.flag, mapq, strand, cigar, nm, as_score, md or '', "discarded_tail_unaligned")
    else:
        # If not emitting protospacer, use full sequence and full span
        query_start, query_end = 0, query_length

    stats.total_kept += 1
    
    return (chrom, bed_start, bed_end, read_name, mapq, strand, nm, as_score), None

def write_mismatch_csv(filename: str, mismatch_rows: List[Dict]) -> None:
    """Write mismatch data to CSV file."""
    with open(filename, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=MISMATCH_FIELDNAMES)
        writer.writeheader()
        writer.writerows(mismatch_rows)

def create_mismatch_plot(mismatch_rows: List[Dict], output_path: str) -> None:
    """Create and save Altair mismatch position plot."""
    if pd is None or alt is None:
        logging.warning("Plotting requested but pandas/altair not installed. Skipping chart.")
        return
        
    df = pd.DataFrame(mismatch_rows)
    if df.empty:
        logging.warning("No mismatches to plot. Skipping chart.")
        return
        
    # Prepare data for plotting
    df["distance_from_pam"] = df["distance_from_pam"].astype(int)
    df = df.sort_values("distance_from_pam", ascending=False)
    ordered_labels = [str(x) for x in df["distance_from_pam"].unique()]
    df["distance_from_pam"] = df["distance_from_pam"].astype(str)

    # Create chart
    base_chart = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X("distance_from_pam:N", sort=ordered_labels,
                    title="Distance from PAM end (0 = last base)"),
            y=alt.Y("count():Q", title="Mismatch count"),
            tooltip=["reason:N", "distance_from_pam:N", "count():Q"],
            color=alt.Color("reason:N", legend=None),
        )
        .properties(width=700, height=150)
    )

    final_chart = (
        base_chart
        .facet(row=alt.Row("reason:N",
                          sort=["kept", "discarded_tail_mismatch", "discarded_tail_unaligned"],
                          title=None))
        .resolve_scale(y='independent')
        .properties(title="Mismatch positions by exact distance from PAM (kept vs discard reasons)")
    )
    
    final_chart.save(output_path)
    logging.info(f"Wrote Altair chart to {output_path}")

def main():
    setup_logging()
    
    ap = argparse.ArgumentParser(
        description="SAM/BAM → BED (kept) + TSV (discarded/unmapped) + mismatch CSV; optional PAM check and Altair factor chart."
    )
    ap.add_argument("input_sam", help="Input SAM or BAM (name-sorted OK).")
    ap.add_argument("--out-bed", default="kept.bed", help="Output BED (ONLY kept alignments).")
    ap.add_argument("--out-discarded", default="discarded.tsv", help="Output TSV for discarded (mapped but failed).")
    ap.add_argument("--out-unmapped", default="unmapped.tsv", help="Output TSV for unmapped alignments.")
    ap.add_argument("--out-mismatch-csv", default="mismatches.csv", help="Output mismatch CSV (kept+discarded).")
    ap.add_argument("--require-pam-clean", action="store_true",
                    help="Keep only alignments where last two query bases are aligned and not mismatches.")
    ap.add_argument("--emit-protospacer", action="store_true",
                    help="BED spans protospacer only (drop last 3 read-end bases).")
    ap.add_argument("--drop-leading-g", action="store_true",
                    help="When --emit-protospacer, also drop 1 leading base (e.g., U6 'G').")
    ap.add_argument("--plot-html", default=None,
                    help="Also render Altair chart (positions as factors) to this HTML file.")
    ap.add_argument("--no-dedup", action="store_true",
                    help="Disable deduplication of identical alignment records.")
    args = ap.parse_args()
    
    logging.info(f"Running with: emit_protospacer={args.emit_protospacer}, drop_leading_g={args.drop_leading_g}, require_pam_clean={args.require_pam_clean}")
    
    # Validate input file exists
    if not Path(args.input_sam).exists():
        logging.error(f"Input file does not exist: {args.input_sam}")
        sys.exit(1)

    # Open input file
    mode = "rb" if args.input_sam.endswith(".bam") else "r"
    try:
        samfile = pysam.AlignmentFile(args.input_sam, mode)
    except Exception as e:
        logging.error(f"Failed to open SAM/BAM file: {e}")
        sys.exit(1)

    # Initialize data structures
    bed_out, disc_out, unmap_out = create_output_files(args)
    seq_cache = SequenceCache()
    mismatch_rows = []
    seen_alignments = set() if not args.no_dedup else None
    stats = AlignmentStats()

    try:
        # Process alignments
        for aln in samfile.fetch(until_eof=True):
            if aln.is_unmapped:
                stats.total_unmapped += 1
                unmap_out.write(f"{aln.query_name}\t{aln.flag}\n")
                continue

            stats.total_seen += 1

            # Optional deduplication
            if seen_alignments is not None:
                alignment_key = (aln.query_name, aln.reference_name, aln.reference_start, 
                               aln.is_reverse, aln.cigarstring)
                if alignment_key in seen_alignments:
                    continue
                seen_alignments.add(alignment_key)

            # Process alignment
            #if aln.query_name != "WTC11_Random_Screen_Crop_14729":
            #    continue
            bed_data, discard_data = process_alignment_record(aln, args, seq_cache, mismatch_rows, stats)
            
            if bed_data is not None:
                bed_out.write("\t".join(map(str, bed_data)) + "\n")
            elif discard_data is not None:
                disc_out.write("\t".join(map(str, discard_data)) + "\n")

    finally:
        # Clean up file handles
        bed_out.close()
        disc_out.close()
        unmap_out.close()
        samfile.close()

    # Write outputs
    write_mismatch_csv(args.out_mismatch_csv, mismatch_rows)
    
    # Report statistics
    logging.info(f"Mapped processed: {stats.total_seen} | Kept(BED): {stats.total_kept} | "
                f"Discarded: {stats.total_discarded} | Unmapped: {stats.total_unmapped}")
    logging.info(f"Mismatches written: {len(mismatch_rows)} to {args.out_mismatch_csv}")

    # Generate plot if requested
    if args.plot_html:
        create_mismatch_plot(mismatch_rows, args.plot_html)

if __name__ == "__main__":
    main()
