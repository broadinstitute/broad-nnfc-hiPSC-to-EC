#! /usr/bin/env python

import pysam
import pandas as pd
import argparse
from collections import Counter

def load_gene_metadata(metadata_file, reverse_complement_primers=False):
    """Load gene metadata with intended targets and inner primers."""
    df = pd.read_csv(metadata_file)
    primer_to_gene = {}  # Maps primer sequence to intended target gene
    all_primers = set()  # All possible primers to search for
    
    for _, row in df.iterrows():
        gene_id = row['intended_target_gene_id']
        inner_primer = row['inner_primer']
        if pd.notna(gene_id) and pd.notna(inner_primer):
            # Remove version from gene ID (e.g., ENSG00000125691.14 -> ENSG00000125691)
            gene_id = gene_id.split('.')[0]
            
            # Reverse complement primer if requested
            if reverse_complement_primers:
                inner_primer = reverse_complement(inner_primer)
            primer_to_gene[inner_primer] = gene_id
            all_primers.add(inner_primer)
    
    return primer_to_gene, all_primers

def find_primer_in_sequence(sequence, primer, search_3_prime=False):
    """Find primer at the beginning or end of sequence with exact match."""
    if not sequence or not primer:
        return False
    
    primer_len = len(primer)
    
    # Check if sequence is long enough
    if len(sequence) < primer_len:
        return False
    
    # Extract the appropriate part of the sequence
    if search_3_prime:
        seq_region = sequence[-primer_len:]  # End of sequence
    else:
        seq_region = sequence[:primer_len]   # Beginning of sequence
    
    return seq_region == primer

# Translation table for reverse complement
DNA_COMPLEMENT_TABLE = str.maketrans("ATGC", "TACG")

def reverse_complement(seq):
    """Return reverse complement of DNA sequence using efficient string translation."""
    return seq.translate(DNA_COMPLEMENT_TABLE)[::-1]

def find_primer_locations(sequence, primer):
    """Find primer in sequence at beginning/end and forward/reverse orientations with exact matches."""
    if not sequence or not primer:
        return {}
    
    results = {
        'start_forward': False,
        'start_reverse': False,
        'end_forward': False,
        'end_reverse': False
    }
    
    primer_len = len(primer)
    if len(sequence) < primer_len:
        return results
    
    primer_rc = reverse_complement(primer)
    
    # Check beginning and end of sequence
    seq_start = sequence[:primer_len]
    seq_end = sequence[-primer_len:]
    
    results['start_forward'] = seq_start == primer
    results['start_reverse'] = seq_start == primer_rc
    results['end_forward'] = seq_end == primer
    results['end_reverse'] = seq_end == primer_rc
    
    return results

# bitmask of flags to EXCLUDE: unmapped, secondary, qc-fail, duplicate, supplementary
EXCLUDE_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800  # = 3844

def get_valid_read_with_genes(bam_path, max_reads=None, threads=0):
    """
    Yield (read, gene_id, seq) for the first `max_reads` salignments that:
      - are mapped
      - are primary (not secondary/supplementary)
      - are not duplicate / qc-fail
      - have a GX tag
      - have a sequence
    Streams linearly; no index required.
    """
    n_valid = 0
    with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam:
        for read in bam.fetch(until_eof=True):
            # fast flag filter
            if (read.flag & EXCLUDE_FLAGS) != 0:
                continue

            # require GX tag
            if not read.has_tag("GX"):
                continue

            # got a valid one
            yield read, read.get_tag("GX"), read.query_sequence
            n_valid += 1
            if max_reads and n_valid >= max_reads:
                break

def get_primer_length_range(all_primers):
    """Get the minimum and maximum primer lengths."""
    lengths = [len(primer) for primer in all_primers]
    return min(lengths), max(lengths)

def find_primer_in_sequence_optimized(sequence, primer_to_gene, min_len, max_len, search_3_prime=False):
    """Find primer using dictionary lookup with variable primer lengths."""
    if not sequence:
        return None
    
    # Try different primer lengths from min to max
    for primer_len in range(min_len, min(max_len + 1, len(sequence) + 1)):
        if search_3_prime:
            seq_region = sequence[-primer_len:]
        else:
            seq_region = sequence[:primer_len]
        
        if seq_region in primer_to_gene:
            return seq_region
    
    return None

def analyze_bam(bam_file, primer_to_gene, all_primers, output_file, search_3_prime=False):
    """Analyze BAM file for inner primer usage in read."""
    
    stats = {
        'total_read': 0,
        'read_with_gene': 0,
        'read_with_primer': 0,
        'read_correct_match': 0,
        'read_incorrect_match': 0
    }
    
    gene_assignment_counts = Counter()
    
    # Get primer length range for optimization
    min_len, max_len = get_primer_length_range(all_primers)
    
    # Buffer for writing results
    write_buffer = []
    buffer_size = 10000  # Write every 10k results
    
    # Open detailed results file for writing
    detailed_file = output_file.replace('.txt', '_detailed.tsv')
    with open(detailed_file, 'w', buffering=8192) as detail_f:
        detail_f.write("Inner_Primer_Found\tIntended_Target_Gene_ID\tObserved_Gene_ID\tMatch_Status\n")
        
        for read, observed_gene, read_sequence in get_valid_read_with_genes(bam_file, threads=4):
            stats['total_read'] += 1
            stats['read_with_gene'] += 1
            gene_assignment_counts[observed_gene] += 1
            
            # Account for read strand - if read is on reverse strand, reverse complement the sequence
            if read.is_reverse:
                search_sequence = reverse_complement(read_sequence)
            else:
                search_sequence = read_sequence
            
            # Search for primer using optimized method
            found_primer = find_primer_in_sequence_optimized(search_sequence, primer_to_gene, 
                                                           min_len, max_len, search_3_prime)
            
            if found_primer:
                stats['read_with_primer'] += 1
                intended_gene = primer_to_gene[found_primer]
                
                # Check if primer matches the gene assignment
                if intended_gene == observed_gene:
                    stats['read_correct_match'] += 1
                    match_status = 'correct'
                else:
                    stats['read_incorrect_match'] += 1
                    match_status = 'incorrect'
                
                # Buffer the write instead of writing immediately
                write_buffer.append(f"{found_primer}\t{intended_gene}\t{observed_gene}\t{match_status}\n")
                
                # Flush buffer when it reaches buffer_size
                if len(write_buffer) >= buffer_size:
                    detail_f.writelines(write_buffer)
                    write_buffer.clear()
        
        # Write any remaining buffered results
        if write_buffer:
            detail_f.writelines(write_buffer)
    
    # Write summary results with buffering
    with open(output_file, 'w', buffering=8192) as f:
        f.write("Inner Primer Usage Analysis Results\n")
        f.write("=" * 40 + "\n\n")
        
        f.write("Summary Statistics:\n")
        f.write(f"Total read processed: {stats['total_read']:,}\n")
        f.write(f"read with gene assignment: {stats['read_with_gene']:,}\n")
        f.write(f"read with primer found: {stats['read_with_primer']:,}\n")
        f.write(f"read with correct primer-gene match: {stats['read_correct_match']:,}\n")
        f.write(f"read with incorrect primer-gene match: {stats['read_incorrect_match']:,}\n")
        
        if stats['read_with_primer'] > 0:
            correct_rate = stats['read_correct_match'] / stats['read_with_primer'] * 100
            f.write(f"Primer-gene match accuracy: {correct_rate:.2f}%\n")
        
        f.write("\nTop 20 Gene Assignments:\n")
        for gene, count in gene_assignment_counts.most_common(20):
            f.write(f"{gene}: {count:,}\n")
        
        f.write(f"\nDetailed results written to: {detailed_file}\n")

def analyze_primer_positions(bam_file, primer_to_gene, all_primers, max_reads=1_000_000):
    """Analyze primer positions and orientations in first N sequences."""
    
    position_stats = {
        'start_forward': 0,
        'start_reverse': 0,
        'end_forward': 0,
        'end_reverse': 0,
        'total_with_gene': 0,
        'total_with_primer': 0
    }
    
    for read, _, read_sequence in get_valid_read_with_genes(bam_file, max_reads, threads=4):
        position_stats['total_with_gene'] += 1
        
        # Account for read strand - if read is on reverse strand, reverse complement the sequence
        if read.is_reverse:
            search_sequence = reverse_complement(read_sequence)
        else:
            search_sequence = read_sequence
        
        # Check all primers using the optimized locations function
        primer_found = False
        for primer in all_primers:
            locations = find_primer_locations(search_sequence, primer)
            if any(locations.values()):
                primer_found = True
                for pos, found in locations.items():
                    if found:
                        position_stats[pos] += 1
                break
        
        if primer_found:
            position_stats['total_with_primer'] += 1
    
    return position_stats, position_stats['total_with_gene']

def main():
    parser = argparse.ArgumentParser(description='Check inner primer usage in BAM file')
    parser.add_argument('bam_file', help='Input BAM file')
    parser.add_argument('metadata_file', help='Gene metadata CSV file')
    parser.add_argument('-o', '--output', default='primer_usage_report.txt', 
                       help='Output report file')
    parser.add_argument('--analyze-positions', action='store_true',
                       help='Analyze primer positions in first 1,000,000 sequences')
    parser.add_argument('--3-prime', action='store_true', dest='search_3_prime',
                       help='Search for primers at the 3\' end (default: 5\' end)')
    parser.add_argument('--reverse-complement', action='store_true', dest='reverse_complement_primers',
                       help='Reverse complement the primers from metadata before searching')

    args = parser.parse_args()

    print("Loading gene metadata...")
    primer_to_gene, all_primers = load_gene_metadata(args.metadata_file, args.reverse_complement_primers)
    print(f"Loaded {len(all_primers)} primers targeting {len(set(primer_to_gene.values()))} genes")
    print(f"Search parameters: position={'3-prime' if args.search_3_prime else '5-prime'}, "
          f"reverse_complement={args.reverse_complement_primers}")

    if args.analyze_positions:
        print("Analyzing primer positions in first 1,000,000 sequences...")
        pos_stats, total_reads = analyze_primer_positions(args.bam_file, primer_to_gene, all_primers)

        print(f"\nPrimer Position Analysis (first {total_reads} sequences):")
        print(f"Total read with gene assignment: {pos_stats['total_with_gene']}")
        print(f"Total read with matching primer found: {pos_stats['total_with_primer']}")
        print(f"Start forward: {pos_stats['start_forward']}")
        print(f"Start reverse complement: {pos_stats['start_reverse']}")
        print(f"End forward: {pos_stats['end_forward']}")
        print(f"End reverse complement: {pos_stats['end_reverse']}")

        with open('primer_position_analysis.txt', 'w') as f:
            f.write(f"Primer Position Analysis (first {total_reads} read sequences)\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total read with gene assignment: {pos_stats['total_with_gene']}\n")
            f.write(f"Total read with matching primer found: {pos_stats['total_with_primer']}\n")
            f.write(f"Start forward: {pos_stats['start_forward']}\n")
            f.write(f"Start reverse complement: {pos_stats['start_reverse']}\n")
            f.write(f"End forward: {pos_stats['end_forward']}\n")
            f.write(f"End reverse complement: {pos_stats['end_reverse']}\n")

        print("Position analysis saved to primer_position_analysis.txt")
    else:
        print("Analyzing BAM file...")
        analyze_bam(args.bam_file, primer_to_gene, all_primers, args.output, args.search_3_prime)
        print(f"Analysis complete. Results written to {args.output}")

if __name__ == "__main__":
    main()
