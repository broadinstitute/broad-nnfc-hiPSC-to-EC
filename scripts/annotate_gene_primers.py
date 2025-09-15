import csv
import pysam
from tqdm import tqdm
from multiprocessing import Pool


# Top-level helper functions for multiprocessing

def process_transcript_batch(args):
    """
    Processes a batch of transcript sequences and finds matches for outer and/or inner primers.

    Args:
        args: Tuple containing:
            transcript_batch: List of (transcript_id, transcript_seq) pairs.
            primer_rows: List of primer rows from the CSV.
            outer_col: Name of the outer primer column.
            inner_col: Name of the inner primer column.
            target_col: Name of the target column.
            transcript_id_fields: List of transcript ID field names to output.
            transcript_id_indices: List of indices for each transcript ID field in the pipe-separated transcript_id.

    Returns:
        List of dictionaries, each representing a match with primer info, distances, and transcript ID fields.
    """
    transcript_batch, primer_rows, outer_col, inner_col, target_col, transcript_id_fields, transcript_id_indices = args
    batch_results = []
    for transcript_id, transcript_seq in transcript_batch:
        transcript_id_parts = transcript_id.split("|")
        for row in primer_rows:
            outer_primer = row.get(outer_col, None)
            inner_primer = row.get(inner_col, None)
            intended_target = row.get(target_col, None)
            outer_pos = transcript_seq.find(outer_primer)
            inner_pos = transcript_seq.find(inner_primer)
            primers_found = "inner"
            if inner_pos == -1:
                continue
            if outer_pos != -1:
                primers_found = "both"
            if primers_found:
                outer_dist = len(transcript_seq) - (outer_pos + len(outer_primer)) if outer_pos != -1 else -1
                inner_dist = len(transcript_seq) - (inner_pos + len(inner_primer))
                rel_dist = abs(outer_pos - inner_pos) if outer_pos != -1 and inner_pos != -1 else -1
                is_intended = intended_target in transcript_id_parts
                result = {
                    "outer_primer": outer_primer,
                    "inner_primer": inner_primer,
                    "outer_distance_from_transcript_end": outer_dist,
                    "inner_distance_from_transcript_end": inner_dist,
                    "distance_between_primers": rel_dist,
                    "primers_found": primers_found,
                    "is_intended_target": is_intended,
                    "intended_target_gene_symbol": intended_target
                }
                for field, idx in zip(transcript_id_fields, transcript_id_indices):
                    result[field] = transcript_id_parts[idx] if idx < len(transcript_id_parts) else None
                batch_results.append(result)
    return batch_results


def find_primers_with_pysam(primer_file, fasta_file, output_file, outer_col=None, inner_col=None, target_col=None, max_workers=8, batch_size=1000):
    """
    Finds outer and/or inner primers in transcript sequences and calculates their distances from the end.
    Reports matches where only outer, only inner, or both primers are found in a transcript.
    Outputs a CSV file with primer info, distances, and selected transcript ID fields.

    Args:
        primer_file (str): Path to the CSV file with primer sequences.
        fasta_file (str): Path to the FASTA file with transcript sequences.
        output_file (str): Path to the output CSV file for matches.
        outer_col (str, optional): The name of the column containing outer primer sequences.
        inner_col (str, optional): The name of the column containing inner primer sequences.
        target_col (str, optional): The name of the column containing intended target gene symbol.
        max_workers (int, optional): Number of parallel processes to use. Default is 8.
        batch_size (int, optional): Number of transcripts to process per batch. Default is 1000.

    Returns:
        None. Writes results to output_file.
    """
    with open(primer_file, 'r', newline='') as f:
        reader = csv.DictReader(f)
        primer_rows = list(reader)

    # Only output these transcript_id fields
    transcript_id_fields = ["transcript_id", "gene_id", "gene_symbol", "gene_type"]
    # Indices for these fields in the pipe-separated transcript_id
    with pysam.FastaFile(fasta_file) as fasta_reader:
        all_transcript_ids = fasta_reader.references
        all_possible_fields = [
            "transcript_id", "gene_id", "otthumg_id", "otthumt_id", "transcript_name", "gene_symbol", "gene_length", "gene_type"
        ]
        transcript_id_indices = [all_possible_fields.index(f) for f in transcript_id_fields]
        transcript_seqs = [(tid, fasta_reader.fetch(reference=tid)) for tid in all_transcript_ids]

    fieldnames = [
        "outer_primer", "inner_primer",
        "outer_distance_from_transcript_end", "inner_distance_from_transcript_end",
        "distance_between_primers", "primers_found", "is_intended_target", "intended_target_gene_symbol"
    ] + transcript_id_fields

    batches = [transcript_seqs[i:i+batch_size] for i in range(0, len(transcript_seqs), batch_size)]
    pool_args = [(batch, primer_rows, outer_col, inner_col, target_col, transcript_id_fields, transcript_id_indices) for batch in batches]
    all_results = []
    with Pool(processes=max_workers) as pool:
        for matches in tqdm(pool.imap(process_transcript_batch, pool_args), total=len(batches), desc="Processing transcript batches"):
            if matches:
                all_results.extend(matches)

    with open(output_file, 'w', newline='') as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames)
        writer.writeheader()
        for match in all_results:
            writer.writerow(match)


# Example usage:
if __name__ == "__main__":
    primer_csv = "../metadata/hiPSC-EC/metadata_genes.csv"
    fasta_file = "/Users/emattei/Annotations/HG38/GENCODE/gencode.v44.transcripts.fa.bgz"
    output_file = "../metadata/hiPSC-EC/primer_matches.csv"

    find_primers_with_pysam(primer_csv, fasta_file, output_file, outer_col="outer_primer", inner_col="inner_primer", target_col="intended_target_gene_symbol")

    print(f"Results written to {output_file}")
