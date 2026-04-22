#!/usr/bin/env python3
"""
Extract the longest transcript per gene from a (Gencode) GTF and write a reduced
GTF that IGV can load, showing gene symbols (gene_name) rather than ENST IDs.

Definition of "longest": exon-union length (sum of non-overlapping exon bases).
Tie-break options can prefer APPRIS principal and/or Ensembl canonical tags.

Usage:
  python3 extract_longest_transcripts_gtf.py \
    --gtf ~/Annotations/HG38/GENCODE/gencode.v43.primary_assembly.annotation.gtf.gz \
    --genes genes.txt \
    --out genes.longest_tx.gtf

Optional:
  bgzip -c genes.longest_tx.gtf > genes.longest_tx.gtf.gz
  tabix -p gff genes.longest_tx.gtf.gz
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'\s*([^ ]+)\s+"([^"]+)"\s*;')


def parse_attrs(attr_str: str) -> dict:
    d = {}
    for m in ATTR_RE.finditer(attr_str):
        d[m.group(1)] = m.group(2)
    return d


def exon_union_length(exons):
    """exons: list of (start,end), 1-based inclusive. returns union length."""
    if not exons:
        return 0
    exons = sorted(exons)
    total = 0
    cs, ce = exons[0]
    for s, e in exons[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            total += ce - cs + 1
            cs, ce = s, e
    total += ce - cs + 1
    return total


def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def appris_rank(appris_str: str):
    """
    APPRIS strings look like:
      principal_1, principal_2, alternative_1, ...
    Lower is better.
    """
    if not appris_str:
        return (999, 999)
    m = re.search(r"principal_(\d+)", appris_str)
    if m:
        return (0, int(m.group(1)))
    m = re.search(r"alternative_(\d+)", appris_str)
    if m:
        return (1, int(m.group(1)))
    return (2, 999)


def main():
    ap = argparse.ArgumentParser(
        description="Extract longest transcript per gene (by exon-union length) from a GTF and output a reduced GTF for IGV."
    )
    ap.add_argument("--gtf", required=True, help="Input GTF (can be .gz)")
    ap.add_argument(
        "--genes",
        required=True,
        help="Gene list (one per line). Match by gene_name OR gene_id.",
    )
    ap.add_argument("--out", required=True, help="Output GTF")
    ap.add_argument(
        "--tie_break",
        choices=["appris", "canonical", "transcript_length", "lex"],
        default="appris",
        help=(
            "Tie-break if equal exon-union length. "
            "'appris' prefers lower appris_principal number; "
            "'canonical' prefers tag 'Ensembl_canonical'; "
            "'transcript_length' uses transcript span; "
            "'lex' uses transcript_id lexical."
        ),
    )
    args = ap.parse_args()

    # Load gene set
    genes = set()
    with open(args.genes) as f:
        for line in f:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.add(g)

    if not genes:
        sys.stderr.write("Gene list is empty.\n")
        sys.exit(2)

    # ---------- First pass: collect exons and transcript metadata for genes of interest ----------
    tx_exons = defaultdict(list)  # transcript_id -> list[(start,end)]
    tx_info = {}  # transcript_id -> dict with gene_id, gene_name, chr, strand, tags, appris
    gene_txs = defaultdict(set)  # gene_key (gene_id and gene_name) -> set(transcript_id)
    tx_span = defaultdict(lambda: [None, None])  # transcript_id -> [min_start, max_end]

    with open_maybe_gz(args.gtf) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue

            chrom, src, feature, start, end, score, strand, frame, attrs = cols
            a = parse_attrs(attrs)
            gene_id = a.get("gene_id")
            gene_name = a.get("gene_name")

            # match either gene_name or gene_id in user's list
            if (gene_name not in genes) and (gene_id not in genes):
                continue

            tid = a.get("transcript_id")
            if not tid:
                continue

            s = int(start)
            e = int(end)

            mn, mx = tx_span[tid]
            tx_span[tid] = [s if mn is None else min(mn, s), e if mx is None else max(mx, e)]

            if feature == "exon":
                tx_exons[tid].append((s, e))

            if tid not in tx_info:
                tags = a.get("tag", "")
                appris = a.get("appris", "")
                tx_info[tid] = {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "chrom": chrom,
                    "strand": strand,
                    "tag": tags,
                    "appris": appris,
                }

            # map gene -> transcripts (key by both id and name so user can supply either)
            if gene_id:
                gene_txs[gene_id].add(tid)
            if gene_name:
                gene_txs[gene_name].add(tid)

    # ---------- Choose longest transcript per gene ----------
    def make_tie_tuple(tid: str):
        info = tx_info.get(tid, {})
        tag = info.get("tag", "") or ""
        appr = info.get("appris", "") or ""
        span = tx_span[tid]
        span_len = 0
        if span[0] is not None and span[1] is not None:
            span_len = span[1] - span[0] + 1

        if args.tie_break == "appris":
            return (appris_rank(appr), 0 if "Ensembl_canonical" in tag else 1, -span_len, tid)
        if args.tie_break == "canonical":
            return (0 if "Ensembl_canonical" in tag else 1, appris_rank(appr), -span_len, tid)
        if args.tie_break == "transcript_length":
            return (-span_len, appris_rank(appr), 0 if "Ensembl_canonical" in tag else 1, tid)
        # lex
        return (tid,)

    selected = {}  # gene -> transcript_id
    for g in genes:
        cand = list(gene_txs.get(g, []))
        if not cand:
            continue

        best_tid = None
        best_sortkey = None
        for tid in cand:
            union_len = exon_union_length(tx_exons.get(tid, []))
            sortkey = (-union_len,) + make_tie_tuple(tid)  # maximize union_len
            if best_sortkey is None or sortkey < best_sortkey:
                best_sortkey = sortkey
                best_tid = tid

        if best_tid:
            selected[g] = best_tid

    selected_tids = set(selected.values())
    if not selected_tids:
        sys.stderr.write(
            "No transcripts selected. Check that your identifiers match gene_name or gene_id in the GTF.\n"
        )
        sys.exit(2)

    # Collect selected genes by gene_id for writing gene feature rows
    selected_gene_ids = set()
    for g, tid in selected.items():
        info = tx_info.get(tid, {})
        gid = info.get("gene_id")
        if gid:
            selected_gene_ids.add(gid)

    # ---------- Second pass: write gene rows + selected transcript rows + all features ----------
    # This ensures IGV shows gene symbols (gene_name) and has a proper gene container.
    with open_maybe_gz(args.gtf) as f, open(args.out, "w") as out:
        out.write(f"# Extracted longest transcript per gene from {args.gtf}\n")
        out.write(f"# Gene list: {args.genes}\n")
        out.write(f"# Longest = exon-union length; tie_break = {args.tie_break}\n")

        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue

            feature = cols[2]
            attrs = cols[8]
            a = parse_attrs(attrs)

            # Always keep gene feature line for selected genes (helps IGV label by gene_name)
            if feature == "gene":
                if a.get("gene_id") in selected_gene_ids:
                    out.write(line)
                continue

            # Keep selected transcript and all its associated features (exon/CDS/UTR/etc)
            tid = a.get("transcript_id")
            if tid and tid in selected_tids:
                out.write(line)

    # ---------- Report mapping ----------
    sys.stderr.write("Selected transcripts:\n")
    for g in sorted(genes):
        tid = selected.get(g)
        if tid:
            info = tx_info.get(tid, {})
            sys.stderr.write(
                f"  {g}\t{tid}\tgene_name={info.get('gene_name','')}\tgene_id={info.get('gene_id','')}\n"
            )
        else:
            sys.stderr.write(f"  {g}\tNOT_FOUND\n")


if __name__ == "__main__":
    main()
