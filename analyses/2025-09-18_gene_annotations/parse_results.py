# Inspect third column of ../../results/2025-09-18_gene_annotations/gene_annotations.csv and extract keyword-like phrases.
import pandas as pd
import re
import os

path = "../../results/2025-09-18_gene_annotations/gene_annotations.csv"
outdir = "../../results/2025-09-18_gene_annotations/"
if not os.path.isfile(path):
    raise FileNotFoundError(f"File not found: {path}")

df: pd.DataFrame = pd.read_csv(path)

# Get third column by position (0-based index 2)
if df.shape[1] < 3:
    raise SystemExit(f"File has only {df.shape[1]} columns; need at least 3.")
colname: str = df.columns[2]
col: pd.Series = df.iloc[:, 2].astype(str).fillna("")
sample_vals: list = col.head(10).tolist()

# Normalize
text: pd.Series = col.str.lower()

# Simple pattern-based categories to look for (expandable)
patterns: dict[str, str] = {
    "mitochondrial": r"\bmitochond\w*",
    "ribosomal": r"\bribosom\w*",
    "pseudogene": r"\bpseudogene\b|\bpseudo[- ]?gene\b",
    "lncrna": r"\blnc[- ]?rna\b|\blong noncoding\b",
    "mirna": r"\bmi[- ]?rna\b|\bmicro[- ]?rna\b",
    "snrna": r"\bsn[- ]?rna\b",
    "snorna": r"\bsno[- ]?rna\b",
    "rrna": r"\br[- ]?rna\b(?!\w)",
    "antisense": r"\bantisense\b",
    "protein_coding": r"\bprotein[- ]?coding\b",
    "transcription_factor": r"\btranscription (factor|regulator)\b|\bhomeobox\b|\bhox\b|\bznf\b|\bzinc finger\b",
    "kinase": r"\bkinase\b|\bpk\b",
    "phosphatase": r"\bphosphatase\b",
    "gpcr": r"\bgpcr\b|\bg protein[- ]?coupled\b",
    "ion_channel": r"\bion channel\b|\bchannel\b",
    "transport": r"\btransporter?\b|\bsolute carrier\b|\bslc\d+",
    "receptor": r"\breceptor\b",
    "collagen": r"\bcollagen\b|\bcol\d",
    "integrin": r"\bintegrin\b|\bitg[a-z]",
    "cytokine": r"\bcytokine\b|\binterleukin\b|\bil[- ]?\d+\b|^il\d+",
    "chemokine": r"\bchemokine\b|\bcxcl\b|\bccl\b|\bxlr\b",
    "histone": r"\bhistone\b|\bh[1234]\.\w",
    "ubiquitin": r"\bubiquitin\b|\bubiq",
    "proteasome": r"\bproteasome\b|\bpsm([abdg]|c)",
    "splicing": r"\bsplice|\bsplicing|\bsnrp|\bsf3\b|\bu2af\b",
    "ribonucl": r"\bribonucleoprotein\b",
    "cytochrome": r"\bcytochrome\b|\bcox\b|\bmt-co",
}

# Count occurrences per category
counts: dict[str, int] = {}
hits_idx: dict[str, list[int]] = {k: [] for k in patterns}
for i, s in enumerate(text):
    for k, rx in patterns.items():
        if re.search(rx, s):
            counts[k] = counts.get(k, 0) + 1
            hits_idx[k].append(i)

# Build summary table
summary_rows: list[dict[str, int]] = []
for k in sorted(counts, key=counts.get, reverse=True):
    summary_rows.append({"category_keyword": k, "n_genes": counts[k]})

summary_df: pd.DataFrame = pd.DataFrame(summary_rows)

# Collect examples for top categories
examples: dict[str, pd.DataFrame] = {}
for k in counts:
    idxs = hits_idx[k][:10]
    examples[k] = df.iloc[idxs, [0,2]].rename(columns={df.columns[0]:"gene_symbol", df.columns[2]: "annotation"})

# Save results
summary_path = os.path.join(outdir, "keyword_category_counts.tsv")
summary_df.to_csv(summary_path, sep="\t", index=False)

# Also build a long table: for each matched category, list genes
long_rows: list[dict[str, str]] = []
for k, idxs in hits_idx.items():
    for i in idxs:
        long_rows.append({
            "category_keyword": k,
            "gene_symbol": str(df.iloc[i,0]),
            "annotation": str(df.iloc[i,2])
        })
long_df: pd.DataFrame = pd.DataFrame(long_rows)
long_path = os.path.join(outdir, "keyword_category_gene_list.tsv")
long_df.to_csv(long_path, sep="\t", index=False)

# Print small preview tables to user
print("\nKeyword category summary (top 20 rows):")
print(summary_df.head(20).to_string(index=False))

print("\nKeyword matches (first 200 rows):")
print(long_df.head(200).to_string(index=False))

print(f"\nAnalyzed column: {colname}")
print("Top categories:", summary_df.head(10).to_dict(orient='records'))
print(f"Saved summary: {summary_path}")
print(f"Saved matches: {long_path}")
print("Sample third-column values:", sample_vals)
