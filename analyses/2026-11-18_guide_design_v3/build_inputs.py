"""Build the element BED and the positive-control genes file for the v2 guide design.

Starts from the 2026-10-13 element set and the 2025-10-24 positive pairs and
adds, on GENCODE v43:

- the three lifted elements that are in no day peak file;
- a C2CD6 promoter element, 500 bp ending at the safe-targeting guide's
  protospacer and covering the MANE TSS, so that guide is a promoter guide;
- the gene each design promoter was designed against, for every design
  promoter the old pairs file left without one and the genome calls a
  promoter (within 2 kb of a protein-coding TSS);
- the gene of each promoter a safe-targeting guide was recovered into
  (SRPK1, FJX1, ZNF609) and of the chr17 lifted element (NUFIP2);
- the gene of every validated enhancer, from the lab's table of pairs
  validated in earlier screens, which the old pairs file covered for 12 of 26.
"""

import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OLD = REPO / "results" / "03-2025-10-24_create_sceptre_input"
ELEMENTS = REPO / "data" / "design_peaks" / "elements_combined.bed"
OLD_EXPORT = Path.home() / "GitHub/crispr-analysis-utils/results/hiPSC-to-EC_guides/igvf/guide_metadata.tsv"
STORE = Path.home() / "Annotations/HG38/GENCODE/gencode.v43.annotation.parquet"
VALIDATED = REPO / "metadata" / "legacy" / "validate_pairs_from_previous_studies_final.tsv"
OUT = REPO / "results" / "2026-11-12_guide_design_v2" / "inputs"

LIFTED = ["chr12:114126970-114129885", "chr12:114148790-114150612", "chr17:29293427-29293781"]
C2CD6_GUIDE_END = 201_619_432
RECOVERED = {  # element id in elements_combined.bed -> gene
    "chr6:35920847-35921147": "SRPK1",
    "chr11:35617933-35618233": "FJX1",
    "chr15:64460290-64460722": "ZNF609",
    "chr17:29293427-29293781": "NUFIP2",
}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    store = pd.read_parquet(STORE)
    genes = store[store.Feature == "gene"]
    ensg_of = genes.groupby("gene_symbol").ensg.agg(list)

    def one_ensg(symbol: str) -> str:
        hits = ensg_of.get(symbol, [])
        if len(hits) != 1:
            raise SystemExit(f"{symbol}: {len(hits)} GENCODE v43 genes")
        return hits[0]

    columns = ["chr", "start", "end", "element_id", "category", "days"]
    elements = pd.read_csv(ELEMENTS, sep="\t", header=None, names=columns, keep_default_na=False)

    old = pd.read_csv(OLD / "all_element_annotations_final.bed", sep="\t", header=None,
                      names=["chr", "start", "end", "element_id", "category"])
    lifted = old[old.element_id.isin(LIFTED)].assign(days="")
    if len(lifted) != len(LIFTED):
        raise SystemExit("a lifted element is missing from the old BED")

    c2cd6 = store[(store.gene_symbol == "C2CD6") & (store.Feature == "transcript") & store.is_reference]
    if len(c2cd6) != 1 or c2cd6.Strand.iloc[0] != "-":
        raise SystemExit("expected one reference C2CD6 transcript on the minus strand")
    tss = int(c2cd6.End.iloc[0]) - 1
    start, end = C2CD6_GUIDE_END - 500, C2CD6_GUIDE_END
    if not start <= tss < end:
        raise SystemExit("the C2CD6 element does not cover its TSS")
    c2cd6_row = pd.DataFrame([{"chr": "chr2", "start": start, "end": end,
                               "element_id": f"chr2:{start}-{end}", "category": "promoter", "days": ""}])
    RECOVERED[f"chr2:{start}-{end}"] = "C2CD6"

    combined = pd.concat([elements, lifted[columns], c2cd6_row], ignore_index=True)
    combined = combined.sort_values(["chr", "start", "end"]).reset_index(drop=True)
    if combined.element_id.duplicated().any():
        raise SystemExit("duplicate element ids")
    combined.to_csv(OUT / "elements_combined_v2.bed", sep="\t", header=False, index=False)

    # Design promoters without a pair, and what the genome calls them.
    pairs = pd.read_csv(OLD / "positive_control_pairs_dedup_final.tsv", sep="\t")
    paired = set(pairs.grna_target)
    promoters = elements[elements.category == "promoter"]

    def keys(row) -> set[str]:
        return {f"{row.chr}:{row.start}-{row.end}", f"{row.chr}:{row.start + 1}-{row.end}", row.element_id}

    missing = promoters[[not (keys(row) & paired) for row in promoters.itertuples()]]
    design = pd.read_csv(OLD / "promoters_overlapping_guides_final.tsv", sep="\t")
    design["key"] = design.element_chromosome + ":" + design.element_start.astype(str) + "-" + design.element_end.astype(str)
    alias = design.groupby("key").element_alias_list.first()

    export = pd.read_csv(OLD_EXPORT, sep="\t", dtype=str, keep_default_na=False)
    export["key"] = export.intended_target_chr + ":" + export.intended_target_start + "-" + export.intended_target_end
    genome_class = export.groupby("key").genomic_element.agg(lambda x: ",".join(sorted(set(x))))

    added, report = [], []
    for row in missing.itertuples():
        key = f"{row.chr}:{row.start}-{row.end}"
        cls = genome_class.get(key, "")
        found = [k for k in keys(row) if k in alias.index]
        if row.element_id in RECOVERED:
            symbol, source = RECOVERED[row.element_id], "recovered_promoter"
        elif found:
            symbol = re.sub(r"_x.*$", "", re.sub(r"_TSS(_\d+)?$", "", alias[found[0]]))
            source = "design_promoter"
        else:
            raise SystemExit(f"{row.element_id}: no design gene")
        report.append({"element_id": row.element_id, "gene_symbol": symbol, "genome_class": cls, "source": source})
        if cls == "promoter":
            added.append({"grna_target": row.element_id, "response_id": one_ensg(symbol),
                          "gene_symbol": symbol, "element_source": source})

    for element_id in ("chr17:29293427-29293781", f"chr2:{start}-{end}"):
        symbol = RECOVERED[element_id]
        added.append({"grna_target": element_id, "response_id": one_ensg(symbol),
                      "gene_symbol": symbol, "element_source": "recovered_promoter"})
        report.append({"element_id": element_id, "gene_symbol": symbol, "genome_class": "new", "source": "recovered_promoter"})

    validated = pd.read_csv(VALIDATED, sep="\t")
    e2g = set(elements.loc[elements.category == "e2g_validated", "element_id"])
    if set(validated.HG38) != e2g:
        raise SystemExit("the validated pairs and the e2g_validated elements disagree")
    for row in validated.itertuples():
        added.append({"grna_target": row.HG38, "response_id": one_ensg(row.Response_ID),
                      "gene_symbol": row.Response_ID, "element_source": "e2g_validated"})

    out = pd.concat([pairs, pd.DataFrame(added)], ignore_index=True).drop_duplicates(["grna_target", "response_id"])
    out.to_csv(OUT / "positive_control_genes_v2.tsv", sep="\t", index=False)
    pd.DataFrame(report).to_csv(OUT / "promoters_added.tsv", sep="\t", index=False)
    print(f"elements {len(elements)} -> {len(combined)}; genes rows {len(pairs)} -> {len(out)}; "
          f"design promoters checked {len(missing)}, added {len(added)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
