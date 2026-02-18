# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "alphagenome>=0.6.0",
#     "biopython>=1.86",
#     "marimo>=0.19.10",
#     "numpy>=2.4.2",
#     "pandas>=3.0.1",
#     "pyfaidx>=0.9.0.3",
#     "pyzmq>=27.1.0",
# ]
# ///

import marimo

__generated_with = "0.19.11"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return


@app.cell
def _():
    import os
    import re
    import io
    import hashlib
    import numpy as np
    import pandas as pd
    from random import shuffle


    from pyfaidx import Faidx
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    return Faidx, dna_client, genome, io, pd, re, shuffle


@app.cell
def _(Faidx, re):
    hg38_fa_fnp = "/Users/emattei/Annotations/HG38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.bgz"

    hg38_fa = Faidx(hg38_fa_fnp)

    wtc11 = "EFO:0009747"
    mesoderm = "CL:0000222"
    endothelial = "CL:0000115"

    ontology_terms = [wtc11, mesoderm, endothelial]

    # Replace this with reading from a file if you want:
    CSV_TEXT = """response_id,grna_target,n_nonzero_trt,n_nonzero_cntrl,pass_qc,p_value,fold_change,se_fold_change,log_2_fold_change,significant
    KIT,chr4:54913795-54914136,200,412620,TRUE,0.00000000000000000000000000000000000000000000078823038618271,0.438489928902994,0.0332811381191506,-1.18938438719697,TRUE
    HMGA1,chr6:34235537-34235838,461,516339,TRUE,0.0,0.426837679748651,0.0164043807107939,-1.22824055699335,TRUE
    RBPMS2,chr15:64774621-64774922,354,512645,TRUE,0.0,0.450306664921255,0.0213601665582408,-1.15102026393715,TRUE
    RAB11A,chr15:65870400-65870701,467,515811,TRUE,0.0,0.564324565117943,0.0154053788290388,-0.825402943245076,TRUE
    RBPMS2,chr15:64773978-64774279,430,512569,TRUE,0.0,0.522459918499036,0.0173659186063174,-0.936607732527732,TRUE
    CDC37,chr19:10402436-10402737,405,510161,TRUE,0.0,0.606513526351769,0.0179339508532012,-0.721388274466534,TRUE
    ZNF827,chr4:145937620-145937921,126,342552,TRUE,0.0,0.339787538725481,0.0308165178004764,-1.55729515039526,TRUE
    EPOR,chr19:11382772-11383073,224,404774,TRUE,0.0,0.424137778109223,0.0270796089165846,-1.23739510484741,TRUE
    RPL19,chr17:39197453-39197953,358,516444,TRUE,0.0,0.921390316994809,0.00683494461530867,-0.118115658326804,TRUE
    CARM1,chr19:10872314-10872623,334,433994,TRUE,0.0,0.60181103700506,0.0269436381628289,-0.732617529448326,TRUE
    """


    def parse_interval(s: str):
        m = re.match(r"^(chr[\w]+):(\d+)-(\d+)$", s)
        if not m:
            raise ValueError(f"Bad interval string: {s}")
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
        if end <= start:
            raise ValueError(f"End <= start in: {s}")
        return chrom, start, end


    return CSV_TEXT, hg38_fa, ontology_terms, parse_interval


@app.cell
def _(
    CSV_TEXT,
    dna_client,
    genome,
    hg38_fa,
    io,
    ontology_terms,
    parse_interval,
    pd,
    shuffle,
):
    api_key = ""
    if not api_key:
        raise RuntimeError("Missing ALPHAGENOME_API_KEY env var.")

    model = dna_client.create(api_key)

    df = pd.read_csv(io.StringIO(CSV_TEXT))
    df["pass_qc"] = df["pass_qc"].astype(str).str.upper().eq("TRUE")

    out_rows = []

    for _, row in df.iterrows():
        if not bool(row["pass_qc"]):
            continue

        gene = str(row["response_id"])
        region = str(row["grna_target"])
        crispri_log2fc = float(row["log_2_fold_change"])

        chrom, start, end = parse_interval(region)

        query_seq = hg38_fa.fetch(chrom, start, end)
        query_seq_list = list(query_seq.seq)
        shuffle(query_seq_list)
        query_seq_scrambled = ''.join(query_seq_list)
        query_seq_scrambled

        query_interval = genome.Interval(chromosome=chrom, start=start, end=end)
        query_interval_resized = query_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

        #out_rows.append(model.predict_interval(
        #    interval=query_interval_resized,
        #    organism=dna_client.Organism.HOMO_SAPIENS,
        #    requested_outputs=[dna_client.OutputType.RNA_SEQ],
        #    ontology_terms=[wtc11, mesoderm, endothelial]
        #))

        query_variant = genome.Variant(
            chromosome=chrom,
            position=start,
            reference_bases=query_seq.seq,
            alternate_bases=query_seq_scrambled,
            name=region
    
        )
        # Make predictions for sequences containing the REF and ALT alleles.
        output = model.predict_variant(
            interval=query_interval_resized,
            variant=query_variant,
            requested_outputs={
                dna_client.OutputType.RNA_SEQ,
            },
            ontology_terms=ontology_terms,
        )

        break
    output

    return output, query_interval_resized


@app.cell
def _(pd):
    from alphagenome.visualization import plot_components
    from alphagenome.data import gene_annotation, track_data, transcript



    # Load gene annotations (from GENCODE).
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )

    # Filter to protein-coding genes and highly supported transcripts.
    gtf_transcript = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1']
    )

    # Extractor for identifying transcripts in a region.
    transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)

    # Also define an extractor that fetches only the longest transcript per gene.
    gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(
        gtf_transcript
    )

    longest_transcript_extractor = transcript.TranscriptExtractor(
        gtf_longest_transcript
    )

    return longest_transcript_extractor, plot_components


@app.cell
def _(longest_transcript_extractor, query_interval_resized):
    longest_transcripts = longest_transcript_extractor.extract(query_interval_resized)
    return (longest_transcripts,)


@app.cell
def _(longest_transcripts, output, plot_components, query_interval_resized):
    # Build plot.
    plot_reference = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.Tracks(
                tdata=output.reference.rna_seq,
                ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
            ),
        ],
        interval=query_interval_resized,
        title='Predicted RNA Expression (RNA_SEQ) - Reference Sequence',
    )
    plot_reference
    return


@app.cell
def _(longest_transcripts, output, plot_components, query_interval_resized):
    # Build plot.
    plot_alternate = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.Tracks(
                tdata=output.alternate.rna_seq.filter_to_positive_strand(),
                ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
            ),
        ],
        interval=query_interval_resized,
        title='Predicted RNA Expression (RNA_SEQ) - Alternate Sequence',
    )
    plot_alternate
    return


@app.cell
def _(longest_transcripts, output, plot_components, query_interval_resized):

    plot2 = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            # RNA-seq tracks.
            plot_components.Tracks(
                tdata=output.alternate.rna_seq.filter_to_positive_strand() - output.reference.rna_seq.filter_to_positive_strand(),
                ylabel_template='{biosample_name} ({strand})\n{name}',
                filled=True,
            ),
        ],
        #annotations=[plot_components.VariantAnnotation([query_variant])],
        interval=query_interval_resized,
        #title=(
        #    'Effect of variant on predicted RNA Expression, DNAse, and ChIP-Histone'
        #    f' in CD34 positive HSC.\n{query_variant=}'
        #),
    )
    plot2
    return


if __name__ == "__main__":
    app.run()
