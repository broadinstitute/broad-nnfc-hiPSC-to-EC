import marimo

__generated_with = "0.13.15"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    x_overlapping_guides = set()
    with open("results/20250613_align_guides/K562_guides_not_overlapping_dnase_peaks_hg38_names.txt", "r") as f:
        for _line in f:
            x_overlapping_guides.add(_line.strip())


    return (x_overlapping_guides,)


@app.cell
def _(x_overlapping_guides):
    x_overlapping_guides
    return


@app.cell
def _(countinue, x_overlapping_guides):
    missed_genes = set()
    count_missing_pairs = 0
    with open("results/20250613_align_guides/hg38_significant.txt", "r") as f:
        for _line in f:
            gene_id, guides = _line.strip().split("\t")
            for guide in guides.split(","):
                if guide in x_overlapping_guides:
                    missed_genes.add(gene_id)
                    count_missing_pairs += 1
                    countinue
        
    return


if __name__ == "__main__":
    app.run()
