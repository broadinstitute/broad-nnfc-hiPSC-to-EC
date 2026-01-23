# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "ipython>=9.9.0",
#     "marimo>=0.19.0",
#     "pyzmq>=27.1.0",
#     "snapatac2>=2.8.0",
# ]
# ///

import marimo

__generated_with = "0.19.4"
app = marimo.App()


@app.cell
def _():
    import snapatac2 as snap

    return (snap,)


@app.cell
def _(snap):
    data = snap.pp.import_fragments(
        "../../data/external/10x_multiome_5_timepoints/atac/10x_5timepoints_channel1.fragments.tsv.gz",
        chrom_sizes=snap.genome.hg38,
        sorted_by_barcode=False,
    )
    data
    return (data,)


@app.cell
def _(data, snap):
    snap.metrics.tsse(data, snap.genome.hg38)

    return


@app.cell
def _(data):
    data.obs.tsse
    return


if __name__ == "__main__":
    app.run()
