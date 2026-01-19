# 10X Multiome 5 Timepoints

## GEX


## ATAC


## CMO analysis

For the CMO analysis I haven't run yet the counts, I used the output that was already generated.

The paths on oak are the following:
 - Channel 1: /oak/stanford/groups/engreitz/Users/dulguun/220519_EC_Profiling/flip_index/CMO1_CITE-seq-Count_outputs
 - Channel 2: /oak/stanford/groups/engreitz/Users/dulguun/220519_EC_Profiling/flip_index/CMO2_CITE-seq-Count_outputs

There are multiple lanes per channel, so they need to be summed up. I also found out that when running CITE-seq-Count, the whitelist was set to the ATAC barcodes passing QC filters and converted to GEX barcodes.

The UMI counts present in the `umi_count` folder are not in a format that can be directly read by `scanpy`. 

The script used to format the CMO metadata is [format_cmo_metadata.py](format_cmo_metadata.py). I fixed the barcode names to make them matche the GEX barcodes from the IGVF single-cell pipeline.
I also added all the lanes together for each channel and annotate the h5ad files with the best CMO per cell. I am also reporting the second best CMO.

The following columns are added to the `obs` dataframe:
 - obs["cmo_best"]
 - obs["cmo_best_count"]
 - obs["cmo_second"]
 - obs["cmo_second_count"]
 - obs["cmo_margin"]

The final output files can be found here:
 - `gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/data/10x_multiome_5_timepoints/cmo_channel1.annotated.h5ad`
 - `gs://fc-76565551-28c3-4b6d-a048-e272103bcbd1/data/10x_multiome_5_timepoints/cmo_channel2.annotated.h5ad`

