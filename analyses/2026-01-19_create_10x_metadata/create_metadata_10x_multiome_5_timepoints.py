# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "google-cloud-storage",
#     "marimo>=0.19.0",
#     "pandas>=2.3.3",
#     "protobuf>=6.33.4",
#     "pyzmq>=27.1.0",
# ]
# ///

import marimo

__generated_with = "0.19.4"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    from google.cloud import storage
    import pandas as pd
    return pd, storage


@app.cell
def _(mo):
    mo.md(r"""
    channel1 =
    GEX_2,GEX_1(with extra sequencing),ATAC_Lib1

    channel2 =
    GEX_1(with extra sequencing),GEX_2,ATAC_Lib2

    Results in table format:
    sample_id,GEX_R1, GEX_R2, ATAC_R1, ATAC_R2, CMO_R1, CMO_R2, channel
    WTC11_channel1,
    """)
    return


@app.cell
def _(storage):
    def list_gcs_objects(bucket_name):
        """Lists all the objects in a given Google Cloud Storage bucket."""
        # Instantiates a client
        storage_client = storage.Client()

        # Get the bucket
        bucket = storage_client.bucket(bucket_name)

        # Lists all blobs in the bucket
        blobs = bucket.list_blobs()

        print(f"Objects in bucket '{bucket_name}':")
        # Iterate over the blobs and print their names
        for blob in blobs:
            print(f"* {blob.name}")

        # You can return the list of names if needed
        return [blob.name for blob in blobs]

    def list_folder_contents(bucket_name, folder_name):
        """Lists all the blobs in a specific folder within a bucket."""
        storage_client = storage.Client()
        # Ensure the folder name ends with a '/' for correct prefix matching
        prefix = folder_name if folder_name.endswith('/') else folder_name + '/'

        blobs = storage_client.list_blobs(bucket_name, prefix=prefix)

        print(f"Files in folder {folder_name}:")
        for blob in blobs:
            print(blob.name)
    return


@app.cell
def _():
    data_bucket = "fc-76565551-28c3-4b6d-a048-e272103bcbd1"

    # GEX library
    gex = "data/10x_multiome_5_timepoints/H2YL5BGXM/"
    # GEX extra sequencing
    gex_extra = "data/10x_multiome_5_timepoints/HF7JFBGXM/"
    # ATAC library
    atac_lib1 = "data/10x_multiome_5_timepoints/HFKNCBGXL/ATAC_Lib_1"
    atac_lib2 = "data/10x_multiome_5_timepoints/HFKNCBGXL/ATAC_Lib_2"
    return atac_lib1, atac_lib2, data_bucket, gex, gex_extra


@app.cell
def _(atac_lib1, atac_lib2, data_bucket, gex, gex_extra, pd, storage):
    import re
    from collections import defaultdict

    def filter_r1_r2_files(file_list):
        """Filter files to only include R1 and R2 reads."""
        return [f for f in file_list if re.search(r'_R[12]_\d+\.fastq\.gz$', f)]

    def filter_atac_files(file_list):
        """Filter files to only include R1 and R3 reads for ATAC."""
        return [f for f in file_list if re.search(r'_R[13]_\d+\.fastq\.gz$', f)]

    def group_files_by_sample_and_read(file_list, read_types=['R1', 'R2']):
        """Group files by sample ID and read type."""
        groups = defaultdict(lambda: {rt: [] for rt in read_types})
    
        for file in file_list:
            match = re.search(r'(GEX_\d+|ATAC_Lib_\d+).*_(R[123])_\d+\.fastq\.gz$', file)
            if match:
                sample_id = match.group(1)
                read_type = match.group(2)
                if read_type in read_types:
                    groups[sample_id][read_type].append(file)
    
        return groups

    # Get file lists and filter for R1/R2
    storage_client = storage.Client()
    gex_blobs = list(storage_client.list_blobs(data_bucket, prefix=gex))
    gex_files = [blob.name for blob in gex_blobs]

    gex_extra_blobs = list(storage_client.list_blobs(data_bucket, prefix=gex_extra))
    gex_extra_files = [blob.name for blob in gex_extra_blobs]

    atac_lib1_blobs = list(storage_client.list_blobs(data_bucket, prefix=atac_lib1))
    atac_lib1_files = [blob.name for blob in atac_lib1_blobs]

    atac_lib2_blobs = list(storage_client.list_blobs(data_bucket, prefix=atac_lib2))
    atac_lib2_files = [blob.name for blob in atac_lib2_blobs]

    # Filter and group GEX
    gex_r1_r2 = filter_r1_r2_files(gex_files)
    gex_extra_r1_r2 = filter_r1_r2_files(gex_extra_files)

    gex_groups = group_files_by_sample_and_read(gex_r1_r2, ['R1', 'R2'])
    gex_extra_groups = group_files_by_sample_and_read(gex_extra_r1_r2, ['R1', 'R2'])

    # Combine GEX_1 files from both folders
    if 'GEX_1' in gex_extra_groups:
        gex_groups['GEX_1']['R1'].extend(gex_extra_groups['GEX_1']['R1'])
        gex_groups['GEX_1']['R2'].extend(gex_extra_groups['GEX_1']['R2'])

    # Filter and group ATAC
    atac_lib1_r1_r3 = filter_atac_files(atac_lib1_files)
    atac_lib2_r1_r3 = filter_atac_files(atac_lib2_files)

    atac_lib1_groups = group_files_by_sample_and_read(atac_lib1_r1_r3, ['R1', 'R3'])
    atac_lib2_groups = group_files_by_sample_and_read(atac_lib2_r1_r3, ['R1', 'R3'])

    # Combine ATAC groups
    atac_groups = {**atac_lib1_groups, **atac_lib2_groups}

    # Sort files to ensure consistent lane order
    for sample in gex_groups:
        gex_groups[sample]['R1'].sort()
        gex_groups[sample]['R2'].sort()

    for sample in atac_groups:
        atac_groups[sample]['R1'].sort()
        atac_groups[sample]['R3'].sort()

    # Create DataFrame with combined GEX and ATAC
    rows = []

    # Define the pairing: GEX_1 -> ATAC_Lib_2, GEX_2 -> ATAC_Lib_1
    gex_atac_pairing = {
        'GEX_1': 'ATAC_Lib_2',
        'GEX_2': 'ATAC_Lib_1'
    }

    for gex_sample, atac_sample in gex_atac_pairing.items():
        # Get GEX files
        gex_r1 = str([f'gs://{data_bucket}/{f}' for f in gex_groups[gex_sample]['R1']]).replace("'", '"')
        gex_r2 = str([f'gs://{data_bucket}/{f}' for f in gex_groups[gex_sample]['R2']]).replace("'", '"')
    
        # Get ATAC files
        atac_r1 = str([f'gs://{data_bucket}/{f}' for f in atac_groups[atac_sample]['R1']]).replace("'", '"')
        atac_r3 = str([f'gs://{data_bucket}/{f}' for f in atac_groups[atac_sample]['R3']]).replace("'", '"')
    
        rows.append({
            'entity:10x_multiome_5_timepoints_id': gex_sample,
            'GEX_R1': gex_r1,
            'GEX_R2': gex_r2,
            'gDNA_R1': atac_r1,
            'gDNA_R3': atac_r3
        })

    df = pd.DataFrame(rows)

    # Save to CSV
    output_file = 'multiome_metadata.tsv'
    df.to_csv(output_file, index=False, sep="\t")
    print(f"Saved metadata to {output_file}")
    print("\nDataFrame:")
    df
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
