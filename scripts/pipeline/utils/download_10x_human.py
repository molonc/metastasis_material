import os
# import click
import subprocess
import sys
import argparse
from dbclients.tantalus import TantalusApi

tantalus_api = TantalusApi()


# @click.command()
# @click.argument("data_dir")
# @click.argument("library")
def download_data(data_dir, library, is_mouse):
    # init directory to download data to
    data_dir = os.path.join(data_dir, library)
    # check if destination path exists
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # init storage client
    # if not is_mouse:
    #     storage_client = tantalus_api.get_storage_client("scrna_rdatarawv3") #scrna_rdatarawmousev3
    #     print('Downloading data from human alignment storage')
    # else:
    #     storage_client = tantalus_api.get_storage_client("scrna_rdatarawmousev3") #
    #     print('Downloading data from mouse alignment storage')
    storage_client = tantalus_api.get_storage_client("scrna_rdatarawmousev3") #
    # list all blobs for library
    blobs = storage_client.list(library)

    for blob in blobs:
        # get flowcell from path
        flowcell = os.path.basename(os.path.dirname(blob))
        # get fastq filename
        filename = os.path.basename(blob)

        # join destination path with flowcell name and create path
        flowcell_path = os.path.join(data_dir, flowcell)
        if not os.path.exists(flowcell_path):
            os.makedirs(flowcell_path)

        # format filepath
        filepath = os.path.join(flowcell_path, filename)
        # check if file already exists with same size from blob storage
        if os.path.exists(filepath) and os.path.getsize(filepath) == storage_client.get_size(blob):
            continue

        # download blob to path
        print(f"downloading {blob} to {filepath}")
        #blob = storage_client.blob_service.get_blob_to_path(container_name="rnaseq", blob_name=blob, file_path=filepath)
        # blob_client = storage_client.blob_service.get_blob_client("rnaseq", blob)
        blob_client = storage_client.blob_service.get_blob_client("rdatarawv3", blob)
        with open(filepath, "wb") as my_blob:
            download_stream = blob_client.download_blob()
            my_blob.write(download_stream.readall())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--library_id', required=True)
    parser.add_argument('--is_mouse', required=False, default='F')
    parser.add_argument('--download_dir', required=True) #, default='.'
    args = vars(parser.parse_args())
    print("Download dir:" + args['download_dir'])
    print("Library:" + args['library_id'])
    print("Is mouse:" + args['is_mouse'])
    if args['is_mouse']=='F':
        is_mouse = False
    else:
        is_mouse = True
        
    
    download_data(args['download_dir'], args['library_id'], is_mouse)


# python Projects/farhia_project/rnaseq/pipeline/utils/download_10x_human.py --library_id SCRNA10X_SA_CHIP0071_000 --download_dir /home/htran/storage/rnaseq_datasets/drug_resistance_RNAseq/testing/
# python ~/Projects/farhia_project/rnaseq/pipeline/utils/download_10x_human.py --library_id TENX_015_SA604X6XB01979_001 --download_dir /home/htran/storage/rnaseq_datasets/drug_resistance_RNAseq/added_pdx_human/
# python ~/Projects/farhia_project/rnaseq/pipeline/utils/download_10x_human.py --library_id TENX_015_SA604X6XB01979_001 --download_dir /home/htran/storage/rnaseq_datasets/drug_resistance_RNAseq/added_pdx_mouse/
