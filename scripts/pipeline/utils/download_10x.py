"""
Modified from download_files.py by J.Pham

"""

import os
import re
import sys
import argparse
import logging
import subprocess
from tabulate import tabulate
from distutils.version import StrictVersion

#import scgenome.db.qc
from dbclients.tantalus import TantalusApi
from dbclients.colossus import ColossusApi
from dbclients.basicclient import NotFoundError

tantalus_api = TantalusApi()
colossus_api = ColossusApi()


# + container='scrna_rdatarawv3' : raw counts data mapped to human genome
# + container='scrna_rdatav3' : filtered counts data using QC script from single cell pipeline, mapped to human genome
# + container='scrna_rdatarawmousev3': raw counts data mapped to mouse genome
# + container='scrna_rdatamousev3': filtered counts data using QC script from single cell pipeline, mapped to mouse genome
def download_data(data_dir, library, is_mouse):
    # init directory to download data to
    data_dir = os.path.join(data_dir, library)
    # check if destination path exists
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # init storage client
    if is_mouse == True:
        print('Downloading mouse alignment output')
        print(library)
        storage_client = tantalus_api.get_storage_client("scrna_rdatarawmousev3")
    else:
        print('Downloading human alignment output')
        print(library)
        storage_client = tantalus_api.get_storage_client("scrna_rdatarawv3")

    # list all blobs for library
    blobs = storage_client.list(library)

    for blob in blobs:
        # get flowcell from path
        flowcell = os.path.basename(os.path.dirname(blob))
        # get fastq filename
        filename = os.path.basename(blob)
        print(filename)

        # join destination path with flowcell name and create path
        flowcell_path = os.path.join(data_dir, flowcell)
        if not os.path.exists(flowcell_path):
            os.makedirs(flowcell_path)

        # format filepath
        filepath = os.path.join(flowcell_path, filename)

        # download blob to path
        print(f"downloading {blob} to {filepath}")
#         blob = storage_client.blob_service.get_blob_to_path(container_name="rnaseq", blob_name=blob, file_path=filepath)
        blob_client = storage_client.blob_service.get_blob_client("rnaseq", blob)  # download bam files: replace rnaseq by bams

        # if(is_mouse == True):
        #     print('Downloading mouse alignment output')
        #     print(library)
        #     blob_client = storage_client.blob_service.get_blob_client("rnaseq", blob)  # download bam files: replace rnaseq by bams
        # else:
        #     print('Downloading human alignment output')
        #     print(library)
        #     blob_client = storage_client.blob_service.get_blob_client("rnaseq", blob)  # download bam files: replace rnaseq by bams
        with open(filepath, "wb") as my_blob:
            download_stream = blob_client.download_blob()
            my_blob.write(download_stream.readall())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--library_id', required=True)
    parser.add_argument('--is_mouse', required=True)
    parser.add_argument('--download_dir', required=True) #, default='.'
    args = vars(parser.parse_args())
    print("Download dir:" + args['download_dir'])
    print("Library:" + args['library_id'])
    print("Is mouse:" + args['is_mouse'])
    if args['is_mouse']=='True':
        is_mouse = True
    else:
        is_mouse = False
    
    download_data(args['download_dir'], args['library_id'], is_mouse)
