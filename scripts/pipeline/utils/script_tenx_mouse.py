import pandas as pd
import os
import azure
import sys
import azure.storage.blob
from dbclients.basicclient import NotFoundError
from dbclients.tantalus import TantalusApi, get_storage_account_key
from datamanagement.transfer_files import TransferProgress
import argparse

tantalus_api = TantalusApi()


def get_azure_blob_service(storage_account):
    client_id = os.environ["CLIENT_ID"]
    secret_key = os.environ["SECRET_KEY"]
    tenant_id = os.environ["TENANT_ID"]
    keyvault_account = os.environ['AZURE_KEYVAULT_ACCOUNT']

    storage_key = get_storage_account_key(
            storage_account,
            client_id,
            secret_key,
            tenant_id,
            keyvault_account)

    blob_service = azure.storage.blob.BlockBlobService(
        account_name=storage_account,
        account_key=storage_key)
    blob_service.MAX_BLOCK_SIZE = 64 * 1024 * 1024

    return blob_service

def main(library_id, cache_folder, container='rdatarawmousev3'):
    print('From container: {0}'.format(container))
    blob_service = get_azure_blob_service("scrnadata")
    block_blob_service = tantalus_api.get_storage_client("scrna_"+container).blob_service

    
    if not os.path.exists(os.path.join(cache_folder, library_id)):
        os.mkdir(os.path.join(cache_folder, library_id))

    #get all the files inside the rdatarawv3 container
    blobs_all = blob_service.list_blobs(container)

    #iterate over the blobs
    for blob in blobs_all:
        #check if the filtering condition is satisfied
        if library_id in blob.name:
        # if (library_id in blob.name) and ("mouse" not in blob.name.lower()):
            print("Downloading: ", blob.name)
            #download the file into the directory
            block_blob_service.get_blob_to_path(
                container,
                blob.name,
                os.path.join(cache_folder, library_id, blob.name),
                progress_callback=TransferProgress().print_progress,
                max_connections=16,
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--library_id', required=True)
    parser.add_argument('--download_dir', default='downloads')

    args = vars(parser.parse_args())
    print(args['download_dir'])
    print(args['library_id'])
    main(args['library_id'], args['download_dir'],'rdatarawmousev3')
    # library_id = sys.argv[1]
    # main(library_id, cache_folder)
