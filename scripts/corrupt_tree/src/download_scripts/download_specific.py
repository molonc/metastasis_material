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

from dbclients.tantalus import TantalusApi
from dbclients.colossus import ColossusApi
from dbclients.basicclient import NotFoundError

tantalus_api = TantalusApi()
colossus_api = ColossusApi()

JIRA_TICKET_REGEX = r".*SC-\d{4}$"

def find_latest_analysis_ticket(pool_id):
    '''
    Return the jira ticket of the newest possible analysis of the given library.
    '''
    analyses = list(colossus_api.list("analysis_information", library__pool_id = pool_id))
    latest_completion_version = None
    analysis_jira_ticket = None
    newest_version = False
    #loop over the list of analyses to find the latest possible one
    for analysis in analyses:
        status = analysis["analysis_run"]["run_status"]
        # print(status)
        current_version = analysis["version"][1:]
        # print(current_version)
        if status in ["complete", "hmmcopy_complete"]:
            if latest_completion_version:
                if StrictVersion(current_version.strip('v')) >= StrictVersion(latest_completion_version.strip('v')):
                    #print("The updated version %s" % (current_version))
                    latest_completion_version = current_version
                    analysis_jira_ticket = analysis["analysis_jira_ticket"]
                    current_analysis = analysis
            else:
                latest_completion_version = current_version
                analysis_jira_ticket = analysis["analysis_jira_ticket"]
                current_analysis = analysis
                newest_version = True
            #print(status, current_version, latest_completion_version, current_version > latest_completion_version)
        else:
            pass
            #if the latest COMPLETE analysis exisis, save the information to the dict.

    if newest_version:
        return analysis_jira_ticket

    return None

def check_user_input(intial_question, options):

    answer = input(intial_question)

    while answer not in options:
        answer = input("Invalid input. Options are {}: ".format(options))

        if answer in options:
            break

    return answer


def check_user_file_input(user_input):
    file_choices = ["metrics", "reads", "segments", "plots"]
    if user_input.lower() == "y":
        return "all"

    elif user_input.lower() == "n":
        file_table = tabulate([[f] for f in file_choices], showindex="always")
        file_question = "\n{}\n\nPlease enter the index of the desired file: ".format(file_table)

        file_index = check_user_input(file_question, list(map(str, range(len(file_choices)))))

    print("\nPreparing to download {} files...\n".format(file_choices[int(file_index)]))
    return file_choices[int(file_index)]


def get_reads(jira_ticket, version, library_id):

    filenames = []
    if StrictVersion(version.strip('v')) < StrictVersion('0.2.23'):
        file_resources = list(tantalus_api.list(
            "file_resource",
            filename__startswith=jira_ticket,
            filename__endswith="{}_multiplier0_reads.csv.gz".format(library_id)
        ))
    else:
        file_resources = list(tantalus_api.list(
            "file_resource",
            filename__startswith=jira_ticket,
            filename__endswith="{}_reads.csv.gz".format(library_id)
        ))

    filenames += [f["filename"] for f in file_resources]

    return filenames

def get_segment(jira_ticket, version, library_id):

    filenames = []
    file_resources = list(tantalus_api.list(
            "file_resource",
            filename__startswith=jira_ticket,
            filename__endswith="{}_segments.csv.gz".format(library_id)
        ))

    filenames += [f["filename"] for f in file_resources]

    return filenames

def get_metrics(jira_ticket, version, library_id):

    filenames = []
    if StrictVersion(version.strip('v')) < StrictVersion('0.2.23'):

        file_resources = list(tantalus_api.list(
            "file_resource",
            filename__startswith=jira_ticket,
            filename__endswith="{}_alignment_metrics.csv.gz".format(library_id)
        ))
        file_resources += list(tantalus_api.list(
            "file_resource",
            filename__startswith=jira_ticket,
            filename__endswith="{}_gc_metrics.csv.gz".format(library_id)
        ))
        file_resources += list(tantalus_api.list(
            "file_resource",
            filename__startswith=jira_ticket,
            filename__endswith="{}_multiplier0_metrics.csv.gz".format(library_id)
        ))
    else:
        file_resources = list(tantalus_api.list(
            "file_resource",
            filename__startswith=jira_ticket,
            filename__endswith="{}_alignment_metrics.csv.gz".format(library_id)
        ))
        file_resources += list(tantalus_api.list(
            "file_resource",
            filename__startswith=jira_ticket,
            filename__endswith="{}_gc_metrics.csv.gz".format(library_id)
        ))
        file_resources += list(tantalus_api.list(
            "file_resource",
            filename__startswith=jira_ticket,
            filename__endswith="{}_metrics.csv.gz".format(library_id)
        ))

    filenames += [f["filename"] for f in file_resources]

    return filenames

def get_cn_states(jira_ticket, version, library_id):

  filenames = []
  results_dataset = tantalus_api.get("resultsdataset",
                                     analysis__jira_ticket=jira_ticket,
                                     results_type="cell_state_prediction",
                                     results_version="v0.0.2"
                                     )

  f_resource = tantalus_api.get("file_resource", id=(results_dataset["file_resources"][0]))
  filenames += [f_resource["filename"]]
  return filenames

def get_all_files(jira_ticket, version, library_id):
    filenames = []
    filenames += get_segment(jira_ticket, version, library_id)
#     filenames += get_reads(jira_ticket, version, library_id)
#     filenames += get_metrics(jira_ticket, version, library_id)
#     if StrictVersion(version.strip('v')) < StrictVersion('0.3.1'):
#         filenames += get_cn_states(jira_ticket, version, library_id)

    return filenames

def main(library_id, download_dir):
    jira_ticket = find_latest_analysis_ticket(library_id) 
    # Hoa, add condition
    if jira_ticket:
        print('Latest Jira ticket  is: {0} correspond to library id: {1}'.format(jira_ticket,library_id))
        # Define source and destination storage clients
        from_storage_name = "singlecellresults"
        to_storage_name = download_dir

        from_storage = tantalus_api.get_storage(from_storage_name)
        from_storage_client = tantalus_api.get_storage_client(from_storage_name)

        # Check if user input is a jira ticket
        # If so, check if its a analysis ticket.
        # Else check if its a library ticket and return analysis tickets

        if re.match(JIRA_TICKET_REGEX, jira_ticket):

            try:
                analysis_object = colossus_api.get("analysis_information", analysis_jira_ticket=jira_ticket)
                analysis_ticket = jira_ticket
                # library_id = analysis_object["library"]["pool_id"]
                print("{} is a valid analysis ticket\n".format(analysis_ticket))

            except NotFoundError:
                raise Exception("{} is not an analysis ticket".format(jira_ticket))

        version = analysis_object["version"]
        filenames = get_all_files(analysis_ticket, version, library_id)
        downloaded_files = []
        for filename in filenames:
            filepath_parsed = filename.split("/")
            cell_state = filepath_parsed[1]
            analysis_type = filepath_parsed[-2]
            file = filepath_parsed[-1]
            print("analysis type: {0};  files: {1}".format(analysis_type,file))
            # REFACTOR
            # This takes care of results from old version of scpipeline
            if analysis_type == "plots":
                if "heatmap" in filename:
                    analysis_type = "hmmcopy_autoploidy"
                else:
                    analysis_type = "alignment"

            # accommadates weird cell_state_prediction naming /storage convention
            if cell_state == "cell_state_prediction":
                subdir = os.path.join(to_storage_name, library_id, cell_state)
                filepath = os.path.join(to_storage_name, library_id, cell_state, "{}_cell_state_prediction.csv".format(library_id))
            else:
                subdir = os.path.join(to_storage_name, library_id, analysis_type)
                filepath = os.path.join(to_storage_name, library_id, analysis_type, file)

            if not os.path.exists(subdir):
                os.makedirs(subdir)

            print("Downloading {} to {}".format(file, subdir))

            try:
                blob = from_storage_client.blob_service.get_blob_to_path(
                    container_name="results",
                    blob_name=filename,
                    file_path=filepath
                )
                downloaded_files.append(filepath)
            except:
                print('{} does not exist in blob.'.format(filename.split("/")[-1]))
                pass

        #clean up dir and unzip
        for root, dirs, files in os.walk(os.path.join(to_storage_name, library_id)):
            for name in files:
                if (os.path.getsize(os.path.join(root, name)) == 0 and re.search(".gz$", os.path.join(root, name))):
                    rm_cmd = ["rm", "-f", os.path.join(root, name)]
                    subprocess.check_call(rm_cmd)
                    continue

                elif re.search(".gz$", os.path.join(root, name)):
                    gunzip_cmd = ["gunzip", "-f", os.path.join(root, name)]
                    subprocess.check_call(gunzip_cmd)


        print("\n********** Download complete **********\n\n")
        for file in downloaded_files:
            print(file)
            
    else:
        print("Jira ticket do not exist for library id: {0}. Please double check the analysis".format(library_id))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--library_id', required=True)
    parser.add_argument('--download_dir', default='downloads')

    args = vars(parser.parse_args())
    print(args['download_dir'])
    print(args['library_id'])
    main(args['library_id'], args['download_dir'])
