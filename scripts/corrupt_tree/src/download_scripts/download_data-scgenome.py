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

import scgenome.db.qc
# from dbclients.tantalus import TantalusApi
from dbclients.colossus import ColossusApi
from dbclients.basicclient import NotFoundError

# tantalus_api = TantalusApi()
colossus_api = ColossusApi()



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




def main(library_id, download_dir):
    jira_ticket = find_latest_analysis_ticket(library_id) 
    # Hoa, add condition
    if jira_ticket:
        print('Latest Jira ticket  is: {0} correspond to library id: {1}'.format(jira_ticket,library_id))

        if not os.path.exists(download_dir): os.makedirs(download_dir)
            
        lib_dir = os.path.join(download_dir,library_id)
        ticket_dir = os.path.join(download_dir,jira_ticket) #,'results'
            
        # Define source and destination storage clients
#         from_storage_name = "singlecellresults"
#         to_storage_name = download_dir
        JIRA_TICKET_REGEX = r".*SC-\d{4}$"
        if re.match(JIRA_TICKET_REGEX, jira_ticket):
            results_tables = scgenome.db.qc.get_qc_data(
                [jira_ticket],
                download_dir,
                do_caching=True
            )
            if not os.path.exists(lib_dir): os.makedirs(lib_dir)
            os.system("mv "+ticket_dir+"/results/*"+" "+lib_dir)
            os.system("rm -fr "+ticket_dir)
            
            
            
    else:
        print("Jira ticket do not exist for library id: {0}. Please double check the analysis".format(library_id))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--library_id', required=True)
    parser.add_argument('--download_dir', default='.')

    args = vars(parser.parse_args())
    print(args['download_dir'])
    print(args['library_id'])
    main(args['library_id'], args['download_dir'])

    
    # download_dir = '/home/htran/storage/datasets/drug_resistance_DLP/SA535'
    # # Download data SA535
    # libs =['A62397A', 'A95625B', 'A98163A', 'A98210A']
    # libs = ['A96233A','A98285A','A108879A','A98168B','A95625A',
    # 'A98259B','A98287A','A62397A','A95625B','A98253B','A98256B',
    # 'A108833B','A108851B','A108870A','A108874A','A108757B','A108768B',
    # 'A108871A','A108874B','A108837A','A108863A','A98168A','A98244A',
    # 'A108732B','A108759B']
    # 
    # for l in libs:
    #   main(l, download_dir)
    
