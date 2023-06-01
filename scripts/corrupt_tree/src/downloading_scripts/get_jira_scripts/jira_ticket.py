
# Hoa Tran
import os
import re
import sys
import argparse
import logging
import pandas as pd

#import subprocess
#from tabulate import tabulate
#from distutils.version import StrictVersion

# import scgenome.db.qc
# from dbclients.tantalus import TantalusApi
#from dbclients.basicclient import NotFoundError
from dbclients.colossus import ColossusApi



# tantalus_api = TantalusApi()
colossus_api = ColossusApi()


def find_latest_analysis_ticket_v2(lib_ids, output_dir):
    jira_tks=[]
    for l in lib_ids:
        tk = ""
        biggest_id=0
    
        obj = list(colossus_api.list('analysis_information',library__pool_id=l))
    
        for info in obj:
            if int(info["id"]) > biggest_id and info["montage_status"]=="Success":
                biggest_id = info["id"]
    
        for info in obj:
            if info["id"] == biggest_id:
                print(f"{l},",info["analysis_jira_ticket"])
                tk = info["analysis_jira_ticket"]
        
        
        jira_tks.append(tk)        
    
    
    df = pd.DataFrame({'library_id':lib_ids,'jira_ticket':jira_tks})
    df.to_csv(os.path.join(output_dir, 'jira_tickets.csv'), index=False)
    return df



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
        print(analysis["analysis_jira_ticket"])
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




def main(library_id):
    
    jira_ticket = find_latest_analysis_ticket(library_id)
    
    if jira_ticket:
        print('Latest Jira ticket  is: {0} correspond to library id: {1}'.format(jira_ticket,library_id))

        #print(jira_ticket)
        JIRA_TICKET_REGEX = r".*SC-\d{4}$"
        if re.match(JIRA_TICKET_REGEX, jira_ticket):
            print('Voila your ticket: {0}'.format(jira_ticket))  
            print('Corresponding to DLP library id: {0}'.format(library_id))
    
    else:
        print("Jira ticket do not exist for library id: {0}. Please double check the analysis".format(library_id))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Description: get the latest compeleted analysis jira from a library""")
    parser.add_argument('--library_id', nargs='+' ,help="Library ids", required=True)
    # parser.add_argument("-l", nargs='+' ,help="Library ids", required=True)
    parser.add_argument('--output_dir', required=True)
    
    args = vars(parser.parse_args())
    # print(args['library_id'])
    # main(args['library_id'])
    # lib_ids = ["A98165A","A98197B"]
    # output_dir = '/home/htran/storage/raw_DLP/metastasis_DLP/SA535/testing/'
    lib_ids = args['library_id']
    output_dir = args['output_dir']
    print(lib_ids)
    find_latest_analysis_ticket_v2(lib_ids, output_dir)
    
    
### How to run script from commandline mode: 
### /home/htran/anaconda3/envs/sisyphus38/bin/python jira_ticket.py --library_id A98165A A98197B --output_dir /home/htran/storage/raw_DLP/metastasis_DLP/SA535/testing/
