 #!/usr/bin/env python

from dbclients.colossus import ColossusApi
import argparse

parser = argparse.ArgumentParser(description="""Description: get the latest compeleted analysis jira from a library""")
parser.add_argument("-l", nargs='+' ,help="Library ids", required=True)
args=parser.parse_args()
colossus_api = ColossusApi()

for l in args.l:
    biggest_id=0

    obj = list(colossus_api.list('analysis_information',library__pool_id=l))

    for info in obj:
        if int(info["id"]) > biggest_id and info["montage_status"]=="Success":
            biggest_id = info["id"]

    for info in obj:
        if info["id"] == biggest_id:
            print(f"{l},",info["analysis_jira_ticket"])

