import sys
import os
# import re
# import gzip
import click
# import subprocess
import pandas as pd
import numpy as np
# from natsort import natsorted
# from dbclients.colossus import ColossusApi

# colossus_api = ColossusApi()

@click.command()
@click.argument("files", nargs = -1)
@click.argument("output_file", nargs = 1)
def merge_libs(**kwargs):
    # Do nothing at this moment
    
    log_df = pd.DataFrame({'lib_id': kwargs["files"]})
    log_df.to_csv(kwargs["output_file"], sep = ",", index_label = False)
        # pb = 0.8
        # get_noisy_cells(total_state_samples, kwargs["output_file"], pb)
        
        

if __name__=='__main__':
    merge_libs()
