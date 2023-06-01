import sys
import os
import re
import gzip
import argparse
import subprocess
import pandas as pd
import numpy as np
from natsort import natsorted
from distutils.version import StrictVersion
from dbclients.colossus import ColossusApi

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

def filter_libs(**kwargs):
    '''
    (1) Given sample_ids and cached directory, filter raw reads and metrics for mappability >=0.99;
                                                                                read quality >=0.75;
                                                                                total mapped reads >= 500000;
                                                                                live cells;
                                                                                non S-phase cells;
                                                                                contaminants of SA501 & SA928;

    (2) Reshape to extract state, copy, and quantile transformed copy data

    Args:
        sample ids ("," sep str)
        cached dir path (path)
        good_cell_dir (path)
        output_dir (path)

    Returns:
        total combined state df
    '''

    best_loci = set()
    total_state_samples = pd.DataFrame()
    library_id = kwargs['library_id']
    jira_id = find_latest_analysis_ticket(library_id)
    analyses = list(colossus_api.list("analysis_information", analysis_jira_ticket=jira_id))

    if len(analyses) == 1:
        print("Only 1 analysis info found. Now processing", jira_id)
        analysis = analyses[0]
        lib_id = analysis['library']['pool_id']

        if 'samples' in kwargs:
            sample_id = list(kwargs['samples'].split(','))
        else:
            sample_id = [analysis['library']['sample']['sample_id']]
        version = analysis['version']

        filtered_metrics = pd.DataFrame()
        filtered_reads = pd.DataFrame()
        tmp_metrics = pd.DataFrame()

        experimental_condition = list(colossus_api.list("experimental_metadata", library__pool_id=lib_id))

        filter_cond = ['A']
        sample_type = 'cell'
        for exp_cond in experimental_condition:
            for chip in exp_cond['chipregionmetadata_set']:
                if chip['metadata_value'] in sample_id and not re.search("NCC|NTC|gDNA", exp_cond['region_code']):
                    filter_cond.append(exp_cond['region_code'])
                if exp_cond['region_code'] in filter_cond and chip['metadata_field'] == "sample_type":
                    sample_type = chip['metadata_value']

        if sample_type == "cell" or sample_type == "cells":
            cell_call_filter = "C1"
        if sample_type == "nuclei":
            cell_call_filter = "C2"

        # print(sample_type, filter_cond)

        # print(cell_call_filter)

        metrics_df = pd.read_csv(kwargs['metrics'])
        filtered_metrics = metrics_df.loc[metrics_df['quality'] >= 0.75]
        # print('quality', filtered_metrics.shape)
        if 'sample_id' not in filtered_metrics.columns:
            # print(filtered_metrics['cell_id'])
            filtered_metrics = filtered_metrics.assign(sample_id = filtered_metrics.apply(lambda row: row['cell_id'].split("-")[0], axis = 1))
            # print(filtered_metrics['sample_id'])
        for sample in sample_id:
            # print(sample)
            tmp_metrics = tmp_metrics.append(filtered_metrics[(filtered_metrics['sample_id'].str.startswith(sample))])
            # print(tmp_metrics.shape)
        filtered_metrics = tmp_metrics
        # print("sample filter", filtered_metrics.shape)
        filtered_metrics = filtered_metrics[(filtered_metrics['cell_call'] == cell_call_filter)]
        # print('cell call filter', filtered_metrics.shape)
        filtered_metrics = filtered_metrics[(filtered_metrics["experimental_condition"].isin(filter_cond))]
        # print('experimental_condition', filtered_metrics.shape)
        filtered_metrics = filtered_metrics['cell_id']

        reads_df = pd.read_csv(kwargs['reads'])
        filtered_reads = reads_df.loc[reads_df['map'] >= 0.99]
        # print('mapping quality', filtered_reads.shape)

        if StrictVersion(version.strip('v')) < StrictVersion('0.3.1'):
            cell_predict_df = pd.read_csv(kwargs['cell_state_predict'])
        else:
            cell_predict_df = pd.read_csv(kwargs['metrics'])

        cell_list_filtered = []
        for sample in sample_id:
            cell_list_filtered += list(filter(lambda x: re.search(sample, x), cell_predict_df['cell_id']))

        cell_predict_df = cell_predict_df[cell_predict_df['cell_id'].isin(cell_list_filtered)]
        cell_predict_df = cell_predict_df[(cell_predict_df['is_s_phase'] == False)]
        not_s_phase = cell_predict_df['cell_id']

        if not (filtered_reads.empty and filtered_metrics.empty and not_s_phase.empty):
            filtered_cnv_data = filtered_reads[filtered_reads['cell_id'].isin(filtered_metrics)]
            filtered_cnv_data = filtered_cnv_data[filtered_cnv_data['cell_id'].isin(not_s_phase)]

            if filtered_cnv_data['state'].dtype == np.float64:
                filtered_cnv_data['state'] = filtered_cnv_data['state'].astype(int)

            filtered_cnv_data['loci'] = filtered_cnv_data[['chr', 'start', 'end']].apply(lambda x: '_'.join(map(str, x)), axis=1)
            filtered_state_data = filtered_cnv_data.pivot(index='loci', columns='cell_id', values='state')
            filtered_state_data = filtered_state_data.dropna(axis='columns', how='all')

        filtered_state_data = filtered_state_data.reindex(index=natsorted(filtered_state_data.index))
        print("dimensions of filtered data", filtered_state_data.shape)
        filtered_state_data.to_csv(kwargs['output_file'], sep=',', index_label=False)

    else:
        print("more than one analysis information associated with " + jira_id + " ticket. Please confirm on Colossus")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--metrics', required=True)
    parser.add_argument('--reads', required=True)
    parser.add_argument('--cell_state_predict', required=False, default=None)
    parser.add_argument('--library_id', required=True)
    parser.add_argument('--samples', required=True)
    parser.add_argument('--output_file', required=True)

    args = vars(parser.parse_args())
    
    library_id = args['library_id']
    samples = args['samples']
    metrics = args['metrics']
    reads = args['reads']
    cell_state_predict = args['cell_state_predict']
    output_file = args['output_file']
    filter_libs(library_id=library_id, samples=samples, metrics=metrics, reads=reads, cell_state_predict=cell_state_predict, output_file=output_file)
