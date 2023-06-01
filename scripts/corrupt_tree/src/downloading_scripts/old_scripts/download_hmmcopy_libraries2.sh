#!/bin/bash
output_dir=/home/htran/storage/raw_DLP/metastasis_DLP/SA535/
input_jira_ticket="/home/htran/Projects/hakwoo_project/corrupt_tree/src/downloading_scripts/jira_tickets.txt"
results_blob_prefix="https://singlecellresults.blob.core.windows.net/results/"
results_blob_suffix_hmm="/results/hmmcopy/"
results_blob_suffix_anno="/results/annotation/"
# results_blob_suffix_anno2="metrics.csv.gz"

while IFS= read -r line
do
  if [ -z "$line" ]
  then
        # echo "\$line is empty"
        echo ""
  else
        echo "-------------------------------------------------"
        echo $line
        jira_ticket=${line%,*}
        library_id=${line#*,}
        echo $jira_ticket
        echo $library_id
        base_dir="${output_dir}${library_id}"
        [ ! -d "$base_dir" ] && mkdir -p "$base_dir"
        # mkdir $base_dir
        
        echo "Downloading hmmcopy data"
        output_hmmcopy_fd="$base_dir/hmmcopy/"
        # mkdir $output_hmmcopy_fd
        [ ! -d "$output_hmmcopy_fd" ] && mkdir -p "$output_hmmcopy_fd"
      
        ## For old convention
        command_script1="az storage copy -s ${results_blob_prefix}${jira_ticket}${results_blob_suffix_hmm}${library_id}_reads.csv.gz -d ${output_hmmcopy_fd}"
        
        ## For new convention, in case the old convention above do not work
        # command_script1="az storage copy -s ${results_blob_prefix}${jira_ticket}${results_blob_suffix_hmm}reads.csv.gz -d ${output_hmmcopy_fd}"
        echo $command_script1
        # az storage copy -s ${results_blob_prefix}${jira_ticket}${results_blob_suffix_hmm}${library_id}_reads.csv.gz -d ${output_hmmcopy_fd}
        # eval "$command_script1"
        # az storage copy -s $results_blob_prefix$jira_ticket$results_blob_suffix_hmm -d $output_hmmcopy_fd --recursive
        
        echo "Downloading hmmcopy annotation data"
        output_hmmcopy_anno="$base_dir/annotation/"
        # mkdir $output_hmmcopy_anno
        [ ! -d "$output_hmmcopy_anno" ] && mkdir -p "$output_hmmcopy_anno"
        
        ## For old convention, file ends with ${library_id}_metrics.csv.gz
        command_script2="az storage copy -s ${results_blob_prefix}${jira_ticket}${results_blob_suffix_anno}${library_id}_metrics.csv.gz -d ${output_hmmcopy_anno}"
        
        ## For new convention, in case the old convention above do not work, file ends with metrics.csv.gz
        # command_script2="az storage copy -s ${results_blob_prefix}${jira_ticket}${results_blob_suffix_anno}metrics.csv.gz -d ${output_hmmcopy_anno}"
        echo $command_script2
        # az storage copy -s ${results_blob_prefix}${jira_ticket}${results_blob_suffix_anno}${library_id}_metrics.csv.gz -d ${output_hmmcopy_anno}
        # eval "$command_script2"
        # az storage copy -s "${results_blob_prefix}${jira_ticket}${results_blob_suffix_anno}" -d $output_hmmcopy_anno --recursive
        # az storage copy -s $results_blob_prefix$jira_ticket$results_blob_suffix_anno -d $output_hmmcopy_anno --recursive
        
  fi    

done < "$input_jira_ticket"




