#!/bin/bash
output_dir=/home/htran/storage/raw_DLP/metastasis_DLP/SA535/
input_jira_ticket="/home/htran/Projects/hakwoo_project/corrupt_tree/src/downloading_scripts/jira_tickets.txt"
results_blob_prefix="https://singlecellresults.blob.core.windows.net/results/"
results_blob_suffix_hmm="/results/hmmcopy/"
results_blob_suffix_anno="/results/annotation/"


library_ids=('A96201B' 'A98232A' 'A98277B' 'A98261A' 'A98197B' 'A98165A' 'A98194A' 'A98169A' 'A98174A' 'A95697A' 'A98163A' 'A98210A')
jira_tickets=('SC-3784' 'SC-3173' 'SC-2716' 'SC-3213' 'SC-3403' 'SC-3377' 'SC-2718' 'SC-3312' 'SC-3553' 'SC-2744' 'SC-3264' 'SC-3580')

for ((n=0;n<${#library_ids[@]};n++));
do
  echo "________________________________________________________"
  echo ${jira_tickets[n]}
  echo ${library_ids[n]} 
  
  base_dir="${output_dir}${library_ids[n]}"
  [ ! -d "$base_dir" ] && mkdir -p "$base_dir"
  # mkdir $base_dir

  echo "Downloading hmmcopy data"
  output_hmmcopy_fd="$base_dir/hmmcopy/"
  # mkdir $output_hmmcopy_fd
  [ ! -d "$output_hmmcopy_fd" ] && mkdir -p "$output_hmmcopy_fd"
  
  out_fn1="${output_hmmcopy_fd}${library_ids[n]}_reads.csv.gz"
  out_fn2="${output_hmmcopy_fd}reads.csv.gz"
  if [ -e $out_fn1 ] || [ -e $out_fn2 ]
  then
      echo "File exist, do nothing"
  else
      echo "Downloading file from azure"
      ## For old convention
      command_script1="az storage copy -s ${results_blob_prefix}${jira_tickets[n]}${results_blob_suffix_hmm}${library_ids[n]}_reads.csv.gz -d ${output_hmmcopy_fd}"
    
      ## For new convention, in case the old convention above do not work
      # command_script1="az storage copy -s ${results_blob_prefix}${jira_tickets[n]}${results_blob_suffix_hmm}reads.csv.gz -d ${output_hmmcopy_fd}"
      echo $command_script1
      eval "$command_script1"

  fi

  
  echo "Downloading hmmcopy annotation data"
  output_hmmcopy_anno="$base_dir/annotation/"
  # mkdir $output_hmmcopy_anno
  [ ! -d "$output_hmmcopy_anno" ] && mkdir -p "$output_hmmcopy_anno"

  out_fn11="${output_hmmcopy_anno}${library_ids[n]}_metrics.csv.gz"
  out_fn22="${output_hmmcopy_anno}metrics.csv.gz"
  if [ -e $out_fn11 ] || [ -e $out_fn22 ]
  then
      echo "ok, file exist, do nothing"
  else
      echo "Downloading file from azure"
      ## For old convention, file ends with ${library_ids[n]}_metrics.csv.gz
      command_script2="az storage copy -s ${results_blob_prefix}${jira_tickets[n]}${results_blob_suffix_anno}${library_ids[n]}_metrics.csv.gz -d ${output_hmmcopy_anno}"
    
      ## For new convention, in case the old convention above do not work, file ends with metrics.csv.gz
      # command_script2="az storage copy -s ${results_blob_prefix}${jira_tickets[n]}${results_blob_suffix_anno}metrics.csv.gz -d ${output_hmmcopy_anno}"
      echo $command_script2
      eval "$command_script2"
  fi

done




