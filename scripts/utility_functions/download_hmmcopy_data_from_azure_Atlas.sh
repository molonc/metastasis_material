#!/bin/bash
output_dir=/home/htran/storage/raw_DLP/metastasis_DLP/SA919/hmmcopy/
results_blob_prefix="https://bccrcprccatlassa.blob.core.windows.net/atlas/"
results_blob_suffix="/results/"

echo "Get list of AtlasId from an input file and downloading hmmcopy data"
input_lib_ticket=/home/htran/storage/raw_DLP/metastasis_DLP/SA919/hmmcopy/atlasIds.txt

IFS=$'\n' read -d '' -r -a lib_tickets < $input_lib_ticket



for ((n=0;n<${#lib_tickets[@]};n++));
do
  echo "________________________________________________________"
  printf "Atlas id 1: %s\n" "${lib_tickets[n]}"
  #echo "${lib_tickets[n]}" | cut -d "," -f 1
  #echo "${lib_tickets[n]}" | cut -d "," -f 2
  ax_id=$(echo "${lib_tickets[n]}" | cut -d "," -f 1 )
  library_id=$(echo "${lib_tickets[n]}" | cut -d "," -f 2 )
  #AT_id=$(echo "${lib_tickets[n]}" | cut -d "," -f 3 )
  echo "$ax_id"
  echo "$library_id"
  #echo "$AT_id"
  
  output_fd="${output_dir}${library_id}/"
  [ ! -d "$output_fd" ] && mkdir -p "$output_fd"
  
  if ls "$output_fd"/*_reads.csv.gz 1>/dev/null 2>&1; then
      echo "File exist, do nothing"
  else
      echo "Downloading file from Atlas azure"
      ## For old convention
      command_script1="az storage copy -s ${results_blob_prefix}${ax_id}${results_blob_suffix} -d ${output_fd} --include-pattern '*hmmcopy_reads.csv.gz;*hmmcopy_metrics.csv.gz' --recursive"
      echo $command_script1
      eval "$command_script1"

  fi

  
done  


## az storage copy -s https://bccrcprccatlassa.blob.core.windows.net/atlas/AX2024/results/F132798/abundance.tsv -d .
