echo "Downloading libraries in input.txt"
cat input.txt | xargs -I % %
cat input.txt | xargs -I % ~/azcopy copy https://singlecellresults.blob.core.windows.net/results/%/results/hmmcopy/ ./%/ --recursive 
