
How to run this function: 

Rscript yourdir/new_encoding_sitka/blob_cn_to_binary.R -i yourdir/new_encoding_sitka/testing_data_sitka/A98207B_filtered_states.csv -o  yourdir/new_encoding_sitka/testing_data_sitka/binary_encoded.csv -f 0

## Note: improve the encoding output, carefully selecting genomic bins
## We're missing synch events across the genome, i.e., those that happen in multiple cells, in more than one bin. Here we add the padding at the left
## to take into account the synch events.
More details are noted in the blob_cn_to_binary.R script. 





