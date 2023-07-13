#!/bin/bash

input_dir=/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_introns/input_data/
save_dir=/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_introns/TreeAlign_result_08/
script_exe=/home/htran/Projects/hakwoo_project/metastasis_material/scripts/method_testing/run_TreeAlign/exe_TreeAlign.py
python_env=/home/htran/anaconda3/envs/sisyphus/bin/python
meta_sample_ids=/home/htran/storage/datasets/metastasis_results/rnaseq_SA535/clonealign_introns/samples_run_list.txt

## Read list of samples into an array
IFS=$'\n' read -d '' -r -a datatags < $meta_sample_ids


nb_iterations=500
min_clone_assign_prob=0.8


for ((n=0;n<${#datatags[@]};n++));
do
  echo "________________________________________________________"
  printf "Testing sample: %s\n" "${datatags[n]}"
  
  
  base_dir="${input_dir}"
  [ ! -d "$base_dir" ] && mkdir -p "$base_dir"
  # mkdir $base_dir
  
  expr_fn="${input_dir}${datatags[n]}_expr_introns.csv.gz"
  
  cnv_fn="${input_dir}${datatags[n]}_clones_cnv.csv.gz"
  
  cell_clones_fn="${input_dir}${datatags[n]}_cell_clones.csv.gz"
  
  [ ! -d "$save_dir" ] && mkdir -p "$save_dir"
  
  echo "Running tree align"
  
  command_script="${python_env} ${script_exe} --expr_fn ${expr_fn} --cnv_fn ${cnv_fn} --cell_clones_fn ${cell_clones_fn} --datatag ${datatags[n]} --save_dir ${save_dir} --iterations ${nb_iterations} --assign_prob ${min_clone_assign_prob}"

 
  echo $command_script
  eval "$command_script"

  
  
  
done  
