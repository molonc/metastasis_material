#!/bin/bash

input_dir=/home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/
save_dir=/home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/test_TreeAlign/
script_exe=/home/htran/Projects/farhia_project/rnaseq/method_testing/TreeAlign/run_TreeAlign/exe_TreeAlign.py
python_env=/home/htran/anaconda3/envs/sisyphus/bin/python
meta_sample_ids=/home/htran/storage/rnaseq_datasets/testing_space/clonealign2/SA535X4XB05462/samples_run_list.txt

IFS=$'\n' read -d '' -r -a datatags < $meta_sample_ids


nb_iterations=500
min_clone_assign_prob=0.7


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