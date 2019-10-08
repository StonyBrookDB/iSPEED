#!/bin/bash
root_folder="/user/teng/test1"
input_1="input1"
input_2="input_vessel"
output_folder="output"
query_type="st_nn_rtree"

clean.sh
queryprocessor_3d -q spjoin -a $root_folder/$input_1 -b $root_folder/$input_2 -h $root_folder/$output_folder \
	 --binpath ./ --compressed_data_path /tmp/allbin -o --compression
runcombiner.sh "$root_folder/$output_folder"
queryprocessor_3d -q spjoin -a $root_folder/$input_1 -b $root_folder/$input_2 -h $root_folder/$output_folder \
	 --binpath ./ --compressed_data_path /tmp/allbin -u fg_3d -n 1 --spatialproc -t $query_type
