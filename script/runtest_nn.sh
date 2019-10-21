#!/bin/bash
root_folder="/user/teng/test1"
input_1="input1"
input_2="input_vessel_1"
output_folder="output"
predicate_type="st_nn_rtree"

clean.sh

#firstly: compress the data
queryprocessor_3d -q compress -o $root_folder/$output_folder --binpath ./ \
	-a $root_folder/$input_1 -b $root_folder/$input_2 

#secondly: combine the compressed data and mbbs
runcombiner.sh "$root_folder/$output_folder"

#thirdly: partition the space 
queryprocessor_3d -q partition -o $root_folder/$output_folder --binpath ./ \
	--compressed_data_path /tmp/allbin -t fg_3d

#forthly: do the join
queryprocessor_3d -q join -o $root_folder/$output_folder --binpath ./ \
	--compressed_data_path /tmp/allbin -p $predicate_type
