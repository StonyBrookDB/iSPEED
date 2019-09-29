#!/bin/bash
sh clean.sh
queryprocessor_3d -q spjoin -a /user/teng/test1/input1 -b /user/teng/test1/input2 -h /user/teng/test1/output \
	 --binpath ./ --compressed_data_path /tmp/allbin -o --compression 
runcombiner.sh /user/teng/test1/output
queryprocessor_3d -q spjoin -a /user/teng/test1/input1 -b /user/teng/test1/input2 -h /user/teng/test1/output \
     --binpath ./ --compressed_data_path /tmp/allbin -u fg_3d -n 1 --spatialproc -t st_intersects
