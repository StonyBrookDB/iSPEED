#!/bin/bash
sh clean.sh
queryprocessor_3d -q spjoin -a /user/teng/test1/input_tiny1 -b /user/teng/test1/input_vessel_tiny1 -h /user/teng/test1/output \
	 --binpath ~/project/iSPEED/build/bin/ --compressed_data_path /ispeed_data/allbin -o --compression
runcombiner.sh /user/teng/test1/output
queryprocessor_3d -q spjoin -a /user/teng/test1/input_tiny1 -b /user/teng/test1/input_vessel_tiny1 -h /user/teng/test1/output \
	 --binpath ~/project/iSPEED/build/bin/ --compressed_data_path /ispeed_data/allbin -u fg_3d -n 1 --spatialproc -t st_nn_voronoi
