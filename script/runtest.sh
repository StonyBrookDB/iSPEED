#!/bin/bash
sh clean.sh
queryprocessor_3d -q spjoin -a /user/teng/test1/input1 -b /user/teng/test1/input2 -h /user/teng/test1/output -u fg_3d -o --compression --binpath ~/project/iSPEED/build/bin/
sh ~/project/iSPEED/script/runcombiner.sh /user/teng/test1/output
queryprocessor_3d -q spjoin -a /user/teng/test1/input1 -b /user/teng/test1/input2 -h /user/teng/test1/output -u fg_3d -n 1 --spatialproc --binpath ~/project/iSPEED/build/bin/
