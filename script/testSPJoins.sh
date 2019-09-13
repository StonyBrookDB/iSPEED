#! /bin/bash

#../build/bin/queryproc --querytype spjoin --geom1 5 --geom2 5 --input1 /user/hoang/osm/rawdata/planet.tsv --input2 /user/hoang/osm/rawdata/europe.tsv -f 1:1,2:1 --outputpath /user/hoang/osmout2 --predicate st_intersects --partitioner fg --s 0.01 --numreducers 200
#samplingrate=0.01

# to test spatial join performance
bucketsize=800
#rm spjoin.log

for numreducers in 120 100 80 60 40 20
do 
    for files in data1 data2 data3 data4 data5 #data6
    do	
	START=$(date +%s)
	../build/bin/queryprocessor_3d -t st_intersects -a /user/yanhui/gis2016/spjoin/${files}/d1 -b /user/yanhui/gis2016/spjoin/${files}/d2 -h /user/yanhui/gis2016/spjoin/${files}/output/ -q spjoin -s 1.0 -n ${numreducers} -u fg_3d --bucket ${bucketsize} -f tileid,1:1,2:1 -i 1 -j 1  #-o overwrite
	succ=$?	
	END=$(date +%s)
	DIFF=$(($END-$START))
	#echo "${START},${END},${numreducers},${files},${DIFF}" >> spjoin.log
	if [[ $succ != 0 ]] ; then
    	echo "${START},${END},${numreducers},${files},${DIFF}, Job Failed" >> spjoin.log
    	#exit $succ
	else
	echo "${START},${END},${numreducers},${files},${DIFF}, Job Successed" >> spjoin.log	
	fi
	hdfs dfs -cat /user/yanhui/gis2016/spjoin/${files}/output_partidx/part* > /scratch/yliang/spjoin/part.${numreducers}.${files}.${bucketsize}.idx
   	hdfs dfs -rm -r /user/yanhui/gis2016/spjoin/${files}/output_joinout
   	hdfs dfs -rm -r /user/yanhui/gis2016/spjoin/${files}/output_mbb
   	hdfs dfs -rm -r /user/yanhui/gis2016/spjoin/${files}/output_partidx
    done
done

