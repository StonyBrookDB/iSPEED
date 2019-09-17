
#! /bin/bash

# parameters: inputpathdir outputdir prefix
# example /tmp /tmp compressed

echo "Running combiner"


if [ $# -ne 1 ]
  then
    echo "Not enough arguments. Need the prefix of the output directory "
    exit 1
fi


node_list="localhost"
data_dir=/ispeed_data

#default option
tmpdir=/tmp/tmpbin

inputbindir=/tmp
commonprefix=compressionoutputtask
outputfile=$data_dir/allbin
outputmbbs=$data_dir/allmbbs

rm -rf $tmpdir
mkdir -p $tmpdir 

echo "Copying data from nodes to the current node" 
# currently hard coding the list of nodes
for i in $node_list
do
	echo "copying data from node ${i}"
	scp ${i}:"${inputbindir}/${commonprefix}*" ${tmpdir}"/"
	# for now do not remove yet
	#ssh ${i} "rm ${inputbindir}/${commonprefix}*" 
done

# getting the total/global index
echo "combining the compression results"
echo "hdfs dfs -cat $1_mbb/0/* | ${HADOOPGIS_BIN_PATH}/compress_combine --inputbin  ${tmpdir} --outputbin ${outputfile} > ${outputmbbs}"
hdfs dfs -cat "$1""_mbb/0/*" | ${HADOOPGIS_BIN_PATH}/compress_combine --inputbin  ${tmpdir} --outputbin ${outputfile} > ${outputmbbs}
echo "done"

#remove the temporary files
#rm -rf ${tmpdir}

# Putting back the file with offset field to HDFS as input for Resque
hdfs dfs -rm -f "$1""_inputresque"
hdfs dfs -put -f ${outputmbbs} "$1""_inputresque"

# Remove all shared memory segments on all nodes
for i in $node_list
do 
	ssh $i "ipcrm -M 0x0000162e"
done

# Run loader
echo "Running loader to load compressed data on individual nodes"

for i in $node_list
do
	echo "loading on node $i"
	if [ "$i"=="localhost" ];
	then
		echo "${HADOOPGIS_BIN_PATH}/compress_load -n ${outputfile}"
		${HADOOPGIS_BIN_PATH}/compress_load -n ${outputfile}
	else
		scp -r ${outputfile} ${i}:${outputfile}
		ssh $i "${HADOOPGIS_BIN_PATH}/compress_load -n ${outputfile}"
	fi
	
	if [ $? -eq 0 ]; then
	    echo "Succeeded loading on $i"
	else
	    echo "Failed to load. Exiting"
	    exit 1
	fi
done

