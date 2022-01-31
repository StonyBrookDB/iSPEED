# iSPEED

## Further documentation:
http://bmidb.cs.stonybrook.edu/ispeed/index

You can try ispeed easily using docker. refer to https://www.docker.com/ for more information about install and use of docker container.

## Build docker images

### build hadoop docker images

Firstly we need a Hadoop environment to store data (hdfs) and conduct computation (map reduce) across nodes. A standalone Hadoop environment can be run on a single machine with the help of docker-compose. 

Files used to build the docker-compose environment can be found in folder "docker-hadoop". The default Hadoop version is 3.1.2, which can be changed by modifying the docker-hadoop/base/Dockerfile. Some parameters of Hadoop can also be changed by modifying file docker-hadoop/hadoop.env . Since the binary executables need be running as mappers or reducers in the node manager, we build the Hadoop images with all the dependencies which can be found in docker-hadoop/base/Dockerfile. 

Then build the docker images for Hadoop with commands:

```
cd docker-hadoop
sh build-image.sh
```

You can now create a "hadoop" network and start the Hadoop docker containers with command:
```
docker network create hadoop
docker-compose up
```
You'll see the servers of namenode, datanode are started in a certain order. Now test the hdfs with command:

```
docker run --rm --env-file=hadoop.env --net=hadoop -it hadoop-base hadoop fs -ls /
```
Note that here we just want to run iSPEED in standalone mode. For running iSPEED with cluster mode, use the docker-compose.full.yml to start the servers, which will also start the resource manager and node manager. Then also modify the framework name in the hadoop.env file from local to yarn. 


### build iSPEED docker image

the dockerfile to build the iSPEED docker image is in the root folder of iSPEED project. It uses the image of the hadoop-base built in last step as the base image, since it contains the whole Hadoop command set and dependencies for iSPEED. Besides, it sets the environment variables for iSPEED and compile the binary executables needed by iSPEED. Build the image with command:

```
docker build -t ispeed .
```
## Upload data

### generate data
we put several .off files in the data folder for testing the functionality of iSPEED. You can use our tool to generate any amount of data for scalability testing (will be added shortly).
### upload data to hdfs
After generating proper amount of testing data, upload those data files into hdfs with commands:
```
hadoop fs -put /path/to/off/file /path/on/hdfs
```

## Run tests
The scripts in the /script folder can be used to run tests by your own. clean.sh delete the output folders and temporary files, you may want to modify it to fulfill your needs. runcombiner.sh is called by runtest.sh and runtest_nn.sh to combine the compressed data which is the output of the preprocessing step. After uploaded data into folders in hdfs, run command:
```
docker run --rm --env-file=hadoop.env --net=hadoop -it ispeed bash
```
It will run a shell which you can run tests. Then you can run tests with command for testing:
```
sh runtest.sh
```

## References:
[1] Yanhui Liang, Hoang Vo, Jun Kong and Fusheng Wang: [iSPEED: an Efficient In-Memory Based Spatial Query System for Large-Scale 3D Data with Complex Structures](https://dl.acm.org/doi/10.1145/3139958.3139961). In Proceedings of SIGSPATIAL 2017. November 7 - 10, 2017, Redondo Beach, California, USA.

[2] Yanhui Liang, Hoang Vo, Jun Kong and Fusheng Wang: [iSPEED: a Scalable and Distributed In-Memory Based Spatial Query System for Large and Structurally Complex 3D Data](http://www.vldb.org/pvldb/vol11/p2078-vo.pdf). Demo Paper. In Proceedings of the 44th International Conference on Very Large Data Bases (VLDB 2018), August 27-31, 2018, Rio de Janeiro, Brazil.
