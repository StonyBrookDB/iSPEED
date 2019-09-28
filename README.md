# iSPEED

You can try ispeed easily using docker. refer to https://www.docker.com/ for more information about install and use of docker container.

## build docker images

### build hadoop docker images

Firstly we need a Hadoop environment to store data (hdfs) and conduct computation (map reduce) across nodes. A standalone Hadoop environment can be run on a single machine with the help of docker-compose. 

Files used to build the docker-compose environment can be found in folder "docker-hadoop". The default Hadoop version is 3.1.2, which can be changed by modifying the docker-hadoop/base/Dockerfile. Some parameters of Hadoop can also be changed by modifying file docker-hadoop/hadoop.env . Since the binary executables need be running as mappers or reducers in the node manager, we build the Hadoop images with all the dependencies which can be found in docker-hadoop/base/Dockerfile. 

Then build the docker images for Hadoop with commands:

```
cd docker-hadoop
sh build-image.sh
```

You can now start the Hadoop docker containers with command:
```
docker-compose up
```
You'll see the servers of namenode, datanode, resourcemanager, nodemanager and historyserver are started in a certain order. Now test the hdfs with command:

```
docker run --rm --env-file=hadoop.env --net=host -it hadoop-base hadoop fs -ls /
```

### build iSPEED docker image

the dockerfile to build the iSPEED docker image is in the root folder of iSPEED project. It uses the image of the hadoop-base built in last step as the base image, since it contains the whole Hadoop command set and dependencies for iSPEED. Besides, it sets the environment variables for iSPEED and compile the binary executables needed by iSPEED. Build the image with command:

```
docker build -t ispeed .
```
Then test the 



