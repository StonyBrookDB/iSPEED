#! /bin/bash

echo "building base"
cd base
docker build -t hadoop-base:latest ./
echo "building datanode"
cd ../datanode
docker build -t hadoop-datanode:latest ./
echo "building namenode"
cd ../namenode
docker build -t hadoop-namenode:latest ./
echo "building resourcemanager"
cd ../resourcemanager
docker build -t hadoop-resourcemanager:latest ./
echo "building nodemanager"
cd ../nodemanager
docker build -t hadoop-nodemanager:latest ./
echo "building historyserver"
cd ../historyserver
docker build -t hadoop-historyserver:latest ./