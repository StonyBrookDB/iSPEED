#!/bin/bash
docker rmi $(docker images --format '{{.Repository}}:{{.Tag}}' | grep 'hadoop-')
