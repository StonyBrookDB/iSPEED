#!/bin/bash
echo "waiting for 30 seconds to ensure the starts of namenode"
sleep 30
$HADOOP_HOME/bin/yarn --config $HADOOP_CONF_DIR resourcemanager
