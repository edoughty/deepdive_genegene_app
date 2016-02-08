#! /usr/bin/env bash

export APP_HOME=`pwd`
export DEEPDIVE_HOME=`cd ${APP_HOME}/../../../deepdive; pwd`

# Database Configuration
export PGUSER=${PGUSER:-`whoami`}
export PGPASSWORD=${PGPASSWORD:-}
export PGDATABASE=${PGDATABASE}
export PGHOST=${PGHOST}
export PGPORT=${PGPORT}
export PARALLELISM=48
export MEMORY="32g"
export JAVA_OPTS="-Xmx"$MEMORY
export SBT_OPTS="-Xmx"$MEMORY

export LD_LIBRARY_PATH=[DEEPDIVE_HOME]/util/dw_linux/lib/protobuf/lib:[DEEPDIVE_HOME]/util/dw_linux/lib64
export PYTHONPATH=$DEEPDIVE_HOME/ddlib/:$PYTHONPATH

cd $DEEPDIVE_HOME

$DEEPDIVE_HOME/sbt/sbt "run -c $APP_HOME/application.conf" $@
