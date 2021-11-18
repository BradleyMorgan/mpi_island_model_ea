#!/bin/bash

#  process_logs.sh
#  mpi_island_model_ea
#
#  Created by Bradley Morgan on 11/18/21.
#  Copyright © 2021 Bradley Morgan. All rights reserved.

DATE=`date "+%m%d%Y_%H%M"`
LOG_FN="logs_${DATE}.tar.gz"

REMOTE_PATH="~/mpi_island_model_ea/mpi_island_model_ea"
LOCAL_PATH="/Users/morgaia/Research/exaevo/results"
REMOTE_FN="${REMOTE_PATH}/${LOG_FN}"
LOCAL_FN="${LOCAL_PATH}/${FN}"

mkdir -p $LOCAL_PATH
ssh morgaia@hpclogin.auburn.edu tar -cvzf $REMOTE_FN $REMOTE_PATH/logs
scp morgaia@hpclogin.auburn.edu:$REMOTE_FN $LOCAL_FN


