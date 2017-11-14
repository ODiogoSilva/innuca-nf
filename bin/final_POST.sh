#!/usr/bin/env sh

st=`cat $(pwd)/.status`

json="{'project_id':'$1','pipeline_id':'$2','process_id':'$3','run_info':'None','run_output':'None','run_stats':'None','log_file':'$(pwd)/.command.log','status':'$st','type':'output'}"

curl -H  "Content-Type: application/json" -L -X POST -d \"$json\" $4
