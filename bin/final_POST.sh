#!/usr/bin/env sh

st=$(cat $(pwd)/.status)
warn=$(cat $(pwd)/.warnings)

json="{'project_id':'$1','pipeline_id':'$2','process_id':'$3','run_info':'None','run_output':'None','warnings':'$warn','log_file':'$(pwd)/.command.log','status':'$st','type':'output'}"

curl -H  "Content-Type: application/json" -L -X POST -d \"$json\" $4
