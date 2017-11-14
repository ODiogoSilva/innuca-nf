#!/usr/bin/env bash

json="{'project_id':'$1','pipeline_id':'$2','process_id':'$3','run_property':'log_file','run_property_value':'$(pwd)/.command.log','type':'output'}"

curl -H  "Content-Type: application/json" -L -X PUT -d \"$json\" $4
