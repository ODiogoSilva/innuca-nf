#!/usr/bin/env sh

projectid=$1
pipelineid=$2
processid=$3
sample=$4
url=$5
username=$6
userid=$7


if [ -s .report.json ];
then
    json_str=$(cat $(pwd)/.report.json)
    echo json_str
    json="{'project_id':'$projectid','pipeline_id':'$pipelineid','process_id':'$processid','sample_name':'$sample','report_json':'$json_str','current_user_name':'$username','current_user_id':'$userid'}"
    curl -H  "Content-Type: application/json" -k -L -X POST -d \"$json\" $url > /dev/null
fi
