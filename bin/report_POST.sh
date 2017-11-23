#!/usr/bin/env sh

projectid=$1
pipelineid=$2
processid=$3
sample=$4
reportjson=$5
url=$6

json="{'project_id':'$projectid','pipeline_id':'$pipelineid','process_id':'$processid', 'sample_name': '$sample', 'report_json': '$reportjson'}"

if [ -s .report.json ];
then
    curl -H  "Content-Type: application/json" -L -X POST -d \"$json\" $url > /dev/null
fi
