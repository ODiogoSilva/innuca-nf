#!/usr/bin/env sh

set -ex

projectid=$1
pipelineid=$2
processid=$3
sample=$4
url=$5
username=$6
userid=$7
task=$8
species=$9


if [ -s .report.json ];
then

    if [ $task = "chewbbaca" ];
    then
        json_str=$(cat $(pwd)/.report.json | sed 's/ //g' | sed s/\"/\'/g)
    else
        json_str=$(cat $(pwd)/.report.json | sed 's/ /_/g' | sed s/\"/\'/g)
    fi

    workdir=$(pwd)
    json="{'project_id':'$projectid','pipeline_id':'$pipelineid','process_id':'$processid','sample_name':'$sample','report_json':$json_str,'current_user_name':'$username','current_user_id':'$userid','workdir':'$workdir','task':'$task','species':'$species'}"
    echo \"${json}\" > .final.json
    cat .final.json | curl -H  "Content-Type: application/json" -k -L -X POST -d @- $url > /dev/null
fi
