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

json_str=""
version_str=""

# If a .report.json file was populated, set the json_str variable
if [ -s .report.json ];
then

    # Modification of the JSON string should be different for chewbbaca
    # output
    if [ $task = "chewbbaca" ];
    then
        json_str=$(cat $(pwd)/.report.json | sed 's/ //g' | sed s/\"/\'/g)
    else
        json_str=$(cat $(pwd)/.report.json | sed 's/ /_/g' | sed s/\"/\'/g)
    fi
fi

# If a .versions file was populated, set the version_str variable
if [ -s .versions ];
then
    version_str=$(cat $(pwd).versions | sed 's/ /_/g' | sed s/\"/\'/g)
fi

# If a .versions OR .report.json file was populated send the request
if [ ! "$json_str" = "" ] || [ ! "$version_str" = "" ];
then
    workdir=$(pwd)
    json="{'project_id':'$projectid','pipeline_id':'$pipelineid','process_id':'$processid','sample_name':'$sample','report_json':$json_str,'current_user_name':'$username','current_user_id':'$userid','workdir':'$workdir','task':'$task','species':'$species', 'versions':'$version_str'}"
    echo \"${json}\" > .final.json
    cat .final.json | curl -H  "Content-Type: application/json" -k -L -X POST -d @- $url > /dev/null
fi
