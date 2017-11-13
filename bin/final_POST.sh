#!/usr/bin/env bash

json="{'project_id': '$1', 'pipeline_id': '$2', 'process_id': '$3', 'run_fino': 'None', 'run_output': 'None, 'run_stats': 'None', 'log_file': '$(pwd)/.command.log', 'status': '$(pwd)/.status'}"

curl -H  "Content-Type: application/json" \
     -X POST -d $json \
     $4
