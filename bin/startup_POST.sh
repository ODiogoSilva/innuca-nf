#!/usr/bin/env bash

#curl -H  "Content-Type: application/json" \
#     -X POST -d $1 \
#     https://

echo "{'project_id': '$1', 'pipeline_id': '$2', 'process_id': '1', " \
     "'status': 'None', 'log': '$(pwd)/.command.log', 'report': 'None', " \
     "'info': 'None', 'stats': 'None'}"