#!/usr/bin/env bash

echo "{'project_id': $1, 'pipeline_id': $2, 'process_id': 1, " \
     "'status': $(pwd)/.status, 'log': $(pwd)/.command.log}, 'report': None, " \
     "'info': None, 'stats': None}"