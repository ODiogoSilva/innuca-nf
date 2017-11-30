#!/usr/bin/env python3

import sys
import json

f1, f2 = sys.argv[1:]

with open(f1) as f1h, open(f2) as f2h:

    j1 = json.load(f1h)
    j2 = json.load(f2h)

    l = {"cagao": [j1, j2]}

    with open(".report.json", "w") as fh:
        fh.write(json.dumps(l))
