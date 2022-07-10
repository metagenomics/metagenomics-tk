#!/usr/bin/env python3

import cobra
import sys

m=cobra.io.read_sbml_model(sys.argv[1], cache=False)
cobra.io.save_json_model(m, sys.argv[2])
