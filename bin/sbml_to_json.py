#!/usr/bin/env python3.7

import cobra
import sys

m=cobra.io.read_sbml_model(sys.argv[1])
cobra.io.save_json_model(m, sys.argv[2])
