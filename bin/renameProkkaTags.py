#!/usr/bin/env python
# This script reads a prokka_tags_tmp.tsv file and generates a dictionary with old tags as keys and new tags as values.
# It finds all files with the specified extensions and applies the substitutions to each file separately.
#
# Usage: python renameProkkaTags.py <BIN_PREFIX>

from __future__ import print_function
from multiprocessing import Pool
import os
import re
import sys
import csv

# Read the BIN_PREFIX value from the command line arguments
BIN_PREFIX = sys.argv[1]
# Number of CPUs to use for parallelization, if not specified use 2 CPUs
CPUS = int(sys.argv[2]) if len(sys.argv) > 2 else 2
locus_tag_prefix = ""

# Read the prokka_tags_tmp.tsv file and generate a dictionary with old tags as keys and new tags as values
# Example of the prokka_tags_tmp.tsv file:
#
# SAMPLE  CONTIG  locus_tag
# test1   test1_19_0ea371 BLGMGGNG_00001
#
with open("prokka_tags_tmp.tsv", "r") as f:
    tag_dict = {}
    reader = csv.DictReader(f, delimiter='\t')
    for line in reader:
        old_tag = line['locus_tag']
        new_tag = "{}_{}_{}".format(BIN_PREFIX, line['CONTIG'].replace(line['SAMPLE'] + '_', ''), line['locus_tag'])
        tag_dict[old_tag] = new_tag
    locus_tag_prefix = line['locus_tag'].split("_")[0]


# Apply the substitutions to a file
def apply_substitutions(args):
    file_path, tag_dict = args
    # Read the file line by line and apply the substitutions.
    # Faster than reading the whole file into memory,
    # as the whole file content is copied multiple times when using the replace function.
    with open(file_path, "r") as f_in, open(file_path + ".tmp", "w") as f_out:
        for line in f_in:
            if locus_tag_prefix in line:
                # Find the exact locus tag with its underscore and digits
                locus_tag = re.search(r"{}_\d+".format(locus_tag_prefix), line).group(0)
                line = line.replace(locus_tag, tag_dict[locus_tag])
            f_out.write(line)
    os.rename(file_path + ".tmp", file_path)


# Find all files with the specified extensions and apply the substitutions to each file separately
file_paths = []
for root, dirs, files in os.walk("."):
    for file_name in files:
        if file_name.endswith((".faa", ".ffn", ".gbk", ".gff", ".sqn", ".tbl")):
            file_path = os.path.join(root, file_name)
            file_paths.append(file_path)

# Use a multiprocessing Pool to parallelize the calls to the apply_substitutions function
pool = Pool(processes=CPUS)
pool.map(apply_substitutions, [(file_path, tag_dict) for file_path in file_paths])
pool.close()
pool.join()
