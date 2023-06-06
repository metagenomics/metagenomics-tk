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

# Read the BIN_PREFIX value from the command line arguments
BIN_PREFIX = sys.argv[1]
locus_tag_prefix = ""

# Read the prokka_tags_tmp.tsv file and generate a dictionary with old tags as keys and new tags as values
with open("prokka_tags_tmp.tsv", "r") as f:
    next(f)  # skip header
    tag_dict = {}
    for line in f:
        fields = line.strip().split("\t")
        old_tag = fields[2]
        new_tag = "{}_{}_{}".format(BIN_PREFIX, fields[1].replace(fields[0] + '_', ''), fields[2])
        tag_dict[old_tag] = new_tag
    locus_tag_prefix = fields[2].split("_")[0]


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
pool = Pool()
pool.map(apply_substitutions, [(file_path, tag_dict) for file_path in file_paths])
pool.close()
pool.join()
