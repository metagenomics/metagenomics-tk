#!/usr/bin/env python3

import multiprocessing as mp
import sys
import csv

kegg = dict()

with open (sys.argv[2]+"/genes_ko.list", "r") as koFile:
    read_tsv = csv.reader(koFile, delimiter="\t")
    for row in read_tsv:
        if row[0] not in kegg:
            kegg[row[0]] = []
        kegg[row[0]].append(row[1])


with open (sys.argv[2]+"/genes_pathway.list", "r") as pathFile:
    read_tsv = csv.reader(pathFile, delimiter="\t")
    for row in read_tsv:
        if row[0] not in kegg:
            kegg[row[0]] = []
        kegg[row[0]].append(row[1])

def searchAppend(line, fieldNames, writer):
    result = [line["query"],line["target"]]
    ko = []
    pathway = []

    if line["target"] in kegg:
        for elem in kegg[line["target"]]:
            if elem.startswith("ko:"):
                ko.append(elem)
            if elem.startswith("path:"):
                pathway.append(elem)
    if not ko:
        ko.append("NONE")
    if not pathway:
        pathway.append("NONE")
    
    result.append(';'.join(ko))
    result.append(';'.join(pathway))
    result.extend([line[f] for f in fieldNames])
    writer.writerow(result)

blast_input = open(sys.argv[1])
read_tsv = csv.DictReader(blast_input, delimiter="\t")


outputHeader = ['GENE','KEGG','KO','PATHWAY'] + read_tsv.fieldnames

with open(sys.argv[3], 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(outputHeader)
    for row in read_tsv:
        searchAppend(row, read_tsv.fieldnames, tsv_writer)
    blast_input.close()
