#!/usr/bin/env python3

import multiprocessing as mp
import sys
import csv

kegg = dict()
outputFile = []
outputFile.append(['GENE','KEGG','KO','PATHWAY'])

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

def searchAppend(line):

    result = [line[0],line[1]]
    ko = []
    pathway = []

    if line[1] in kegg:
        for elem in kegg[line[1]]:
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
    outputFile.append(result)



diamond_input = open(sys.argv[1])
read_tsv = csv.reader(diamond_input, delimiter="\t")

for row in read_tsv:
    searchAppend(row)
diamond_input.close()

with open(sys.argv[3], 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    for elem in outputFile:
        tsv_writer.writerow(elem)