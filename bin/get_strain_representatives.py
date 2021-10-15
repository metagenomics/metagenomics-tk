#!/usr/bin/env python3

import argparse as ap
import pandas as pd
import sys
from select_representative import get_representative

sys.dont_write_bytecode=True

def get_params(args):
    p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--genome_attributes', type=str, required=True,
                   help="Genome Attributes (CONTAMINATION, N50, ...)")
    p.add_argument('-c', '--genome_clusters', type=str, required=True,
                   help="Strain Clusters")
    p.add_argument('-o', '--out_dir', type=str, required=True,
                   help="Output Path")
    return p.parse_args(args)


def main(argv):
    args = get_params(argv)
    attributes = pd.read_csv(args.genome_attributes,  sep="\t")
    clusters = pd.read_csv(args.genome_clusters, sep="\t")

    merged_clusters = clusters.merge(attributes, on='BIN_ID', how='left')
    merged_clusters.STRAIN_CLUSTER = merged_clusters.CLUSTER.astype(str) + merged_clusters.STRAIN_CLUSTER.astype(str)
    representatives = list(map(lambda c: (c, get_representative(merged_clusters[merged_clusters.STRAIN_CLUSTER == c].copy())),
                               merged_clusters["STRAIN_CLUSTER"].unique()))

    # set representative column
    representatives_df = pd.DataFrame(list(representatives), columns=["REPRESENTATIVE_CLUSTER", "BIN_ID"])
    representatives_df["REPRESENTATIVE"] = 1
    representatives_df = pd.merge(merged_clusters, representatives_df, how="left", on="BIN_ID")
    representatives_df.REPRESENTATIVE.fillna(0, inplace=True)

    # rename columns to be compliant with the dereplication standard
    representatives_df[["BIN_ID", "STRAIN_CLUSTER", "REPRESENTATIVE"]]\
        .rename(columns={"BIN_ID": "GENOME", "STRAIN_CLUSTER": "CLUSTER"})\
        .to_csv(args.out_dir + "/strain_representatives.tsv", index=False, sep="\t")

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
