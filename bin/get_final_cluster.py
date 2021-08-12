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
                   help="Clusters produced by fastcluster")
    p.add_argument('-r', '--representative_comparison', type=str, required=True,
                   help="An ANI comparison of all representatives")
    p.add_argument('-a', '--representative_ani', type=float, default=0.95, required=False,
                   help="ANI cutoff for representative comparison (Default > 0.95)")
    p.add_argument('-o', '--out_dir', type=str, required=True,
                   help="Output Path")
    return p.parse_args(args)


def update_cluster(clusters, attributes, genome_a, genome_b):
    cluster_representative_b = clusters[clusters["CLUSTER"] == clusters[clusters["GENOME"] == genome_b].CLUSTER.iloc[0]]
    cluster_representative_a = clusters[clusters["CLUSTER"] == clusters[clusters["GENOME"] == genome_a].CLUSTER.iloc[0]]

    genome_rep_a = cluster_representative_a[ cluster_representative_a["REPRESENTATIVE"] == 1.0].GENOME.iloc[0]
    genome_rep_b = cluster_representative_b[ cluster_representative_b["REPRESENTATIVE"] == 1.0].GENOME.iloc[0]

    cluster_rep_a = cluster_representative_a.iloc[0].CLUSTER
    cluster_rep_b = cluster_representative_b.iloc[0].CLUSTER

    cluster_a = clusters[clusters["CLUSTER"] == cluster_rep_a]
    cluster_b = clusters[clusters["CLUSTER"] == cluster_rep_b]

    clusters.loc[clusters.CLUSTER == cluster_b.CLUSTER.iloc[0],"CLUSTER"] = cluster_a.CLUSTER.iloc[0]
    rep_a_attr = attributes.loc[attributes.BIN_ID == genome_rep_a,:]
    rep_b_attr = attributes.loc[attributes.BIN_ID == genome_rep_b,:]

    rep_genome = get_representative(pd.concat([rep_a_attr, rep_b_attr]))

    if rep_genome == genome_rep_a:
        clusters.loc[clusters.GENOME == genome_rep_b, "REPRESENTATIVE"] = 0
    else:
        clusters.loc[clusters.GENOME == genome_rep_a, "REPRESENTATIVE"] = 0
    clusters.loc[clusters.GENOME == rep_genome,"REPRESENTATIVE"] = 1

    return clusters


def main(argv):
    args = get_params(argv)
    attributes = pd.read_csv(args.genome_attributes,  sep="\t")
    clusters = pd.read_csv(args.genome_clusters, sep="\t")
    representatives = pd.read_csv(args.representative_comparison, sep="\t")
    similar_representatives = representatives[representatives.ANI > args.representative_ani]

    for i in range(0, similar_representatives.shape[0]):
        genome_a = similar_representatives.iloc[i]["GENOME_A"]
        genome_b = similar_representatives.iloc[i]["GENOME_B"]
        clusters = update_cluster(clusters, attributes, genome_a, genome_b)

    clusters.to_csv(args.out_dir + "/representatives.tsv",sep='\t', index=False)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
