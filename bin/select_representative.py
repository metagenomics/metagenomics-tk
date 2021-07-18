#!/usr/bin/env python3

import argparse as ap
import pandas as pd
import numpy as np
import sys

sys.dont_write_bytecode=True

def get_params(args):
    p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--genome_attributes', type=str, required=True,
                   help="Genome Attributes (CONTAMINATION, N50, ...)")
    p.add_argument('-d', '--distance_matrix', type=str, required=True,
                   help="Distance Matrix")
    p.add_argument('-c', '--genome_clusters', type=str, required=True,
                   help="Clusters produced by fastcluster")
    p.add_argument('-o', '--out_dir', type=str, required=True,
                   help="Output Path")
    return p.parse_args(args)


def get_representative(genome_attributes):
    ranks = select(genome_attributes)
    ranks["SUM_RANKS"] = ranks.filter(regex=(".*RANK")).sum(axis=1)
    return ranks.sort_values("SUM_RANKS").iloc[0]["PATH"]


def select(attributes):
    attributes['CONTAMINATION_RANK'] = attributes['CONTAMINATION'].rank(ascending=True)
    attributes['COMPLETENESS_RANK'] = attributes['COMPLETENESS'].rank(ascending=False)
    attributes['COVERAGE_RANK'] = attributes['COVERAGE'].rank(ascending=False)
    attributes['N50_RANK'] = attributes['N50'].rank(ascending=False)
    attributes['HETEROGENEITY_RANK'] = attributes['HETEROGENEITY'].rank(ascending=True)
    return attributes


def write_representative(cluster_reps, path):
    with open(path, 'w') as fp:
        fp.write('\n'.join('{}\t{}'.format(x[0],x[1]) for x in cluster_reps))
        fp.write('\n')


def remove_duplicates(tuples):
    return list(set([tuple(sorted(list(x))) for x in tuples if x[0]!=x[1]]))


def get_closest_representative_comparisons(representatives, distances):
    genomes = [x[1] for x in representatives]
    filtered_distances = distances.filter(items=genomes,axis=1).filter(items=genomes, axis=0)
    np.fill_diagonal(filtered_distances.values.astype(float), np.nan )
    idx = filtered_distances.idxmin()
    ls = list(idx.to_dict().items())
    return remove_duplicates(ls)


def get_multiple_representative_comparisons(representatives, distances):
    genomes = [x[1] for x in representatives]
    filtered_distances = distances.filter(items=genomes,axis=1).filter(items=genomes, axis=0)
    triangle = np.triu(1-filtered_distances,k=+1)
    tr_df = 1-pd.DataFrame(triangle, columns=filtered_distances.columns.values, index=filtered_distances.index.values)
    tr_df_min = tr_df[tr_df < 1.1]
    tr_dr = tr_df_min.stack().dropna()

    return remove_duplicates(list(tr_dr.index))


def main(argv):
    args = get_params(argv)
    attributes = pd.read_csv(args.genome_attributes,  sep="\t")
    distances = pd.read_csv(args.distance_matrix, index_col=False, sep='\t', header=None)
    distances = distances.pivot(index=0, columns=1, values=2)
    clusters = pd.read_csv(args.genome_clusters, sep="\t")
    clusters_attr  = pd.merge(attributes, clusters, left_on=["PATH"], right_on=["GENOME"]).rename(columns={"GENOME_x": "GENOME"})
    representatives = list(map(lambda c: (c, get_representative(clusters_attr[clusters_attr.CLUSTER == c].copy())), clusters_attr["CLUSTER"].unique()))

    # set representative column
    representatives_df = pd.DataFrame(list(representatives), columns=["CLUSTER", "GENOME"])
    representatives_df["REPRESENTATIVE"] = 1
    representatives_df = pd.merge(clusters, representatives_df, how="left")
    representatives_df.REPRESENTATIVE.fillna(0, inplace=True)

    representatives_to_compare = get_closest_representative_comparisons(representatives, distances)
    representatives_df.to_csv(args.out_dir + "/representatives.tsv", index=False, sep="\t")
    write_representative(representatives_to_compare, args.out_dir + "/representatives_to_compare.tsv")

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
