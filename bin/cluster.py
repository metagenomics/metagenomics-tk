#!/usr/bin/env python3

import argparse as ap
import sys

from fastcluster import *
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

sys.dont_write_bytecode=True


def get_params(args):
    p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--distance_matrix', type=str, required=True,
                   help="Distance Matrix")
    p.add_argument('-o', '--out_dir', type=str, required=True,
                   help="Output Path")
    p.add_argument('-c', '--cutoff', type=float, required=False,
                   help="Mash distance cutoff", default=0.05)
    p.add_argument('-p', '--png_path', type=str, required=False,
                   help="PNG Path")
    return p.parse_args(args)


def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = hierarchy.dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata


def run_clustering(mat, mash_distance_cutoff = 0.05):
    mat = average(scipy.spatial.distance.squareform(mat))
    final_clusters = hierarchy.fcluster(mat, mash_distance_cutoff, criterion='distance')
    return mat, final_clusters


def build_dendrogram(mat, path, labels):
    plt.figure()
    plt.gcf().subplots_adjust(bottom=0.45)
    fancy_dendrogram(mat, labels=labels, leaf_rotation=90, max_d=0.05)
    plt.savefig(path)


def write_file(mat, path):
    np.savetxt(path, mat, delimiter='\t')


def write_clusters(clusters, labels, output_path):
    header = "CLUSTER\tGENOME\n"
    outputStr = header
    for i in range(0, len(clusters)):
        outputStr += ' '.join(map(str, clusters[i].flatten())) + "\t" + str(labels[i]) + "\n"

    with open(output_path, "w") as output:
        output.write(outputStr)


def main(argv):
    args = get_params(argv)
    mat1 = pd.read_csv(args.distance_matrix, index_col=False, sep='\t', header=None)
    distances = mat1.pivot(index=0, columns=1, values=2)
    dendrogram_mat, clusters = run_clustering(distances.values, args.cutoff)
    write_file(dendrogram_mat, args.out_dir + "/distances.tsv")
    write_file(clusters, args.out_dir + "/clusters.txt")
    write_clusters(clusters, distances.index.values, args.out_dir + "/clusters.tsv")
    if args.png_path is not None:
        build_dendrogram(dendrogram_mat, args.png_path, distances.index.values)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
