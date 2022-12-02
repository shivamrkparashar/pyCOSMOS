import config
import numpy as np
from clustering_functions import *
from cluster_bin import *
import time
import periodic_distance
import sys
from plot_psd_bar_pore_type import plot_psd_bar_pore_type

if __name__ == "__main__":
    # x, y, z are the geometric points and diameter is their  geometric diameters
    x, y, z, diameter = read_vpsdpts(config.INPUT)
    nbins = 20
    print("nbins = ", nbins)
    plot_histogram(diameter, nbins = nbins)
    np.savetxt('diameters_of_geometric_points.dat', diameter)

    pore_type_histogram = np.zeros((config.Npores, nbins))
    bin_freq, bin_edges, bin_index_of_points = make_histogram(diameter, nbins)
    # np.argsort returns the indices that would sort an array
    # bin_order is the order in which the bins are analyzed
    bin_order = np.argsort(bin_freq)[::-1]
    bin_analyzed = []

    for b, boi in enumerate(bin_order):
        # xk, yk, zk are the points inside bin_number = boi
        if config.pore_type_count < config.Npores-1:
            #xk, yk, zk, pore_type_label = \
            new_pore_type_found = cluster_bin(x, y, z, bin_index_of_points, boi)

            if new_pore_type_found: # bins which coorespond to new pore types
                bin_analyzed.append(boi)

    all_cluster_pore_type_labels_flatten = [j for sub in config.all_cluster_pore_type_labels for j in sub]
    all_cluster_center_list_flatten = [j for sub in config.all_cluster_center_list for j in sub]
    all_cluster_diameter_list_flatten = [j for sub in config.all_cluster_diameter_list for j in sub]

    fill_pore_type_matrix()
    show_pore_type_matrix()


    for xi, yi, zi, di in zip(x, y, z, diameter):
        # which bin does each xyz point belong to
        bin_index = np.digitize(di, bin_edges)- 1

        # which pore type each xyz point belong to
        i, j, k = int(floor(xi)),int(floor(yi)), int(floor(zi))
        pore_type = int(config.pore_type_matrix_2[i, j, k])

        # store above two information in pore_type_histogram
        pore_type_histogram[pore_type][bin_index] += 1

    #print(pore_type_histogram)
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
    plot_psd_bar_pore_type(bin_centers, pore_type_histogram)
    show()

    sys.exit()