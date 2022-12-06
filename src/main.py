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
    # number of bins = max diameter in psd
    # Each bin will be 1 A wide
    nbins = int(ceil(diameter.max()))
    print("nbins = ", nbins)
    plot_histogram(diameter, nbins = nbins)
    np.savetxt('diameters_of_geometric_points.dat', diameter)

    bin_freq, bin_edges, bin_index_of_points = make_histogram(diameter, nbins)
    # np.argsort returns the indices that would sort an array
    # bin_order is the order in which the bins are analyzed
    # bins are analyzed from the highest bin to the lowest
    bin_order = np.argsort(bin_freq)[::-1]
    bin_analyzed = []

    # Identify and analyze primary bins
    for b, boi in enumerate(bin_order):
        new_pore_type_found = cluster_bin(x, y, z, bin_index_of_points, boi)

        # if a secondary bin is found, it is more likely that other bins smaller than it
        # will also be secondary
        #if not new_pore_type_found:
            #break

    print('--------------------------')
    print('All Primary bins identified')

    all_cluster_pore_type_labels_flatten = [j for sub in config.all_cluster_pore_type_labels for j in sub]
    all_cluster_center_list_flatten = [j for sub in config.all_cluster_center_list for j in sub]
    all_cluster_diameter_list_flatten = [j for sub in config.all_cluster_diameter_list for j in sub]

    # Analyze secondary bins
    print('Analyzing the secondary bins')
    fill_pore_type_matrix()
    show_pore_type_matrix()


    pore_type_histogram = np.zeros((config.pore_type_count, nbins))
    for xi, yi, zi, di in zip(x, y, z, diameter):
        # which bin does each xyz point belong to
        bin_index = np.digitize(di, bin_edges)- 1

        # which pore type each xyz point belong to
        i, j, k = int(floor(xi)),int(floor(yi)), int(floor(zi))
        pore_type = int(config.pore_type_matrix_2[i, j, k])

        # store above two information in pore_type_histogram
        pore_type_histogram[pore_type-1][bin_index] += 1

    #print(pore_type_histogram)
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
    plot_psd_bar_pore_type(bin_centers, pore_type_histogram)
    show()

    sys.exit()