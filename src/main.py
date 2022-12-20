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
    a, b, c, geom_dia = read_vpsdpts(config.INPUT)

    # number of bins = max diameter in psd
    # Each bin will be 1 A wide
    nbins = int(ceil(geom_dia.max())- floor(geom_dia.min()))
    print("nbins = ", nbins)
    plot_histogram(geom_dia, nbins = nbins)
    np.savetxt('diameters_of_geometric_points.dat', geom_dia)

    bin_freq, bin_edges, bin_index_of_points = make_histogram(geom_dia, nbins)
    # np.argsort returns the indices that would sort an array
    # bin_order is the order in which the bins are analyzed
    # bins are analyzed from the highest bin to the lowest
    bin_order = np.argsort(bin_freq)[::-1]
    bin_analyzed = []

    # Identify and analyze primary bins
    for boi in bin_order:

        # If the number of pore types identified = NumberOfPoreTypes, break the loop
        if config.pore_type_count == config.NumberOfPoreTypes:
            break

        # 1. Extract the points within bin boi
        ak, bk, ck, diak = extract_points_in_a_bin(a, b, c, geom_dia, bin_index_of_points, boi)

        # 2. Cluster the points within bin boi
        Ncluster, Npoints_per_cluster, cluster_labels, cluster_center_list, cluster_diameter_list, fraction_of_noisy_points =\
        cluster_points_within_a_bin(ak, bk, ck)

        # 3, Classify the points as primary or secondary
        #primary_bin = 1
        primary_bin = classify_bin(diak, Npoints_per_cluster, fraction_of_noisy_points, boi)

        if primary_bin:

            config.pore_type_count += primary_bin
            print('New pore type identified, pore type count = ', config.pore_type_count)
            config.all_cluster_center_list.append(cluster_center_list)
            config.all_cluster_diameter_list.append(cluster_diameter_list)
            config.all_cluster_pore_type_labels.append(Ncluster * [config.pore_type_count])

            # update the pore type matrix
            update_pore_type_matrix(ak, bk, ck, cluster_labels)

            # x, y, z coordinates in 3d
            xk, yk, zk = abc_to_xyz(ak, bk, ck)
            # save xyz file
            save_as_xyz(boi, xk, yk, zk, cluster_labels)
            # plot x, y, z coordinates
            plot_xyz(xk, yk, zk, cluster_labels)
        else:
            print('%d is a secondary bin' %boi)

        # new_pore_type_found = cluster_bin(a, b, c, geom_dia, bin_index_of_points, boi)


        # Comment below is not true
        # if a secondary bin is found, it is more likely that other bins smaller than it
        # will also be secondary
        #if not new_pore_type_found:
            #break

    print('--------------------------')
    print('All Primary bins identified')
    print('Number of pore types = %d' %config.pore_type_count)

    if config.pore_type_count == 0:
        print('Exiting, change epsilon and minpts params of dbscan')
        exit(1)


    #show_pore_type_matrix()

    # Analyze secondary bins
    print('Analyzing the secondary bins')
    fill_pore_type_matrix()
    show_pore_type_matrix()


    pore_type_histogram = np.zeros((config.pore_type_count, nbins))

    for xi, yi, zi, di in zip(a, b, c, geom_dia):

        # which bin does each xyz point belong to
        bin_index = np.digitize(di, bin_edges)- 1

        # which pore type each xyz point belong to
        i, j, k = int(floor(xi)),int(floor(yi)), int(floor(zi))
        pore_type = int(config.pore_type_matrix_with_pore_type_labels[i, j, k])

        # store above two information in pore_type_histogram
        pore_type_histogram[pore_type-1][bin_index] += 1

    #print(pore_type_histogram)
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
    plot_psd_bar_pore_type(bin_centers, pore_type_histogram)
    show()

    sys.exit()