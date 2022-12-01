import config
import numpy as np
from clustering_functions import *
from cluster_bin import *
import time
import periodic_distance
import sys
#print(sys.path)
#print(sys.executable)
from plot_psd_bar_pore_type import plot_psd_bar_pore_type

if __name__ == "__main__":
    # x, y, z coordinates of the geometric probes and the diameter of largest spheres at those ponits
    x, y, z, diameter = read_vpsdpts(config.INPUT)
    nbins = 20
    print("nbins = ", nbins)
    plot_histogram(diameter, nbins = nbins)
    pore_type_histogram = np.zeros((config.Npores, nbins))
    bin_freq, bin_edges, bin_index_of_points = make_histogram(diameter, nbins)
    # np.argsort returns the indices that would sort an array
    # bin_order is the order in which the bins are executed
    bin_order = np.argsort(bin_freq)[::-1]
    bin_analyzed = []
    for b, boi in enumerate(bin_order):
        # xk, yk, zk are the points inside bin_number = boi
        if config.pore_type_count < config.Npores-1:
            #xk, yk, zk, pore_type_label = \
            new_pore_type_found = cluster_bin(x, y, z, bin_index_of_points, boi, compute_pbc_matrix = 1)

            if new_pore_type_found: # bins which coorespond to new pore types
                bin_analyzed.append(boi)

    all_cluster_pore_type_labels_flatten = [j for sub in config.all_cluster_pore_type_labels for j in sub]
    all_cluster_center_list_flatten = [j for sub in config.all_cluster_center_list for j in sub]
    all_cluster_diameter_list_flatten = [j for sub in config.all_cluster_diameter_list for j in sub]

    update_pore_type_matrix_2()
    show_pore_type_matrix_2()


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


    # After All pore type centers identified
    for b, boi in enumerate(bin_order):
        if boi not in bin_analyzed:
            # For the rest of the bins, assign each point to the pore type of nearest cluster center weighted by pore size

            for i, index in enumerate(bin_index_of_points):
                if index == b:
                    # distance of the point from each cluster center
                    distance_from_cluster_center = []

                    # calculate distance of the point to each of the cluster centers
                    for j, cluster_center in enumerate(all_cluster_center_list_flatten):
                        # distance weighted by the pore size
                        distance_from_cluster_center.append(distance(np.array([x[i], y[i], z[i]]), cluster_center)
                                                            /all_cluster_diameter_list_flatten[j])

                    # index of the distance_from_cluster_center list where the distance is minimum
                    distance_from_cluster_center = np.array(distance_from_cluster_center)
                    index_min = np.argmin(distance_from_cluster_center)

                    # Converting the coordinate of the points to integer
                    xint, yint, zint = int(floor(x[i])), int(floor(y[i])), int(floor(z[i]))
                    # update the pore_type_matrix by the pore type of the cluster which is closest to the point
                    config.pore_type_matrix[xint, yint, zint] = all_cluster_pore_type_labels_flatten[index_min]
#%%
    from importlib import reload
    import clustering_functions
    reload(clustering_functions)
    from clustering_functions import show_pore_type_matrix
    show_pore_type_matrix()

#%%
    #print(len(xk))
    """
    from periodic_distance import distance_matrix_periodic as cython_periodic
    from importlib import reload
    import sys
    reload(periodic_distance)
    import periodic_distance
    importlib.reload(sys.modules['periodic_distance'])
    from periodic_distance import distance_matrix_periodic as cython_periodic
    start = time.time()
    pbc_matrix_cython = cython_periodic(xk, yk, zk)
    end = time.time()
    print(end-start)
    """
