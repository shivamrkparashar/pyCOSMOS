import numpy as np
from pylab import *
# for 3d plotting
import plotly.graph_objs as go
import plotly.express as px
import plotly.io as pio
pio.renderers.default ="browser"
import pandas as pd
from sklearn.cluster import *
from math import sqrt
from numpy import floor, ceil
import config

def plot_histogram(arr, nbins):
    figure()
    # The min and max (range) of the histogram are the nearest integers given by floor or ceil
    range = (floor(arr.min()), ceil(arr.max()))
    hist(arr, density =True, bins = nbins, range = range, color ='orange', edgecolor='black', linewidth=1)
    xlabel('Diameter ($\AA$)')
    ylabel('Normalized Distribution')
    xlim(0,)
    ylim(0,)
    savefig('PSD.jpg', format='jpg',dpi=300 ,bbox_inches='tight')
    np.savetxt('psd_diameter.dat', arr)
def read_vpsdpts(INPUT):
    """
    Reads the vpsdpts file generated by Zeo++
    
    Parameters
    ------------------------
    INPUT: string
    filename with vpsdpts extension
    
    Returns
    ------------------------
    x: numpy array
    y: numpy array
    z: numpy array
    diameter: numpy array
    """
    x, y, z, diameter = [], [], [], []
    linenumber = 0

    with open(INPUT) as inp:
        for line in inp:
            linenumber += 1
            if linenumber > 2:
                s = line.split()
                x.append(float(s[1]))
                y.append(float(s[2]))
                z.append(float(s[3]))
                diameter.append(2*float(s[4]))

    # Convert lists into numpy arrays
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    diameter = np.array(diameter)

    return x, y, z, diameter

def make_histogram(diameter_arr, nbins):
    """
    Makes a histogram from an array and number of bins

    Parameters
    ------------------------
    diameter_arr: numpy array
    nbins: integer

    Returns
    ------------------------
    bin_freq: numpy array
    bin_edges: numpy array
    bin_index: numpy array
        Tells in which bin does each of the point goes
    """
    # The min and max (range) of the histogram are the nearest integers given by floor or ceil
    range = (floor(diameter_arr.min()), ceil(diameter_arr.max()))
    bin_freq, bin_edges =  np.histogram(diameter_arr, bins =nbins, range = range) # Length nbins + 1
    #bin_freq, bin_edges =  np.histogram(diameter_arr, bins =nbins, range = (0, np.ceil(diameter_arr.max()))) # Length nbins + 1

    bin_index = np.zeros(len(diameter_arr))
    for i, di in enumerate(diameter_arr):
        for j, _ in enumerate(bin_edges):
            if di > bin_edges[j] and di <= bin_edges[j+1]:
                bin_index[i] = j

    return bin_freq, bin_edges, bin_index

def approx_equal(lst1, lst2):
    """
    Checks if the two lists are approximately equal 
    
    Parameters
    ------------------------
    lst1: list
    lst2: list
    
    lst1 and lst2 are like [c1, c2]
    where ci = np.array([xi, yi, zi])

    Returns
    ------------------------
    True or False
    """
    #result = [] # len(lst1) can be True or False depending whether each element of lst1 is close to lst2
    
    if len(lst1) == len(lst2):
        for i in range(len(lst1)):
            result = np.allclose(lst1[i], lst2[i], atol = 1)
            if not result: # if any element is false
                return False
            
        # If all elements are True
        return True
    else:
        #print('length not equal')
        return False

def distance(v1, v2):
    """
    Calculated distance between two vectors v1 and v2.
    
    Parameters
    ------------------------
    v1: numpy array
    v2: numpy array
    
    Returns
    ------------------------
    distance: double 
    """
    dx = v1[0] - v2[0]
    dy = v1[1] - v2[1]
    dz = v1[2] - v2[2]
    
    if abs(dx) > config.Lx/2:
        dx = config.Lx - abs(dx)
    
    if abs(dy) > config.Ly/2:
        dy = config.Ly - abs(dy)
        
    if abs(dz) > config.Lz/2:
        dz = config.Lz - abs(dz)
    
    squared_distance = dx**2 + dy**2 + dz**2
    
    return sqrt(squared_distance)

def distance_matrix_periodic(x_arr, y_arr, z_arr):
    
    # Precompute distances using periodic boundary conditions
    xyz_matrix = np.array([x_arr, y_arr, z_arr]).T
    periodic_distance_matrix = np.zeros((shape(xyz_matrix)[0], shape(xyz_matrix)[0]))
    
    for i, vi in enumerate(xyz_matrix):
        for j,vj in enumerate(xyz_matrix):
            periodic_distance_matrix[i][j] = distance(vi, vj) 
            
    return periodic_distance_matrix

def update_pore_type_matrix(x_arr, y_arr, z_arr, pore_type_labels):
    """
    pore_type_label = can be range(Npores) = 0 or 1
    Assigns the pore_type_labels array into a 3d matrix grid
    
    """
    for i, li in enumerate(pore_type_labels):
        #for _ in np.unique(pore_type_labels):
        xi, yi, zi = int(floor(x_arr[i])), int(floor(y_arr[i])), int(floor(z_arr[i]))
        config.pore_type_matrix[xi, yi, zi] = li
        # updating the pore type matrix according to the bins first
        # So that the grid points corresponging to the particles belonging in these bins
        # are not decideded based on the distance from the cluster surfaces
        config.pore_type_matrix_2[xi, yi, zi] = li

def update_pore_type_matrix_2():
    """
    based on only the cluster centers and cluster sizes
    does not require information of other bins
    :return:
    """
    print('Updating pore type matrix')
    print('-------------------------')

    all_cluster_pore_type_labels_flatten = [j for sub in config.all_cluster_pore_type_labels for j in sub]
    all_cluster_center_list_flatten = [j for sub in config.all_cluster_center_list for j in sub]
    all_cluster_diameter_list_flatten = [j for sub in config.all_cluster_diameter_list for j in sub]

    x = np.arange(0.5, (int(ceil(config.Lx))) + 0.5)
    y = np.arange(0.5, (int(ceil(config.Ly))) + 0.5)
    z = np.arange(0.5, (int(ceil(config.Lz))) + 0.5)

    for i, xi in enumerate(x):
        for j, yi in enumerate(y):
            for k, zi in enumerate(z):

                # proceed only if the pore type matrix has -1
                # This is done to make sure the grid points calculated using the largest and
                # second largest bins, and third largest bins do not get overwritten
                #if config.pore_type_matrix_2[i, j, k] != -1: continue


                # calculate distance of the point to each of the cluster centers
                distance_from_cluster_center = []
                for cc, cluster_center in enumerate(all_cluster_center_list_flatten):
                    # distance weighted by the pore size
                    #distance_from_cluster_center.append(distance(np.array([xi, yi, zi]), cluster_center))
                    #distance_from_cluster_center.append(distance(np.array([xi, yi, zi]), cluster_center)
                    #                                    /all_cluster_diameter_list_flatten[cc])
                    # distance from cluster surface
                    distance_from_cluster_center.append(distance(np.array([xi, yi, zi]), cluster_center)
                                                        -0.5*all_cluster_diameter_list_flatten[cc])

                # index of the distance_from_cluster_center list where the distance is minimum
                distance_from_cluster_center = np.array(distance_from_cluster_center)
                index_min = np.argmin(distance_from_cluster_center)

                config.pore_type_matrix_2[i][j][k] = all_cluster_pore_type_labels_flatten[index_min]
                config.pore_type_matrix_3[i][j][k] = index_min



def count_pore_type_labels(i, j, k):
    """
    Counts the number of neighbors with pore_type = 0,1, or 2 around i, j, k index
    """
    # Nx, Ny, Nz are the numnber of elements in x, y, and z dirctions
    Nx, Ny, Nz = shape(config.pore_type_matrix)
    count_pore_0, count_pore_1, count_pore_2 = 0, 0, 0

    for di in [0, -1, 1]:
        for dj in [0, -1, 1]:
            for dk in [0, -1, 1]:
                ni = i + di
                nj = j + dj
                nk = k + dk

                # if the new index exceeds the size of the matrix, the new index is 0
                if ni == Nx: ni = 0
                if nj == Ny: nj = 0
                if nk == Nz: nk = 0

                if config.pore_type_matrix[ni][nj][nk] == 0:
                    count_pore_0 += 1
                if config.pore_type_matrix[ni][nj][nk] == 1:
                    count_pore_1 += 1
                if config.pore_type_matrix[ni][nj][nk] == 2:
                    count_pore_2 += 1

    return count_pore_0, count_pore_1, count_pore_2
def smoothen_pore_type_matrix():

    X = range(0, (int(floor(config.Lx))))
    Y = range(0, (int(floor(config.Ly))))
    Z = range(0, (int(floor(config.Lz))))

    for i in X:
        for j in Y:
            for k in Z:
                if config.pore_type_matrix[i][j][k] == -1:
                    count_pore_0, count_pore_1, count_pore_2 = count_pore_type_labels(i, j, k)

                    if count_pore_0 > 9:
                        config.pore_type_matrix[i][j][k] = 0

                    if count_pore_1 > 9:
                        config.pore_type_matrix[i][j][k] = 1

                    if count_pore_2 > 9:
                        config.pore_type_matrix[i][j][k] = 2


def show_pore_type_matrix():
    """
    plot x, y, z coordinates in 3d
    Assumed grid of 1A width in all directions
    Also saves the pore type matrix in a csv file
    :return:
    """
    X = range(0, (int(floor(config.Lx))))
    Y = range(0, (int(floor(config.Ly))))
    Z = range(0, (int(floor(config.Lz))))

    xx, yy, zz, pore_type_label_1d = [], [], [], []
    for i in X:
        for j in Y:
            for k in Z:
                xx.append(i)
                yy.append(j)
                zz.append(k)
                pore_type_label_1d.append(config.pore_type_matrix[i][j][k])

    #xx, yy, zz = np.meshgrid(X, Y, Z)
    # Smoothen the pore type matrix
    smoothen_pore_type_matrix()

    # Here .flatten() vectorizes the coordinates column-wise.
    #d = {'x':xx.flatten(), 'y':yy.flatten(), 'z':zz.flatten(), 'color':config.pore_type_matrix.flatten()}
    d = {'x':xx, 'y':yy, 'z':zz, 'color':pore_type_label_1d}
    df = pd.DataFrame(data = d)
    #df.to_csv(config.File_pore_type_matrix, index =False)

    df["color"] = df["color"].astype(str)
    fig = px.scatter_3d(df, x = 'x', y = 'y', z='z', color = 'color',)
    fig.update_traces(marker=dict(size=6,
                         line=dict(width=0.5,
                         #color='Black')),
                      color='Black')),
                    selector=dict(mode='markers'),
                      marker_symbol = 'square')

    fig.show()

def show_pore_type_matrix_2():
    """
    plot x, y, z coordinates in 3d
    Assumed grid of 1A width in all directions
    Also saves the pore type matrix in a csv file
    :return:
    """
    X = range(0, (int(ceil(config.Lx))))
    Y = range(0, (int(ceil(config.Ly))))
    Z = range(0, (int(ceil(config.Lz))))

    xx, yy, zz, pore_type_label_1d = [], [], [], []
    cluster_center_label_1d = []
    for i in X:
        for j in Y:
            for k in Z:
                xx.append(i)
                yy.append(j)
                zz.append(k)
                pore_type_label_1d.append(config.pore_type_matrix_2[i][j][k])
                cluster_center_label_1d.append(config.pore_type_matrix_3[i][j][k])

    d = {'x':xx, 'y':yy, 'z':zz, 'color':pore_type_label_1d}
    df = pd.DataFrame(data = d)
    df.to_csv(config.File_pore_type_matrix, index =False)

    d = {'x':xx, 'y':yy, 'z':zz, 'color':cluster_center_label_1d}
    df = pd.DataFrame(data = d)
    df.to_csv('pore_type_matrix_with_cluster_center_label.csv', index =False)

    df["color"] = df["color"].astype(str)
    fig = px.scatter_3d(df, x = 'x', y = 'y', z='z', color = 'color',)
    fig.update_traces(marker=dict(size=6,
                                  line=dict(width=0.5,
                                            #color='Black')),
                                            color='Black')),
                      selector=dict(mode='markers'),
                      marker_symbol = 'square')

    fig.show()
def dbscan(x_arr, y_arr, z_arr,periodic_distance_matrix, eps, min_samples):
    """
    DBSCAN to cluster points
    """
    sol = DBSCAN(eps = eps, min_samples = min_samples,metric = "precomputed", n_jobs =8).fit(periodic_distance_matrix)
    labels = sol.labels_ # can be -1, 0, 1, 2 ...

    fraction_of_noisy_points = np.count_nonzero(labels == -1)/len(labels)

    # Number of clusters
    #Ncluster = max(np.unique(labels)) + 1
    
    # calculate cluster centers
    cluster_center_list, cluster_diameter_list = best_cluster_center(x_arr, y_arr, z_arr, labels)

    # plot x, y, z coordinates in 3d
    if fraction_of_noisy_points < 0.5:
        d = {'x':x_arr, 'y':y_arr, 'z':z_arr, 'color':labels}
        df = pd.DataFrame(data = d)
        df["color"] = df["color"].astype(str)
        fig = px.scatter_3d(df, x = 'x', y = 'y', z='z', color = "color")
        fig.update_traces(marker=dict(size=4,
                             line=dict(width=0.8,
                             color='Black')),
                  selector=dict(mode='markers'))
        fig.show()
    return labels, cluster_center_list, cluster_diameter_list, fraction_of_noisy_points

def best_cluster_center(x_arr, y_arr, z_arr, labels):
    """
    Inputs
    --------------------------
    x_arr, y_arr, z_arr: numpy array of data points
    labels: numpy array of length x_arr. Element of this array can be 0, 1, 2 ...

    Calculates the center of the cluster
    """
    cluster_center_list = []
    cluster_diameter_list = []

    Ncluster = max(np.unique(labels))+1
    #print('Number of clusters = ', Ncluster)

    if Ncluster ==0:
        print('All data points classified as noise')
        return

    # Initialize the center of each cluster
    # x_center = [x_cluster_0, x_cluster_1 ....]
    x_center = np.zeros(Ncluster)
    y_center = np.zeros(Ncluster)
    z_center = np.zeros(Ncluster)

    # Center = sum over all the points with same labels
    #for i, li in enumerate(labels):
    for x, y, z, li in zip(x_arr, y_arr, z_arr, labels):
        if li != -1:  # Be very careful, as label = -1 means noise. But x_center[-1] = last element
            x_center[li] += x
            y_center[li] += y
            z_center[li] += z

    # Divide center by the number of points with same labels
    for li in range(Ncluster):
        x_center[li] /= np.count_nonzero(labels==li)
        y_center[li] /= np.count_nonzero(labels==li)
        z_center[li] /= np.count_nonzero(labels==li)
        #print('Cluster centers %1.3g, %1.3g, %1.3g' %(x_center[li], y_center[li], z_center[li]))

    """
    # To calculate actual center of the cluster
    For each cluster:
        For each mirror image of center:
            calculate avg_distance_from_center
    """

    for li in range(Ncluster):
        sum_distance_from_center_list = [] # one element for each image of cluster center

        for dx in [0, config.Lx/2, -config.Lx/2]:
            for dy in [0, config.Ly/2, -config.Ly/2]:
                for dz in [0, config.Lz/2, -config.Lz/2]:
                    x_new_center = x_center[li] + dx
                    y_new_center = y_center[li] + dy
                    z_new_center = z_center[li] + dz

                    # new center should be inside the unit cell
                    x_new_center, y_new_center, z_new_center = put_point_in_box(x_new_center, y_new_center, z_new_center)
                    new_center = np.array([x_new_center, y_new_center, z_new_center])

                    # the sum of squares of the distance
                    sum_distance2_from_center = 0
                    for j, lj in enumerate(labels):
                        if lj == li:
                            point = np.array([x_arr[j], y_arr[j], z_arr[j]])
                            sum_distance2_from_center += distance(new_center, point)**2
                    sum_distance_from_center_list.append(sum_distance2_from_center)

        index_min = np.argmin(sum_distance_from_center_list)
        #index_min = min(range(len(sum_distance_from_center_list)), key=sum_distance_from_center_list.__getitem__)

        #print(li, index_min)
        index = 0
        #Npoints_in_cluster = 0
        for dx in [0, config.Lx/2, -config.Lx/2]:
            for dy in [0, config.Ly/2, -config.Ly/2]:
                for dz in [0, config.Lz/2, -config.Lz/2]:
                    if index == index_min:
                        xx = x_center[li] + dx
                        yy = y_center[li] + dy
                        zz = z_center[li] + dz

                        # new center should be inside the unit cell
                        xx, yy, zz = put_point_in_box(xx, yy, zz)
                        # cluster_diameter based on radius of gyration
                        cluster_diameter = 2*sqrt(sum_distance_from_center_list[index]/np.count_nonzero(labels==li))

                        cluster_center_list.append(np.array([xx, yy, zz]))
                        cluster_diameter_list.append(cluster_diameter)
                        #print('Original center = %1.3g, %1.3g, %1.3g' %(x_center[li], y_center[li], z_center[li]))
                        print('-----------Cluster # %d ---------------' %li)
                        print('Center  = %1.3g, %1.3g, %1.3g' %(xx, yy, zz))
                        print('Average size = %1.3g A' %cluster_diameter)
                        print('---------------------------------------')

                    index += 1


        cluster_center = cluster_center_list[-1]
        radius_list = [] # distance of each point within a cluster to its center
        for x, y, z, lj in zip(x_arr, y_arr, z_arr, labels):
            if lj == li: # only those points that belong to the cluster
                radius = distance(cluster_center, np.array([x, y, z]))
                radius_list.append(radius)
        """        
        figure()
        hist(radius_list, density =True, color ='magenta', edgecolor='black')
        title('Pore type = %d, Cluster number = %d' %(config.pore_type_count+1, li+1))
        xlabel('Distance from cluster center ($\AA$)')
        ylabel('Normalized distribution')
        show()
        """

    return(cluster_center_list, cluster_diameter_list)

def put_point_in_box(x, y, z):
    if x > config.Lx:
        x -= config.Lx
    if x < 0:
        x += config.Lx

    if y > config.Ly:
        y -= config.Ly
    if y < 0:
        y += config.Ly

    if z > config.Lz:
        z -= config.Lz
    if z < 0:
        z += config.Lz

    return x, y, z