import numpy as np
import config
from clustering_functions import periodic_distance, abc_to_xyz, put_point_in_box, abc_to_xyz_arrays

def save_cluster_as_xyz(bin_number, cluster_number, x_arr, y_arr, z_arr, name=''):
    """
    Saves the x, y, z arrays as an xyz file for visualization in vmd
    A different element is chosen for each clusters

    :param boi: bin number
    :param x_arr:
    :param y_arr:
    :param z_arr:
    :return:
    """
    xyz_file = '%sbin_%d_cluster_%d.xyz' %(name, bin_number, cluster_number)
    Natoms = len(x_arr)

    with open(xyz_file, 'w') as out:
        out.write('%d \n' %Natoms)
        out.write('bin # %d, cluster # %d \n' %(bin_number, cluster_number))
        for xi, yi, zi in zip(x_arr, y_arr, z_arr):
            out.write("0 %1.3f %1.3f %1.3f\n" %(xi, yi, zi))

def classify_cluster_shape(bin_number, cluster_number, ac, bc, cc, center):
    """
    Classifies the cluster to be either spherical or channel type
    The input cluster coordinates are centered at the unit cell

    :param bin_number:
    :param cluster_number:
    :param ac:
    :param bc:
    :param cc:
    :param center:
    :return:
    """
    x, y, z = abc_to_xyz_arrays(ac, bc, cc)
    save_cluster_as_xyz(bin_number, cluster_number, x, y, z)

    orientation = np.array([0, 0, 0])

    ac, bc, cc = centralize_cluster(ac, bc, cc, center)

    # check if the cluster of points is a channel
    atol = 1
    if np.isclose(np.amin(ac), 0, atol=atol) and np.isclose(np.amax(ac), config.Lx, atol=atol):
        hist, bin_edges = np.histogram(ac, bins=20, density=True)
        if np.any(np.absolute(hist) < 1e-2):
            orientation += np.array([0, 0, 0])
        else:
            orientation += np.array([1, 0, 0])

    if np.isclose(np.amin(bc), 0, atol=atol) and np.isclose(np.amax(bc), config.Ly, atol=atol):
        hist, bin_edges = np.histogram(bc, bins=20, density=True)
        if np.any(np.absolute(hist) < 1e-2):
            orientation += np.array([0, 0, 0])
        else:
            orientation += np.array([0, 1, 0])

    if np.isclose(np.amin(cc), 0, atol=atol) and np.isclose(np.amax(cc), config.Lz, atol=atol):
        hist, bin_edges = np.histogram(cc, bins=20, density=True)
        if np.any(np.absolute(hist) < 1e-2):
            orientation += np.array([0, 0, 0])
        else:
            orientation += np.array([0, 0, 1])

    if (orientation == np.array([0, 0, 0])).all():
        shape = 'sphere'
    else:
        shape = 'channel'
        if not np.isclose(np.linalg.norm(orientation), 0, atol=1e-4):
            orientation = orientation/np.linalg.norm(orientation)

    Npoints = len(ac)
    print('--------------- Cluster # %d ---------------' %cluster_number)
    print('Center (a, b, c) = %1.3g, %1.3g, %1.3g' %(center[0], center[1], center[2]))
    print('shape = %s' % shape)
    length = 0
    if shape == 'sphere':
        # based on the number of points within a cluster
        #orientation = np.array([0, 0, 0])
        size = 2 * (3.* Npoints/(4* np.pi * config.rho)) ** (1 / 3.) # diameter
        print('Diameter = %1.2f \n' %size)

    elif shape == 'channel':
        #orientation = calculate_orientation_of_channel(ac, bc, cc)
        da, db, dc = np.multiply(orientation, np.array([config.Lx, config.Ly, config.Lz]))
        length = np.sqrt(da**2 + db**2 + dc**2 + 2*da*db*np.cos(config.gamma)\
                       + 2*db*dc*np.cos(config.alpha) + 2*dc*da*np.cos(config.beta))
        size = (4*Npoints/(np.pi*length*config.rho)) ** (1/2.) # diameter

        print('Orientation (a, b, c) = ', np.round(orientation,2))
        print('Length = %1.2f, Diameter = %1.2f \n' %(length, size))

    return shape, size, length, orientation

def calculate_orientation_of_channel(ac, bc, cc):

    points_xyz = []
    points_abc = []
    for ai, bi, ci in zip(ac, bc, cc):
        points_xyz.append(abc_to_xyz(ai, bi, ci))
        points_abc.append(np.array([ai, bi, ci]))
    points_xyz = np.array(points_xyz)
    points_abc = np.array(points_abc)


    cov = np.cov(points_abc, rowvar =False)
    #cov = np.cov(points_xyz, rowvar=False)

    # Check the eigenvalues of the covariance matrix
    eigvals, eigenvacs = np.linalg.eig(cov)

    #sort the eigenvalues, largest eigenvalues first
    index = np.argsort(eigvals)[::-1]
    eigvals = eigvals[index]
    eigenvacs = eigenvacs[:,index]

    # return first eigenvector
    return eigenvacs[:, 0]



def centralize_cluster(ac, bc, cc, center):
    """
    Shift the coordinates of the cluster so that its center is at the center of the unit cell

    :param ac:
    :param bc:
    :param cc:
    :return:
    """
    an, bn, cn = [], [], []

    for ai, bi, ci in zip(ac, bc, cc):
        ai += (config.Lx/2 - center[0])
        bi += (config.Ly/2 - center[1])
        ci += (config.Lz/2 - center[2])
        ai, bi, ci = put_point_in_box(ai, bi, ci)
        an.append(ai)
        bn.append(bi)
        cn.append(ci)

    # convert to numpy arrays
    an = np.array(an)
    bn = np.array(bn)
    cn = np.array(cn)

    return an, bn, cn

def identify_cluster_shape(bin_number, cluster_number, ac, bc, cc, center):
    """

    :param ac: numpy array of all points of a particle cluster
    :param bc:
    :param cc:
    :param center: numpy array, best cluster center
    :return:
    """
    # Check for 0 or 1 points
    points_abc = []

    # Some clusters are divided into 2, 4 or 8 parts because of the periodic boundary conditions
    # Here I shift the coordinates such that the points are all together
    #
    for ai, bi, ci in zip(ac, bc, cc):
        ai += (config.Lx/2 - center[0])
        bi += (config.Ly/2 - center[1])
        ci += (config.Lz/2 - center[2])
        ai, bi, ci = put_point_in_box(ai, bi, ci)
        points_abc.append([ai, bi, ci])

    if len(points_abc) == 0:
        return "empty"
    elif len(points_abc) == 1:
        return "point"

    # Check for a line
    elif len(points_abc) == 2:
        return "line"

    # Check for a sphere
    else:
        # Compute the mean of the points
        mean = np.mean(points_abc, axis=0)

        # Compute the distances of the points from the mean
        #distances = [np.linalg.norm(p - mean) for p in points]
        distances = [periodic_distance(p, mean) for p in points_abc]

        # Check if the distances are all similar
        #print(distances)
        if np.isclose(distances, distances[0], rtol=0.1).all():
            return "sphere"

        # If the distances are not all similar, check for a cylinder
        else:
            # Compute the covariance matrix of the points
            points_xyz = []
            for points in points_abc:
                a, b, c = points
                points_xyz.append(abc_to_xyz(a, b, c))

            #print(np.shape(points_xyz))
            points_xyz = np.array(points_xyz)
            save_cluster_as_xyz(bin_number, cluster_number, points_xyz[:,0], points_xyz[:,1], points_xyz[:,2])
            # each
            cov = np.cov(points_xyz, rowvar=False)

            # Check the eigenvalues of the covariance matrix
            eigvals, _ = np.linalg.eig(cov)

            #sort the eigenvalues, largest eigenvalues first
            index = np.argsort(eigvals)[::-1]
            eigvals = eigvals[index]

            print(eigvals)
            # If two eigenvalues are much larger than the third, the points lie on a cylinder
            if not np.isclose(eigvals[0], eigvals[1], rtol=1) and not np.isclose(eigvals[0], eigvals[2], rtol=1)\
                    and np.isclose(eigvals[2], eigvals[1], rtol=1):
                return "cylinder"
            else:
                return "3d cluster"