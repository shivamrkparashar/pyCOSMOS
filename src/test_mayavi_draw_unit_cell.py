import mayavi.mlab as mlab
import numpy as np
from matplotlib import colors
import mayavi
from mayavi.api import Engine

"""
This is a testing script
"""


black = (0,0,0)
white = (1,1,1)
mlab.figure(bgcolor=white, size=(1024, 1024))

def draw_line(A, B):
    """
    radius of the tubes used to represent the lines, If None, simple lines are used
    :param A:
    :param B:
    :return:
    """
    mlab.plot3d([A[0],B[0]], [A[1], B[1]], [A[2], B[2]], color = black, tube_radius = None)

def draw_unit_cell(): 
    Lx, Ly, Lz = 10, 10, 10
    alpha, beta, gamma = 90, 90, 90

    vx = [0.00, Lx]
    vy = [0.00, Ly]
    vz = [0.00, Lz]

    cord = []
    for x in vx:
        for y in vy:
            for z in vz:
                cord.append([x, y, z])

    draw_line(cord[0], cord[1])
    draw_line(cord[0], cord[2])
    draw_line(cord[0], cord[4])
    draw_line(cord[5], cord[4])
    draw_line(cord[5], cord[1])
    draw_line(cord[5], cord[7])
    draw_line(cord[3], cord[1])
    draw_line(cord[3], cord[2])
    draw_line(cord[3], cord[7])
    draw_line(cord[6], cord[4])
    draw_line(cord[6], cord[7])
    draw_line(cord[6], cord[2])


def plot_sphere(x0, y0, z0, r, color):
    """
    Draws a sphere with (x0, y0, z0) as center and radius = r
    :param x0: x coordinate of the sphere center
    :param y0: y coordinate of the sphere center
    :param z0: z coordinate of the sphere center
    :param r: radius of the sphere
    :param color: color in rgb format. For example (1,0,0) is red
    :return: mayavi mesh
    """
    # phi, theta = (0, 2pi), (0, pi)
    [phi, theta] = np.mgrid[0:2*np.pi:30j, 0:np.pi:30j]
    x = np.cos(phi) * np.sin(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(theta)
    return mlab.mesh(r * x + x0, r * y + y0, r * z + z0, color=color)

def read_pdb(pdbfile):

    name, X, Y, Z = [], [], [], []
    with open(pdbfile, 'r') as inp:
        for line in inp:
            if 'HETATM' in line:
                s = line.split()
                n, x, y, z = s[-1], float(s[5]), float(s[6]), float(s[7])
                name.append(n)
                X.append(x)
                Y.append(y)
                Z.append(z)

    return name, X, Y, Z


def visualize_pdb(pdbfile):
    """
    visualize pdb file in mayavi
    draws atoms as sphere
    bonds are not plotted
    :param pdbfile: path to pdb file
    :return: None
    """
    name, X, Y, Z = read_pdb(pdbfile)
    radius = {'C': 1.2, 'H': 0.8, 'N': 1.2, 'O': 1.2, 'Zn': 1.7, 'Zr': 1.7}
    color = {'C': 'silver', 'H': 'snow', 'N': 'royalblue', 'O': 'red', 'Zn': 'springgreen', 'Zr':
             'magenta'}
    color_rgb = {k: colors.to_rgb(v) for k, v in color.items()}

    for ni, x, y, z in zip(name, X, Y, Z):
        plot_sphere(x, y, z, radius[ni], color_rgb[ni])

    # set angle and save images
    scene = mlab.gcf()
    scene.scene.parallel_projection = True
    scene.scene.x_minus_view()
    scene.scene.save('img_x_minus.png')

    scene.scene.isometric_view()
    scene.scene.save('img_isometric.png')

    #scene.Savesc

#draw_unit_cell()

dir1 = '/home/shivam/Desktop/Automatic_compartmentalization/Program/2_testing_NporeType_input/BAWROI_freeONLY'
visualize_pdb(dir1 + '/BAWROI_freeONLY.pdb')
#plot_sphere(Lx, 0, 0, 2.4, (1, 0, 0))

#mlab.plot3d([0, 1000], [0, 0], [0, 0], color=black, tube_radius=10.)
#mlab.plot3d([0, 0], [0, 1500], [0, 0], color=black, tube_radius=10.)
#mlab.plot3d([0, 0], [0, 0], [0, 1500], color=black, tube_radius=10.)


mlab.show()
