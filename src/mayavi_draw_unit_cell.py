
import mayavi.mlab as mlab
import numpy as np

black = (0,0,0)
white = (1,1,1)
mlab.figure(bgcolor=white)

def draw_line(A, B):
    """


    radius of the tubes used to represent the lines, If None, simple lines are used
    :param A:
    :param B:
    :return:
    """
    mlab.plot3d([A[0],B[0]], [A[1], B[1]], [A[2], B[2]], color = black, tube_radius = None)

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
    return mlab.mesh(r * x + x0, r * y + y0, r * z + z0, color = color)

plot_sphere(Lx, 0, 0, 2.4, (1, 0, 0))

#mlab.plot3d([0, 1000], [0, 0], [0, 0], color=black, tube_radius=10.)
#mlab.plot3d([0, 0], [0, 1500], [0, 0], color=black, tube_radius=10.)
#mlab.plot3d([0, 0], [0, 0], [0, 1500], color=black, tube_radius=10.)


mlab.show()
