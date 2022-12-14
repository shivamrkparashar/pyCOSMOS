import config
from math import pi

def epsilon(rho):
    return (1/(rho))**(1/3.)


def number_of_points_in_cluster(cluster_diameter, rho):
    Nmax = 50000
    rho = Nmax/config.Box**3
    return rho*pi/6*cluster_diameter**3


def calculate_dbscan_parameters(avg_geom_dia):
    Nmax = 50000
    rho = Nmax/config.Volume_of_uc

    epsilon = (1/(rho))**(1/3.)
    number_of_points_in_cluster = rho*pi/6*avg_geom_dia**3

    min_points = 0.05*number_of_points_in_cluster
    return epsilon*3, min_points


Nmax = 50000
Box = 28.9058
Vuc = Box**3

rho = Nmax/Vuc

for di in [5, 10, 12, 15, 20, 25]:
    print(di, 0.1*number_of_points_in_cluster(di, rho))
