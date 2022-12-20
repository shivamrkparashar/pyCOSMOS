# pyCOSMOS: python for Compartmentalization Of Solid Metal-Organic framework Structures
This algorithm divides the unit cell of a MOF into pore compartments


## How to use?
Build the cython file:  
`python3 setup_periodic_distance.py build_ext --inplace`


## Required: Perform Zeo++ pore size distribution calculation
1. Install Zeo++ from: http://zeoplusplus.org/. After successful installation, "network" executable is generated
2. Perform psd calculation using: `./network -ha -vpsd 1.657 1.657 50000 Structure.cif`. This will generate *vpsdpts files.


## Inputs
1. *vpspts file from Zeo++ pore size distribution calculation
2. Number of pore types. (If you don't know this, run this code with guess value of 1, 2, and 3 in this order)

## Outputs
1. pore_type_matrix_with_cluster_labels.csv
2. pore_type_matrix_with_pore_type_labels.csv
3. pore_type_matrix_with_cluster_center_labels.html
4. pore_type_matrix_with_pore_type_labels.html
5. geometric_points_with_cluster_labels_for_pore_type_*Npores.html 

