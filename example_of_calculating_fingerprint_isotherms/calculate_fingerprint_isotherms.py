#! /usr/bin/env python

import pandas as pd
import numpy as np
import glob, re


# INPUTS
pore_type_matrix_csv = 'Cu-BTC_pore_type_matrix.csv'
raspa_movie_folder = 'RASPA_GCMC_Simulations'
Outfile = 'fingerprint_isotherms.dat'
Temperature = 87.29

# extract out pdb files to read
pdbfiles = glob.glob(f'{raspa_movie_folder}/P_*/Movies/System_0/*all*')

# read pore_type_matrix_csv file
df = pd.read_csv(pore_type_matrix_csv)
# x, y, z, and pore_type_labels are 1d numpy arrays
x, y, z = df.iloc[:, 0].values, df.iloc[:, 1].values, df.iloc[:, 2].values
pore_type_label_list = df.iloc[:, 3].values
matrix_len = int(np.cbrt(len(pore_type_label_list)))

P = []  # Pressure in kPa
# Three fingerprint isotherms for Cu-BTC
fingerprint_0, fingerprint_1 , fingerprint_2 = [], [], []
uncounted_list = []

for pdbfile in pdbfiles:
    # N0, N1, N2 are number of molecules adsorbed in compartments 0, 1, 2, respectively
    N0, N1, N2 = 0, 0, 0
    uncounted = 0
    frame_number = 0

    # extract pressure from pdb file
    pattern = "%3.6f_(.*?)_allcomponents.pdb" % Temperature
    Pa = float(re.search(pattern, pdbfile).group(1))  # in Pascals
    P.append(Pa / 1000.0)  # in kPa

    # reading each adsorbate pdb file
    with open(pdbfile) as inp:
        for line in inp:
            if 'CRYST1' in line:
                words = line.split()
                # Lx, Ly, Lz are the lengths of 1x1x1 unit cell
                Lx, Ly, Lz = float(words[1]) / 2, float(words[2]) / 2, float(words[3]) / 2
            if 'ENDMDL' in line:
                frame_number += 1
                continue
            elif 'ATOM' in line:
                words = line.split()
                Xcom, Ycom, Zcom = float(words[4]), float(words[5]), float(words[6])
                if Xcom > Lx:
                    Xcom -= Lx
                if Ycom > Ly:
                    Ycom -= Ly
                if Zcom > Lz:
                    Zcom -= Lz
                xint, yint, zint = int(np.floor(Xcom)), int(np.floor(Ycom)), int(np.floor(Zcom))

                index = xint * matrix_len ** 2 + yint * matrix_len + zint
                pore_type_label = pore_type_label_list[index]

                if pore_type_label == -1:
                    with open('uncounted.pdb', 'a') as out:
                        out.write(line)
                    uncounted += 1

                if pore_type_label == 0:
                    N0 += 1
                if pore_type_label == 1:
                    N1 += 1
                if pore_type_label == 2:
                    N2 += 1

    fingerprint_0.append(N0 / frame_number)
    fingerprint_1.append(N1 / frame_number)
    fingerprint_2.append(N2 / frame_number)
    uncounted_list.append(uncounted / frame_number)

# save to file
df = pd.DataFrame({'Pressure (kPa)': P, 'N0': fingerprint_0, 'N1': fingerprint_1, 'N2': fingerprint_2})
df = df.sort_values(by='Pressure (kPa)')
df.to_csv(Outfile, sep='\t', index=False, float_format='%15.3f')
print(f'Written {Outfile}')
