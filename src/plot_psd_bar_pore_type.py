import pandas as pd
from pylab import *
import numpy as np

import config

fs = 12
rcParams['font.size'] = '%s' %fs


pore_type_histogram = np.array([[0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
  2.0000e+00, 2.4000e+02, 1.2100e+02, 8.0000e+01, 8.9000e+01, 2.3500e+02,
  1.3100e+02, 1.4900e+02, 1.1200e+02, 1.1100e+02, 8.2000e+01, 7.5500e+02,
  5.6500e+02, 1.1235e+04],
 [0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 6.0000e+00, 7.7000e+01,
  4.9000e+01, 3.8900e+02, 4.4000e+01, 3.4000e+01, 3.2000e+01, 2.0800e+02,
  2.6000e+01, 8.4030e+03, 8.3600e+02, 5.7100e+02, 4.3700e+02, 7.5500e+02,
  7.3900e+02, 9.5400e+02]])


def plot_psd_bar_pore_type(bin_centers, pore_type_histogram):

    if config.Npores == 2:
        d = {'bin_centers': bin_centers, 'pore_type_hist_1': pore_type_histogram[0],
         'pore_type_hist_2':pore_type_histogram[1]}
    elif config.Npores == 3:
        d = {'bin_centers': bin_centers, 'pore_type_hist_1': pore_type_histogram[0],
             'pore_type_hist_2':pore_type_histogram[1], 'pore_type_hist_3': pore_type_histogram[2]}
    df = pd.DataFrame(data =d)
    df.to_csv('pore_type_histogram.dat', index = False)

    width_bar = bin_centers[1] - bin_centers[0]
    norm = np.sum(pore_type_histogram)*width_bar # calculates sum over all element across dimensions
    pore_type_histogram /= norm

    #colors = ['teal', 'gold', 'skyblue']
    #colors = ['green', 'yellow', 'skyblue'] # PCN-224
    #colors = ['yellow', 'green', 'orange'] # ZIF-412
    #colors = ['magenta',  'green', 'orange' ] # Cu-BTC
    colors = ['orange', 'green', 'skyblue'] # UiO-66

    figure()
    for i in range(config.Npores):
        if i == 0:
            bar(bin_centers, pore_type_histogram[i], width = width_bar,
                label = 'pore %d' %(i+1), color = colors[i], edgecolor = 'black')
        if i == 1:
            bar(bin_centers, pore_type_histogram[i], width = width_bar,
                bottom= pore_type_histogram[0],
                label = 'pore %d' %(i+1), color = colors[i], edgecolor = 'black')
        if i == 2:
            bar(bin_centers, pore_type_histogram[i], width = width_bar,
                bottom= pore_type_histogram[0] + pore_type_histogram[1],
                label = 'pore %d' %(i+1), color = colors[i], edgecolor = 'black')

    legend(loc = 'best')
    xlabel('Diameter $(\AA)$')
    xlim(0,)
    ylim(0, )
    ylabel('Normalized Distribution')
    savefig('ptd_psd.jpg' ,format='jpg',dpi=300,bbox_inches='tight')

#show()
