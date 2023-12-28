import pandas as pd
from pylab import *
import numpy as np

import config

fs = 15
rcParams['font.size'] = '%s' %fs


def plot_psd_bar_pore_type(bin_centers, pore_type_histogram):

    if config.pore_type_count == 1:
       d = {'bin_centers': bin_centers, 'pore_type_hist_1': pore_type_histogram[0]}

    elif config.pore_type_count == 2:
        d = {'bin_centers': bin_centers, 'pore_type_hist_1': pore_type_histogram[0],
         'pore_type_hist_2':pore_type_histogram[1]}

    elif config.pore_type_count == 3:
        d = {'bin_centers': bin_centers, 'pore_type_hist_1': pore_type_histogram[0],
             'pore_type_hist_2':pore_type_histogram[1], 'pore_type_hist_3': pore_type_histogram[2]}

    elif config.pore_type_count == 4:
        d = {'bin_centers': bin_centers, 'pore_type_hist_1': pore_type_histogram[0],
             'pore_type_hist_2':pore_type_histogram[1], 'pore_type_hist_3': pore_type_histogram[2],
             'pore_type_hist_4': pore_type_histogram[3]}

    df = pd.DataFrame(data =d)
    df.to_csv('pore_type_histogram.dat', index = False, float_format = '%1.3f')

    width_bar = bin_centers[1] - bin_centers[0]
    norm = np.sum(pore_type_histogram)*width_bar # calculates sum over all element across dimensions
    pore_type_histogram /= norm

    colors = config.color_list 

    figure()
    for i in range(config.pore_type_count):
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
        if i == 3:
            bar(bin_centers, pore_type_histogram[i], width=width_bar,
                bottom=pore_type_histogram[0] + pore_type_histogram[1] + pore_type_histogram[2],
                label='pore %d' % (i + 1), color=colors[i], edgecolor='black')

    legend(loc = 'best')
    xlabel('Diameter $(\AA)$')
    xlim(0,)
    ylim(0, )
    ylabel('Normalized Distribution')
    savefig('ptd_psd.jpg' ,format='jpg',dpi=300,bbox_inches='tight')


#show()
