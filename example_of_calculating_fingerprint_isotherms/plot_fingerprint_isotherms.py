#! /usr/bin/env python

from pylab import *
import pandas as pd

fs = 15
rcParams['font.size'] = '%s' %fs
P0 = 0.974*101.325


def Column(filename, *args):
    """
    Reads a column formatted file separated by spaces
    """
    df = pd.read_csv(filename, delim_whitespace=True)
    return [df.iloc[:, a] for a in args]


# Adsorption
P, N0, N1, N2 = Column('fingerprint_isotherms.dat', 0, 1, 2, 3)
#P /= P0
semilogx(P, N0, '--', label = 'L3', color = 'orange')
semilogx(P, N1, '--', label = 'L2', color = 'green')
semilogx(P, N2, '--', label = 'S1', color = 'magenta')
semilogx(P, N0+N1+N2, '-ko', label = 'total', markerfacecolor = 'black')


title('Ar adsorption in Dehy Cu-BTC at 87.3 K')
xlabel('Fugacity (kPa)')
#xlabel('$P/P_0$')
ylabel('Number of molecules per uc')
legend(loc = 'best')
savefig('fig_fingeprints_with_total_fugacity.jpg', format = 'jpg', dpi=400 ,bbox_inches='tight')
show()

