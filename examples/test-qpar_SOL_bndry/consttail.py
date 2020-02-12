""" generate const tail Te profile. """

__author__ = "J.G.Chen"
__email__ =  "cjgls@pku.edu.cn"
__date__ = "10/13/2018"
__version__ = "0.0.1"

import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as nd

plt.clf()
ny = 128
dy = 1.
constny = 55
sigma = 2

yind = np.arange(ny) * dy   # [cm]
ymean = ny/2    # np.mean(yind)
te = 1.0 + 10 * np.exp(-((yind - ymean) / 10.0) ** 2)  # [eV]
plt.plot(te, label='Gaussian')

# set const tail
te[:constny] = te[constny]
te[ny-constny:] = te[ny-constny]
assert te[constny] == te[ny-constny]
print "Tail: ", te[constny], te[ny-constny]
plt.plot(te, label='Const Tail')

# smooth
te = nd.gaussian_filter(te, sigma=sigma)
plt.plot(te, label='Smooth')

gradte = np.gradient(te)
plt.plot(gradte, label='dte/dy')
plt.plot(gradte == 0, '*', label='dte/dy == 0')
plt.legend(fontsize=22, title='te')
plt.xlabel('Y')
plt.show()

# write to netcdf file
# save2nc('qparwcoll.grd.nc', 'a', te=np.tile(te, (5, 1)))
