#-*-Python-*-

plt.ion()
from boutdata import collect # import all module in boutdata
from boututils import deriv
from boututils import save2nc
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# set constant parameters
g=boutgrid('cbm18_dens4_68x64_couple.nc')
grid='cbm18_dens4_68x64_couple.nc'
ee=1.602e-19    # in SI unit
Zi=1.0
Nbar=1.0e20     # density unit in /m^3
Tibar=10.            # in unit eV
path='data-elm'

# collect data from the 3 field run
Div_Flux_couple=collect('Div_Flux_couple',path=path,tind=[-2,-1])
Div_Flux_couple=np.squeeze(Div_Flux_couple)
dcflux=Div_Flux_couple.mean(axis=2)
#dcflux=Div_Flux_couple
#plt.ion()
#plt.plot(Div_Flux_couple[:,32,-1])
#plt.plot(dcflux[:,32])

#save flux to grid file
print 'save_flux'
save2nc(grid, 'a', Div_Flux_couple=dcflux)

