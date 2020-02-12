#-*-Python-*-

plt.ion()
from boutdata import collect # import all module in boutdata
from boututils import save2nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# set up the grid file
g=boutgrid('cbm18_dens4_68x64_couple.nc')
grid='cbm18_dens4_68x64_couple.nc'
path='data-trans'

# collect data from transport run
Te_x=collect('Te_x',path=path)
Ni_x=collect('Ni_x',path=path)
Bbar=collect('bmag',path=path)
tbar=collect('tbar',path=path)
Lbar=collect('Lbar',path=path)
pei=collect('pei',path=path,tind=[-1,-1])
density_unit=1.e20  #1.E20/m^3,
ee = 1.602e-19        # electron charge

pei=np.squeeze(pei)*Ni_x*Te_x
P_mean = pei*density_unit*ee    ## calculate P_mean in unit Pascals

# save pressure to grid file
print 'save_P'
save2nc(grid, 'a', pressure_couple=P_mean)

