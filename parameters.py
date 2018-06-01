from __future__ import division
from past.utils import old_div
import matplotlib
from matplotlib import rc
from matplotlib import axes
from numpy import *
from pylab import *
from scipy.integrate import *
from scipy.interpolate import *
from scipy.optimize import curve_fit
from scipy.optimize import *
import scipy.optimize
import scipy.special

#physical parameters
mu=1.1 #magnetic moment in units 10^30 G cm^3
ps=5.04 # spin period in seconds
mdotglobal=11.1 #accretion rate in units 4 pi G M/ (c\kappa)
#coefficients
alpha=0.1 #viscosity
eta=0.0 #accretion efficiency 
kt=0.5  # Btor/Bpol
epsilon=0.5 #some fraction of the radiation that goes to accelerate the outflows
psi=1. #accounts for the net angular momentum lost in the wind
n=1. #vertical structure parameter, \rho \propto (1-(z/H)^2)^n
mmean=0.6 #mmean=m_particle/m_proton (mmean=0.6 for completely ionized Solar-metallicity matter)
#dimensionless coefficients
lam=3.9848e10
chi=old_div(8.8e-6,mmean)
pstar=4.33e-5
#vertical structure
tvert=old_div(219.,1024.) #vertical structure n=1
hvert=sqrt(5.)  #vertical structure n=1
#program coefficients
tol=1e-6
toff=1e-12
varold=0. #old temperature profile, n=1
varnew=1. #new temperature profile with any n
qeqest=0.905 # initial estimate for qeq 
xiest=1. # initial estimate for xi
defac=0.9 # grid non-linearity parameter
routfactor=100. # outer, starting, radius in inner radius units
htorcr=1000.0
