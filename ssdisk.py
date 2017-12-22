#THIS REPRODUCE STANDARD DISK
import matplotlib
from matplotlib import rc
from matplotlib import axes
from numpy import *
from pylab import *
from scipy.integrate import *
from scipy.interpolate import *
from scipy.optimize import curve_fit

from parameters import all

#Uncomment the following if you want to use LaTeX in figures 
rc('font',**{'family':'serif','serif':['Times']})
rc('mathtext',fontset='cm')
rc('mathtext',rm='stix')
rc('text', usetex=True)
# #add amsmath to the preamble
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amssymb,amsmath}"] 

import numpy as np
import matplotlib.pyplot as plt
import numpy.random
import time
import os

#tau in A zone
def ftaussa(r, dotm,alpha, rin): #r is the radius in GM/c^2
    fcorr=1.-sqrt(rin/r)
    fcorr=(fcorr+fabs(fcorr))/2.
    taua=4.07/alpha/dotm*r**1.5*0.34/fcorr
#    print taua
# 4.6/alpha/dotmss*rss**1.5/2.5
    return taua

#tau in B zone
def ftaussb(r, mdot,alpha, rin):
#    dotmss=mdot*0.07
#    rss=r/6.
    fcorr=1.-sqrt(rin/r)
    fcorr=(fcorr+fabs(fcorr))/2.
    taub=93200.*alpha**(-0.8)*mdot**0.6*r**(-0.6)*0.35*fcorr**0.6
#1.7e5*alpha**(-0.8)*dotmss**0.6*rss**(-0.6)/2.5
    return taub


#H/R in A zone
def hssa(mdot,r, rin):
    fcorr=1.-sqrt(rin/r)
    fcorr=(fcorr+fabs(fcorr))/2.
    return 1.5*mdot*fcorr*sqrt(5.)

#H/R in B zone
def hssb(r, mdot,alpha, rin):
    fcorr=1.-sqrt(rin/r)
    fcorr=(fcorr+fabs(fcorr))/2.
    h=0.00747*alpha**(-0.1)*mdot**0.2*r**(1./20.)*fcorr**0.2
 #   print h
#    jji=raw_input()
    return 0.00747*alpha**(-0.1)*mdot**0.2*r**(21./20.)*fcorr**0.2*sqrt(5.)

#non-Keplerianity in Standard disk
def ost(r, rin, mdot, alpha, tau):
    return 0.5*mdot/alpha/tau/r*(3.*sqrt(r)-5.*sqrt(rin))

