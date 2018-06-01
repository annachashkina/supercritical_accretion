from __future__ import print_function
from __future__ import division
from builtins import str
from past.utils import old_div
import matplotlib
from matplotlib import rc
from matplotlib import axes
from numpy import *
from pylab import *
from scipy.integrate import *
from scipy.interpolate import *
from scipy.optimize import curve_fit

#Uncomment the following if you want to use LaTeX in figures
rc('font',**{'family':'serif','serif':['Times']})
rc('mathtext',fontset='cm')
rc('mathtext',rm='stix')
rc('text', usetex=True)
# #add amsmath to the preamble
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amssymb,amsmath}"]

import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import numpy.random
import time
import os
from physics import G
from parameters import *

bound=1 #1 - local pressure balance, 0 - equatorial pressure balance

def parset(**kwargs):
#(newmu=newmu, neweta=neweta, newp=newp, newmdot=newmdot, newalpha=newalpha):
    global mu, eta, ps, mdotglobal, alpha

    mu=kwargs['newmu']
    eta=kwargs['neweta']
    ps=kwargs['newps']
    mdotglobal=kwargs['newmdot']
    alpha=kwargs['newalpha']

    if(ps<0.):
        ps=-ps*peq()

    print("bocon parset:")
    print("  now mu = "+str(mu))
    print("  now mdot = "+str(mdotglobal))
    print("  now eta = "+str(eta))
    print("  now P = "+str(ps)+"s")
    print("  now alpha = "+str(alpha))
   
def peq():
    return pstar*(lam*mu**2/mdotglobal)**(old_div(3.,7.))*2.**(old_div(-3.,14.))
   
def rafun():
    return (lam*mu**2/mdotglobal)**(old_div(2.,7.))*2.**(old_div(-1.,7.))
   
def oin(rin, hin,mdotin):
    beta=old_div(hin,rin)
#    print "bocon.oin: mu = "+str(mu)
#    print "beta= "+str(beta)+'\n'
#    print "rin= "+str(rin)+'\n'
    if((beta*eta)>=0.9):
        beta=old_div(0.9,eta)
    return rin**1.5/(1.-eta*beta)*(old_div(pstar,ps)+2.*kt*lam*mu**2.*hin/rin**6./mdotin)

def fhin(rin, tau,mdotin):
# aprobado: formula especial para espesor vertical del disco interior
    return 2.*hvert**2.*rin*(eta*mdotin+lam*mu**2./rin**4.)/tau

def fwrfin(rin, hin,mdotin):
   coef=G(n+1)
   if(bound==1):
      return 2.*G(n+1)*alpha*(eta*mdotin+lam*mu**2./rin**4.)/rin**2.*hin
   else:
      return 2.*alpha*(eta*mdotin+lam*mu**2./rin**4.)/rin**2.*hin

def quartic(c1,c2):
    # solves x^4+c2x+c1=0 for c1<0, c2>0, when there is !1 solution
    p=-4./3.*c1
    q=-c2**2
    ap=old_div((sqrt(q**2+4.*p**3)-q),2.)
    am=old_div(-(sqrt(q**2+4.*p**3)+q),2.)
    t=ap**(old_div(1.,3.))+sign(am)*(abs(am))**(old_div(1.,3.))

    if(t<=toff):
        return sqrt(sqrt(-c1))
    x=sqrt(t)*(sqrt(2.*c2/t**1.5-1.)-1.)/2.
    if(x<=toff):
        return old_div(-c1,c2)
    else:
        return x

def ctemp(h, wrf, tau):
    c1=-45./256./pi/alpha * wrf/h
    c2=9./64./pi*chi * tau/h
    t=quartic(c1,c2)
    return t


def ABCfun(wrf, rin, hin, tau, omega,mdot):
    omegain=omega
    wrfin=wrf
    r=rin
    h=hin
    a=2.*wrf*r 
    b=16.*pi*mdot*r*ctemp(h,wrf,tau)**4/wrf/(1.+tvert*tau)
    c=4.*mdot/r**0.5
    k=16.*pi*r**2.5*ctemp(h,wrf,tau)**4/(1.+tvert*tau)*omega*epsilon*(1.-psi)
    return -a-b+c+k+alpha*tau*(1.-omega**2)

def tausolve(rrin,mdotin):
    ra=rafun()
    rin=ra*rrin
   
    tau1=1.
    tau2=1.e10
    
    hin1=fhin(rin,tau1,mdotin)
    omega1=oin(rin,hin1,mdotin)
    wrf1=fwrfin(rin, hin1,mdotin)
    hin2=fhin(rin,tau2,mdotin)
    omega2=oin(rin,hin2,mdotin)
    wrf2=fwrfin(rin, hin2,mdotin)

    f1=ABCfun(wrf1, rin, hin1, tau1, omega1,mdotin)
    f2=ABCfun(wrf2, rin, hin2, tau2, omega2,mdotin)
    if((f1*f2)>=0.):
        tau1=taumin(rrin,mdotin) 
        hin1=fhin(rin,tau1,mdotin)
        omega1=oin(rin,hin1,mdotin)
        wrf1=fwrfin(rin, hin1,mdotin)
        f1=ABCfun(wrf1, rin, hin1, tau1, omega1,mdotin)
        if((f1*f2)>=0.):
            return sqrt(-1.),sqrt(-1.),sqrt(-1.),sqrt(-1.)
    while(abs(old_div((tau1-tau2),(tau1+tau2)))>tol):
        tau=sqrt(tau1*tau2)
        hin=fhin(rin,tau,mdotin)
        omega=oin(rin,hin,mdotin)
        wrf=fwrfin(rin, hin,mdotin)
        f=ABCfun(wrf, rin, hin, tau, omega,mdotin)
        
        if((f1*f)<0.):
            tau2=tau
        else:
            tau1=tau

    return tau, hin, omega,wrf

def taumin(rrin,mdotin):
    ra=rafun()
    rin=ra*rrin

    tau1=1.
    tau2=1.e10
    
    hin1=fhin(rin,tau1,mdotin)
    omega1=oin(rin,hin1,mdotin)
    wrf1=fwrfin(rin, hin1,mdotin)
    hin2=fhin(rin,tau2,mdotin)
    omega2=oin(rin,hin2,mdotin)
    wrf2=fwrfin(rin, hin2,mdotin)
    f1=ABCfun(wrf1, rin, hin1, tau1, omega1,mdotin)
    f2=ABCfun(wrf2, rin, hin2, tau2, omega2,mdotin)

    while(abs(old_div((tau1-tau2),(tau1+tau2)))>tol):
        tau=sqrt(tau1*tau2)
        hin=fhin(rin,tau,mdotin)
        omega=oin(rin,hin,mdotin)
        wrf=fwrfin(rin, hin,mdotin)
        f=ABCfun(wrf, rin, hin, tau, omega,mdotin)

        if(f<0.):
            return tau
        else:
            if(f1<f2):
                tau2=tau
                f2=f
            else:
                tau1=tau
                f1=f
    return tau


