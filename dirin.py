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
import bocon as b
import physics as f
# import ssdisk as ss
# import readkit as rk
# import readkit as rk

from parameters import *
# functions:
from physics import G, beta, comegaga, ctau, cwrf, ctemp, fh, ftau
from physics import Scal, Pcal, Qcal # SPQR
from physics import xiinf, raffugen, peqgen
# constants
#from physics import mu, ps, mdotglobal, alpha, eta, kt, epsilon, psi, n, mmean, lam, chi, tvert, hvert
# more constants:
#from physics import tol, toff, varold, varnew, qeqest, xiest, defac
# visualization functions:
# from plots import

def XQset(newxi, newqeq):
    global xiest, qeqest
    print "XQ set"
    xiest=newxi
    qeqest=newqeq

def parset(**kwargs):
#(newmu=newmu, neweta=neweta, newp=newp, newmdot=newmdot, newalpha=newalpha):
    global mu, eta, ps, mdotglobal, alpha

    mu=kwargs['newmu']
    eta=kwargs['neweta']
    ps=kwargs['newps']
    mdotglobal=kwargs['newmdot']
    alpha=kwargs['newalpha']

    if(ps<0.):
        ps=-ps*b.peq()

    print "dirin parset:"
    print "  now mu = "+str(mu)
    print "  now mdot = "+str(mdotglobal)
    print "  now eta = "+str(eta)
    print "  now P = "+str(ps)+"s"
    print "  now alpha = "+str(alpha)

# dOmega/dr = 
def domega(omega,r,tau,mdot,wrf,tc):
    return -0.5*omega/r+alpha*tau*(omega**2.-1.)/(2.*mdot*r**0.5)+64.*pi/3.*G(n)*(n+1.)*epsilon*(psi-1.)*omega*tc**4.*r**2./tau/mdot+wrf*r**0.5/mdot

# dmdot/dr = 
def dmdot(r,tc,tau):
    return 64.*pi/3.*G(n)*(n+1.)*epsilon*r**2*tc**4./tau

# dwrf/dr = 
def dwrf(tau,omega,r):
    return alpha*tau*(omega**2.-1)/r**2.

# dtau/dr = 
def dtau(tau,tc,wrf,r,omega,mdot):
    bet=beta(tau,tc,wrf)
    csig=ctau(tau,tc,wrf,r,mdot)
    ctem=ctemp(tau,tc,wrf,r,mdot)
    cw=cwrf(tau,tc,wrf,r,mdot)

    part1=csig+ctem/8.*(1.-3.*bet)/(1.-0.75*bet)
    part2=comegaga(tc,tau,r,wrf)-1.5*omega/r**2.5
    part3=r**(-1.5)*domega(omega,r,tau,mdot,wrf,tc)
    part4=alpha*tau/(wrf*r*r)*(omega*omega-1.)*(cw+ctem/8.*(1.-3.*bet)/(1.-0.75*bet))
    part5=3.*ctem/(8.*r)*(1.-bet)/(1.-0.75*bet)-3.*cw/r
    t=tau/part1*(part2+part3-part4+part5)

    return t

# dTc/dr = 
def dtemp(tau,tc,wrf,r,omega,mdot):
    bet=beta(tau,tc,wrf)
    part1=(8.-6.*bet)
    part2=(1.-3.*bet)/tau*dtau(tau,tc,wrf,r,omega,mdot)
    part3=alpha*tau*(omega*omega-1.)*(1.+bet)/wrf/r/r
    part4=3.*(1.-bet)/r
    return tc/part1*(part2+part3-part4)

# main disc integration outside in
def rastr(rrin, qeq):
    # rrin = xi, the inner disc radius in rA units, qeq
    ra=b.rafun()
    rin=ra*rrin
    t=1
    rout=100.*rin
    ddr=-1.e-4
    mdot=mdotglobal
    
    if(mdot<1.5):
        ddr*=(mdot/1.5)
    ddr/=(1.+0.5*(b.peq()/ps))
#    defac=0.99

    r=rout
    omega=1.
    oprev=omega
    rprev=r

    htormax=0.
    wrf=2.*mdot/r**2*omega*(sqrt(r)-qeq*sqrt(rin))
    wprev=wrf
#    tc=0.00486297109754*tcqeq
    tau=4./chi**0.8*(pi/9.)**0.2*mdot**0.6/alpha**0.8/rout**0.6*(1.-qeq*(rin/rout)**0.5)**0.6
 #   tau=ftau(wrf,r,tc)
    h=hvert*sqrt(r**3*wrf/alpha/tau)
    tc=b.ctemp(h,wrf,tau)
    
    hprev=h
    tauprev=tau
    tprev=tc

    while(r>=rin):
        dr=ddr*(r-rin*defac)
        r1=r+dr/2.
        h=hvert*sqrt(r**3*wrf/alpha/tau)
        
        if(isnan(tau)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.)
        omega1=omega+domega(omega,r,tau,mdot,wrf,tc)*dr/2.
        wrf1=wrf+dwrf(tau,omega,r)*dr/2.
        tau1=tau+dtau(tau,tc,wrf,r,omega,mdot)*dr/2.
        tc1=tc+dtemp(tau,tc,wrf,r,omega,mdot)*dr/2.
        beta1=beta(tau1,tc1,wrf1)
        
        if((r*r)>(9.*tau/(4.*64.*pi*tc**4.))):
            mdot1=mdot+dmdot(r,tc,tau)*dr/2.
        else:
            mdot1=mdot
 
        if(wrf1<=0.):
            print "negative stress! wrf = "+str(wrf1)
            return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.),sqrt(-1.), sqrt(-1.)  

     
#        h1=fh(wrf1,r1,tc1)
        h1=hvert*sqrt(r1**3*wrf1/alpha/tau1)
        
        if(isnan(tau1)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.)
        oprev=omega
        omega=omega+domega(omega1,r1,tau1,mdot1,wrf1,tc1)*dr
        wprev=wrf
        wrf=wrf+dwrf(tau1,omega1,r1)*dr
        tcprev=tc
        tc=tc+dtemp(tau1,tc1,wrf1,r1,omega1,mdot1)*dr
 #       tau=ftau(wrf,r,tc)
        tauprev=tau
        tau=tau+dtau(tau1,tc1,wrf1,r1,omega1,mdot1)*dr

        hprev=h
##        h=fh(wrf,r,tc)
        if((r1*r1)>(9.*tau1/(4.*64.*pi*tc1**4.))):
            mdot=mdot+dmdot(r1,tc1,tau1)*dr
        rprev=r
        r+=dr

        if((h/r)>htormax):
            htormax=h/r
            rmax=r
    mdotin=mdot
    wrfin=wrf

    if(isnan(omega)|(r>rin)):
        return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.),sqrt(-1.), sqrt(-1.) 
    oint=(omega-oprev)*(rin-rprev)/(r-rprev)+oprev
    hint=(h-hprev)*(rin-rprev)/(r-rprev)+hprev
    tint=(tau-tauprev)*(rin-rprev)/(r-rprev)+tauprev
    wint=(wrf-wprev)*(rin-rprev)/(r-rprev)+wprev

    omega_in=b.oin(rin, hint,mdotin)
    wrf_in=b.fwrfin(rin, hint,mdotin)
#    print "rastr: omega = "+str(omega)+"; h = "+str(h)+"; htormax = "+str(htormax)+"; mdotin = "+str(mdot)+"; wrf = "+str(wrf)+"\n"
#    print "rastr: omega_in = "+str(omega_in)+"; h_in = "+str(hint)+" ; wrf_in = "+str(wrf)+"\n"

    return oint, hint, tint, htormax, mdotin, wint
# criterion for the boundary conditions -- do we fit them?

def doffwrfin(xi,qeq):
    ra=b.rafun()
    rin=ra*xi
    oin,hin,tauin,hrmax,mdotin,wrfin=rastr(xi,qeq)
    omegaBC=b.oin(rin, hin,mdotin)
    wrfBC=b.fwrfin(rin, hin,mdotin)
    print "oin-omegaBC "+str(oin-omegaBC)+" wrfin-wrfBC "+str(wrfin-wrfBC)+'\n' 
 #   tauBC, hBC, omegaBC, wrfBC = b.tausolve(xi, mdotin)
 #   print 'here2'
    return oin-omegaBC, wrfin-wrfBC
# wrapper for root calculation
def vrapper(arg):
    xi=arg[0]
    qeq=arg[1]
    a,b=doffwrfin(xi,qeq)
#    b=doffo(xi,qeq,tc)
    print "vrap a= "+str(a)+" vrap b= "+str(b)+'\n'
    print "xi= "+str(xi)+" qeq= "+str(qeq)+'\n'

    return (a,b)
# main procedure searching the root

def ordiv_smart(newmu, newmdot, newps):
    b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.0,newalpha=0.1)
    parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.0,newalpha=0.1)
    #    print 'here'
    print "ordiv_smart mu = "+str(mu)
    print "ordiv_smart: xiest = "+str(xiest)+", qeqest = "+str(qeqest)
    tstart=time.time()
    co=scipy.optimize.root(vrapper,(xiest,qeqest),method='hybr',jac=None,tol=1e-4,callback=None,options=None)
    tend=time.time()
    print "co.x[0]= "+str(co.x[0])+' \n'
    print "co.x[1]= "+str(co.x[1])+' \n'
    print "calculation took "+str(tend-tstart)+"s = "+str((tend-tstart)/60.)+"min"
    if(not(co.success)):
        print "ordiv_smart not coverged"
    return co.x[0], co.x[1], co.success

def corot(newmu, newmdot, newps,rrin,qeq):
    b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
    parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
    pstar=4.33e-5

    oint, hint, tint, htormax, mdotin, wint=rastr(rrin, qeq)
    ral=rrin*(lam*mu**2/mdotin)**(2./7.)*2.**(-1./7.)
    rco=(newps/pstar)**(2./3.)

    print ral/rco

###############################
