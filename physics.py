import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import time
import os
import scipy.optimize
import scipy.special

#physical parameters
mu=1.
ps=10.
mdotglobal=100.
#coefficients
alpha=0.1
eta=0.0
kt=0.5
epsilon=0.5
psi=1.
n=1.
mmean=0.5
#dimensionless coefficients
lam=3.9848e10
chi=8.8e-6/mmean
pstar=4.33e-5
#vertical structure
tvert=219./1024.
hvert=sqrt(5.)
#program coefficients
tol=1e-7
toff=1e-12
varold=0. #old temperature profile
varnew=1. #new temperature profile
qeqest=0.905
xiest=0.63
defac=0.99

# vertical integration factor, int(1-x^2)^n dx
def G(n):
    return sqrt(pi)*scipy.special.gamma(n+1.)/scipy.special.gamma(n+3./2.)

# no one knows what for...
def par0(xs, ys):
    n=ys[0]*(xs[1]**2-xs[2]**2)+ys[1]*(xs[2]**2-xs[0]**2)+ys[2]*(xs[0]**2-xs[1]**2)
    d=ys[0]*(xs[1]-xs[2])+ys[1]*(xs[2]-xs[0])+ys[2]*(xs[0]-xs[1])
    return n/2./ d

# gas to total pressure
def beta(tau,tc,wrf):
 #   print "beta= "+str(8.*chi/5.*alpha*tau*tc/wrf)
#    jjii=raw_input()
    return chi*alpha*tau*tc/wrf*2.*(n+1.)/(2.*n+3.)

def comegaga(tc,tau,r,wrf):
    return 64.*pi/3.*G(n)*(n+1.)*tc**4./(tau*r*wrf)

def ctau(tau,tc,wrf,r,mdot):
    bet=beta(tau,tc,wrf)
    return 2./(G(n+1.))*mdot/(alpha*r*r*tau)*Scal(bet)

def cwrf(tau,tc,wrf,r,mdot):
    bet=beta(tau,tc,wrf)
    return 2./(G(n+1.))*mdot/(alpha*r*r*tau)*Pcal(bet)

def ctemp(tau,tc,wrf,r,mdot):
    bet=beta(tau,tc,wrf)
    return 2./(G(n+1.))*mdot/(alpha*r*r*tau)*Qcal(bet)

def Scal(beta):
    part1=sqrt(pi)*(5.*n-3.)*beta/16.*scipy.special.gamma((5.*n+1.)/4.)/scipy.special.gamma((5.*n+7.)/4.)
    part2=sqrt(pi)*(3.-n)*(1.-beta)/2.*scipy.special.gamma(n+1.)/scipy.special.gamma(n+2.5)
    part3=1.5*sqrt(pi)*beta*scipy.special.gamma((5.*n+5.)/4.)/scipy.special.gamma((5.*n+7.)/4.)
    part4=6.*sqrt(pi)*(1-beta)*scipy.special.gamma(n+2.)/scipy.special.gamma(n+2.5)

    return part1-part2-part3-part4

def Pcal(beta):
    part1=sqrt(pi)*(5.*n-3.)*beta/16.*scipy.special.gamma((5.*n+1.)/4.)/scipy.special.gamma((5.*n+7.)/4.)
    part2=sqrt(pi)*(3.-n)*(1.-beta)/2.*scipy.special.gamma(n+1.)/scipy.special.gamma(n+2.5)
    part3=0.5*sqrt(pi)*beta*scipy.special.gamma((5.*n+5.)/4.)/scipy.special.gamma((5.*n+7.)/4.)
    part4=2.*sqrt(pi)*(1-beta)*scipy.special.gamma(n+2.)/scipy.special.gamma(n+2.5)

    return -part1+part2+part3+part4

def Qcal(beta):
    part1=1.5*sqrt(pi)*beta*scipy.special.gamma((5.*n+5.)/4.)/scipy.special.gamma((5.*n+7.)/4.)
    part2=12.*sqrt(pi)*(1-beta)*scipy.special.gamma(n+2.)/scipy.special.gamma(n+2.5)

    return part1+part2

# previous version of SPQR coefficients:
def Scal_pr(beta):
    return 104./15.*beta-17.*pi/32.*beta-104./15.

def Pcal_pr(beta):
    return 5./32.*pi*beta-8./3.*beta+8./3.

def Qcal_pr(beta):
    return 9.*pi/16.*beta+64./5.-64.*beta/5.

# some unused vertical structure function:
def ex(wrf,r,tc,h):
    a=256.*pi*alpha*tc**4.
 #   print "a= "+str(a)+'\n'
    b=-45.*wrf
#    print "b= "+str(b)+'\n'
    d=360.*chi*tc*wrf*r**3.
#    print "d= "+str(d)+'\n'
    return a*h**3.+b*h**2.+d

def fh1(wrf,r,tc):
    h=scipy.optimize.fsolve(ex, 0.1, args=(wrf,r,tc))
    print h
#

# disc thickness
def fh(wrf,r,tc):
    h=np.zeros(3, dtype=double)
    a=256.*pi*alpha*tc**4.
    b=-45.*wrf
    d=360.*chi*tc*wrf*r**3.
    p=-b*b/(3.*a*a)
    q=(2.*b**3.+27.*a*a*d)/(27.*a**3.)
    k=-1.

    if((p<0.)&((4.*p**3+27.*q*q)<=0)):
       t1=2.*sqrt(-p/3.)*cos(1./3.*arccos(3.*q/(2.*p)*sqrt(-3./p)))-b/(3.*a)
       if(t1>0):
           h[0]=t1
       else:
           h[0]=1000000.
       t2=2.*sqrt(-p/3.)*cos(1./3.*arccos(3.*q/(2.*p)*sqrt(-3./p))-2.*pi/3.)-b/(3.*a)
       if(t2>0):
            h[1]=t2
       else:
           h[1]=1000000.
       t3=2.*sqrt(-p/3.)*cos(1./3.*arccos(3.*q/(2.*p)*sqrt(-3./p))-4.*pi/3.)-b/(3.*a)
       if(t3>0):
           h[2]=t3
       else:
           h[2]=1000000.
       k=h.min()
    if((p<0.)&((4.*p**3+27.*q*q)>0)):
        k=-2.*fabs(q)/q*sqrt(-p/3.)*cosh(1./3.*arccosh(-3.*fabs(q)/(2.*p)*sqrt(-3./p)))-b/(3.*a)
    if(p>0.):
        k=-2.*sqrt(p/3.)*sinh(1./3.*arcsinh(3.*q/(2.*p)*sqrt(3./p)))-b/(3.*a)
    return k

# surface density
def ftau(wrf,r,tc):
    h=fh(wrf,r,tc)
    return wrf*r**3.*hvert**2./(alpha*h**2.)

def xiinf(hin, muu, mdott):
    return (alpha*lam*muu**2/mdott*hin)**(2./9.)/raffugen(muu, mdott)

def raffugen(muu, mdott):
    return (lam*muu**2/mdott)**(2./7.)*2.**(-1./7.)

def peqgen(muu, mdott):
    return pstar*(lam*muu**2/mdott)**(3./7.)*2.**(-3./14.)
