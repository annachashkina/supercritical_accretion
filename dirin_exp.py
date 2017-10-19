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
import bocon_general as b
import ssdisk as ss
import readkit as rk
#import readkit as rk

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

def G(n):
    return sqrt(pi)*scipy.special.gamma(n+1.)/scipy.special.gamma(n+3./2.)

def par0(xs, ys):
    n=ys[0]*(xs[1]**2-xs[2]**2)+ys[1]*(xs[2]**2-xs[0]**2)+ys[2]*(xs[0]**2-xs[1]**2)
    d=ys[0]*(xs[1]-xs[2])+ys[1]*(xs[2]-xs[0])+ys[2]*(xs[0]-xs[1])

    return n/2./ d

def pl(x, norm, pslope):
    return norm*x**(-pslope)

def lgmin(x):
    y=(log10(x)).min()
    y=floor(y)
    return 10.**y

def lgmax(x):
    y=(log10(x)).max()
    y=ceil(y)
    return 10.**y

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


def domega(omega,r,tau,mdot,wrf,tc):
    return -0.5*omega/r+alpha*tau*(omega**2.-1.)/(2.*mdot*r**0.5)+64.*pi/3.*G(n)*(n+1.)*epsilon*(psi-1.)*omega*tc**4.*r**2./tau/mdot+wrf*r**0.5/mdot

def dmdot(r,tc,tau):
    return 64.*pi/3.*G(n)*(n+1.)*epsilon*r**2*tc**4./tau

def dwrf(tau,omega,r):
    return alpha*tau*(omega**2.-1)/r**2.


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

def dtemp(tau,tc,wrf,r,omega,mdot):
    bet=beta(tau,tc,wrf)
    part1=(8.-6.*bet)
    part2=(1.-3.*bet)/tau*dtau(tau,tc,wrf,r,omega,mdot)
    part3=alpha*tau*(omega*omega-1.)*(1.+bet)/wrf/r/r
    part4=3.*(1.-bet)/r
    return tc/part1*(part2+part3-part4)

def beta(tau,tc,wrf):
 #   print "beta= "+str(8.*chi/5.*alpha*tau*tc/wrf)
#    jjii=raw_input()
    return chi*alpha*tau*tc/wrf*2.*(n+1.)/(2.*n+3.)

def comegaga(tc,tau,r,wrf):
    return 64.*pi/3.*G(n)*(n+1.)*tc**4./(tau*r*wrf)
 #   return 512.*pi/9.*tc**4./(tau*r*wrf)

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

def Scal_pr(beta):
    return 104./15.*beta-17.*pi/32.*beta-104./15.

def Pcal_pr(beta):
    return 5./32.*pi*beta-8./3.*beta+8./3.

def Qcal_pr(beta):
    return 9.*pi/16.*beta+64./5.-64.*beta/5.


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


def func1(wrf,r,tc):
    print "solution 1 "+str(fh(wrf,r,tc))+'\n'
    print "solution 2 "+str(fh1(wrf,r,tc))+'\n'
    return 0.

def fh(wrf,r,tc):
    h=np.zeros(3, dtype=double)
    a=256.*pi*alpha*tc**4.
 #   print "a= "+str(a)+'\n'
    b=-45.*wrf
#    print "b= "+str(b)+'\n'
    d=360.*chi*tc*wrf*r**3.
#    print "d= "+str(d)+'\n'
    p=-b*b/(3.*a*a)
 #   print "p= "+str(p)+'\n'
    q=(2.*b**3.+27.*a*a*d)/(27.*a**3.)
#    print "q= "+str(q)+'\n'
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
 #      print "4p^3+27q^2 = "+str(4.*p**3+27.*q*q)+'\n'
 #      print "t1= "+str(h[0]/r)+" t2= "+str(h[1]/r)+" t3= "+str(h[2]/r)+'\n'
#       print "minimum is "+str(k/r)+'\n'
 #      koko=raw_input()
    if((p<0.)&((4.*p**3+27.*q*q)>0)):
        k=-2.*fabs(q)/q*sqrt(-p/3.)*cosh(1./3.*arccosh(-3.*fabs(q)/(2.*p)*sqrt(-3./p)))-b/(3.*a)
 #       print "2 "+str(k/r)+'\n'
    if(p>0.):
        k=-2.*sqrt(p/3.)*sinh(1./3.*arcsinh(3.*q/(2.*p)*sqrt(3./p)))-b/(3.*a)
 #       print "3 "+str(k/r)+'\n'
#    print "k "+str(k/r)+'\n'
 #   if(k<0.):
 #       print"AAAAAAAAAAAAAAAAAAAAAAAAAAA "+str(k)+'\n'
 #       print "4p^3+27q^2 = "+str(4.*p**3+27.*q*q)+'\n'
 #       huhi=raw_input()
    return k


def ftau(wrf,r,tc):
    h=fh(wrf,r,tc)
    return wrf*r**3.*hvert**2./(alpha*h**2.)

def rastr(rrin, qeq):
    ra=b.rafun()
    rin=ra*rrin
    t=1
    rout=100.*rin
    ddr=-1.e-4
    mdot=mdotglobal
    
    if(mdot<1.5):
        ddr*=(mdot/1.5)
    ddr/=(1.+0.5*(b.peq()/ps))
    defac=0.99

    r=rout
    omega=1.
    oprev=omega
    rprev=r

    htormax=0.
#    mdot=mdotglobal
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
#    h=fh(wrf,r,tc)
#    tau=ftau(wrf,r,tc)    

    while(r>=rin):
        dr=ddr*(r-rin*defac)
        r1=r+dr/2.
 #       print r/rin
        
#        h=fh(wrf,r,tc)
        h=hvert*sqrt(r**3*wrf/alpha/tau)
        
        
        if(isnan(tau)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.)
        omega1=omega+domega(omega,r,tau,mdot,wrf,tc)*dr/2.
        wrf1=wrf+dwrf(tau,omega,r)*dr/2.
        tau1=tau+dtau(tau,tc,wrf,r,omega,mdot)*dr/2.
        tc1=tc+dtemp(tau,tc,wrf,r,omega,mdot)*dr/2.
#        tau1=ftau(wrf1,r1,tc1)
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
#            print r1/rin
#            print mdot
#            jiu=raw_input()
#        else:
#            mdot=mdotglobal
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
#    print omega_in
 #   jjoij=raw_input()
    wrf_in=b.fwrfin(rin, hint,mdotin)
#    print wrf_in
#    jjoij=raw_input()
    
#    print "rastr: omega = "+str(omega)+"; h = "+str(h)+"; htormax = "+str(htormax)+"; mdotin = "+str(mdot)+"; wrf = "+str(wrf)+"\n"
#    print "rastr: omega_in = "+str(omega_in)+"; h_in = "+str(hint)+" ; wrf_in = "+str(wrf)+"\n"

    return oint, hint, tint, htormax, mdotin, wint

#    oint=interp1d(rar, oar, bounds_error=False, fill_value=(oar[rar.argmin()], 1.))\#
#    hint=interp1d(rar, har, bounds_error=False, fill_value=(oar[rar.argmin()], 1.))#
#    tint=interp1d(rar, tauar, bounds_error=False, fill_value=(oar[rar.argmin()], 1.))   
#    return oint(rin), hint(rin), tint(rin), htormax, mdotin,wrfin


def rastr1(rrin, qeq):
    ra=b.rafun()
    rin=ra*rrin
#    print rin
#    kok=raw_input()
    t=1
    rout=100.*rin
    ddr=-3.e-5
    mdot=mdotglobal
    print "rastr1: mdot= "+str(mdot)+" mu= "+str(mu)+" \n"
#    huh=raw_input()
    
    if(mdot<1.5):
        ddr*=(mdot/1.5)
    ddr/=(1.+0.1*(b.peq()/ps))
    defac=0.99

    r=rout
    omega=1.
    drout=1e-2
    rlast=rout
    oprev=omega
    rprev=r
    

    rl=[]
    ol=[]
    tl=[]
    hl=[]
    bb=[]
    taussa=[]
    taussb=[]
    hssa=[]
    hssb=[]
    dt=[]
    do=[]
    dw=[]
    dta=[]
    mdot123=[]
    P1m=[]
    P3m=[]
    qqra=[]
    qqad=[]
    qqpl=[]
    wwrf=[]
    ttc=[]

    htormax=0.
#    mdot=mdotglobal
    wrf=2.*mdot/r**2*omega*(sqrt(r)-qeq*sqrt(rin))
    tau=4./chi**0.8*(pi/9.)**0.2*mdot**0.6/alpha**0.8/rout**0.6*(1.-qeq*(rin/rout)**0.5)**0.6
 #   tau=ftau(wrf,r,tc)
    h=hvert*sqrt(r**3*wrf/alpha/tau)
    tc=b.ctemp(h,wrf,tau)
#    print h/r
    hprev=h
    wprev=wrf
    tauprev=tau
 
    f=open('rastr_new_1.txt','w')

    while(r>=rin):
        dr=ddr*(r-rin*defac)
        r1=r+dr/2.
#        print r/rin

        h=hvert*sqrt(r**3*wrf/alpha/tau)
#        print h/r

        if(isnan(tau)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.)
        omega1=omega+domega(omega,r,tau,mdot,wrf,tc)*dr/2.
        wrf1=wrf+dwrf(tau,omega,r)*dr/2.
        tau1=tau+dtau(tau,tc,wrf,r,omega,mdot)*dr/2.
        tc1=tc+dtemp(tau,tc,wrf,r,omega,mdot)*dr/2.

        beta1=beta(tau1,tc1,wrf1)

        
        if((r*r)>(9.*tau1/(4.*64.*pi*tc**4.))):

            mdot1=mdot+dmdot(r,tc,tau)*dr/2.
        else:
            mdot1=mdot
 
        if(wrf1<=0.):
            print "negative stress! wrf = "+str(wrf1)
            return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.),sqrt(-1.), sqrt(-1.)  

        h1=hvert*sqrt(r1**3*wrf1/alpha/tau1)
        
        if(isnan(tau1)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.)
        omega=omega+domega(omega1,r1,tau1,mdot1,wrf1,tc1)*dr
        wrf=wrf+dwrf(tau1,omega1,r1)*dr
        tc=tc+dtemp(tau1,tc1,wrf1,r1,omega1,mdot1)*dr

        tau=tau+dtau(tau1,tc1,wrf1,r1,omega1,mdot1)*dr

        P4=2./3.
        P3=1.77526e-5*alpha*tc1*tau1/wrf1
        P1=56329.9*h1*h1/tc1*r1**(-3.)

 #       Qrad=128./9.*tc1**4./tau1
 
        bet=beta(tau1,tc1,wrf1)

        Stau=Scal(bet)
        Swrf=Pcal(bet)
        Stemp=Qcal(bet)

        part1=Stau/tau1*dtau(tau1,tc1,wrf1,r1,omega1,mdot1)
        part2=Swrf/wrf1*dwrf(tau1,omega1,r1)
        part3=Stemp/tc1*dtemp(tau1,tc1,wrf1,r1,omega1,mdot1)
        part4=3./r*Swrf

        Qrad=-16./3.*G(n)*(n+1.)*tc1**4.*4.*pi/(tau1)


        Qadv=2./(G(n+1.))*mdot1*wrf1/(alpha*r1*tau1)*(part1+part2+part3+part4)

        Qplus=wrf1*r1**(-1./2.)*(domega(omega1,r1,tau1,mdot1,wrf1,tc1)-3./2.*omega1/r1)

        r+=dr
        
        if((r*r)>(9.*tau1/(4.*64.*pi*tc**4.))):
            mdot=mdot+dmdot(r1,tc1,tau1)*dr

        if((h/r)>htormax):
            htormax=h/r
            rmax=r
        if(r<(rlast/(1.+drout))):
            ttc.append(tc)
            qqra.append(Qrad)
            wwrf.append(wrf)
            qqpl.append(Qplus)                                                     
            qqad.append(Qadv)                                                    
            P1m.append(P1)
            P3m.append(P3)
            rl.append(r/rin)
            ol.append(omega)
            tl.append(tau)
            hl.append(h/r)
            dt.append(abs(dtemp(tau1,tc1,wrf1,r1,omega1,mdot1))/tc1)
            do.append(abs(domega(omega1,r1,tau1,mdot1,wrf1,tc1)))
            dw.append(abs(dwrf(tau1,omega1,r1))/wrf1)
            dta.append(abs(dtau(tau1,tc1,wrf1,r1,omega1,mdot1))/tau1)
            bb.append(beta1)
            rlast=r
            taussa.append(ss.ftaussa(r,mdot, alpha, rin))
            taussb.append(ss.ftaussb(r,mdot, alpha, rin))
            hssa.append(ss.hssa(mdot,r, rin))
            hssb.append(ss.hssb(r,mdot, alpha, rin))
            mdot123.append(mdot)
            f.write(str(r/rin)+' '+str(omega)+' '+str(h/r)+' '+str(tau)+' '+str(tc)+' '+str(wrf)+' '+str(ss.ftaussa(r,mdot, alpha, rin))+' '+str(ss.ftaussb(r,mdot, alpha, rin))+' '+str(Qadv)+' '+str(Qrad)+' '+str(Qplus)+' '+str(ss.hssa(mdot,r, rin))+' '+str(ss.hssb(r,mdot, alpha, rin))+' '+str(mdot)+'\n')


    mdotin=mdot
    tauin=tau
#    print "tau= "+str(tauin)+'\n'
#    print mdotin
#    jjoij=raw_input()
    wrfin=wrf
#    print wrfin
#    jjoij=raw_input()
    f.close()

    tcar=asarray(ttc, dtype=double)
    qplus=asarray(qqpl, dtype=double)
    wrfar=asarray(wwrf, dtype=double)
    qadv=asarray(qqad, dtype=double)
    qrad=asarray(qqra, dtype=double)
    p1ar=asarray(P1m, dtype=double)
    p3ar=asarray(P3m, dtype=double)
    rar=asarray(rl, dtype=double)
    oar=asarray(ol, dtype=double)
    tauar=asarray(tl, dtype=double)
    har=asarray(hl, dtype=double)
    taua=asarray(taussa, dtype=double)
    taub=asarray(taussb, dtype=double)
    hrssa=asarray(hssa, dtype=double)
    hrssb=asarray(hssb, dtype=double)
    mdar=asarray(mdot123, dtype=double)
    dtem=asarray(dt, dtype=double)
    dome=asarray(do, dtype=double)
    dwr=asarray(dw, dtype=double)
    dtaau=asarray(dta, dtype=double)
    bbeta=asarray(bb, dtype=double)


    oint=(omega-oprev)*(rin-rprev)/(r-rprev)+oprev
    hint=(h-hprev)*(rin-rprev)/(r-rprev)+hprev
    tint=(tau-tauprev)*(rin-rprev)/(r-rprev)+tauprev
    wint=(wrf-wprev)*(rin-rprev)/(r-rprev)+wprev

    pp=1./4.*(n+1.)*G(n)
    
    plt.clf()
    fig=figure()
    plt.subplot (3, 3, 1)
    plot(rar*rin*rrin*206746., p1ar, color='black',label='P1')
    plot(rar*rin*rrin*206746., p3ar, color='green',label='P3')
    plot(rar*rin*rrin*206746., p1ar*0.+2./3., color='red',label='P2')
    plot(rar*rin*rrin*206746., p3ar*0.+pp, color='blue',label='P4')
    ylim(0.,7.)
    ylabel('$P1,2,3,4$')
    xlabel('$r$')
    legend()
#    yscale('log')
#    xscale('log')


    plt.subplot (3, 3, 2)
    plot(rar*rin*rrin*206746., qadv/qplus, color='red',label='Qadv/Qvis')
    plot(rar*rin*rrin*206746., qrad/qplus, color='green',label='Qrad/Qvis')
 #   xlim(0.01,5.e10)
    ylabel('$P$')
    xlabel('$r$')
#    yscale('log')
    xscale('log')
    legend()

    plt.subplot (3, 3, 3)
    plot(rar*rin*rrin*206746., taua, color='blue',label='rad')
    plot(rar*rin*rrin*206746., taub, color='green',label='gas')
    plot(rar*rin*rrin*206746., tauar, color='red',label='calc')
    ylabel(r'$\tau$')
    xlabel('$r$')
    yscale('log')
    xscale('log')
    legend()

    plt.subplot (3, 3, 4)
    plot(rar*rin*rrin*206746., har, color='red',label='calc')
    plot(rar*rin*rrin*206746., hrssa/rar/rin, color='blue',label='rad')
    plot(rar*rin*rrin*206746., hrssb/rar/rin, color='green',label='gas')
    ylabel('$h/r$')
    xlabel('$r$')
    xscale('log')
    yscale('log')
    legend()

    plt.subplot (3, 3, 5)
    plot(rar*rin*rrin*206746., wrfar*2.56787e+21, color='red')
    ylabel(r'$W_{rf}$')
    xlabel('$r$')
    xscale('log')
    yscale('log')

    plt.subplot (3, 3, 6)
    plot(rar*rin*rrin*206746., mdar/mdotglobal, color='red')
    ylabel(r'$\dot M/\dot M_{0}$')
    xlabel('$r$')
    xscale('log')
#    yscale('log')

    plt.subplot (3, 3, 7)
    plot(rar*rin*rrin*206746., oar, color='red')
    ylabel(r'$\Omega/\Omega_{\rm K}$')
    xlabel('$r$')
    xscale('log')
#    yscale('log')

    plt.subplot (3, 3, 8)
    plot(rar*rin*rrin*206746., tcar*9.6e7, color='red')
    ylabel(r'$T_{c}$')
    xlabel('$r$')
    xscale('log')
    yscale('log')
    fig.set_size_inches(15, 15)
    savefig('all.eps')
    

    plt.clf()
    plot(rar, fabs(qadv/qplus), color='blue')
    ylabel('$Qadv/Qplus$')
    xlabel('$r/rin$')
 #   ylim(0.,10.)   
    yscale('log')
    xscale('log')
    savefig('Qadv.eps')

    plt.clf()
    plot(rar, qrad/qplus, color='green')
    ylabel('$Qrad/Qplus$')
    xlabel('$r/rin$')   
 #   yscale('log')
    xscale('log')
    savefig('Qrad.eps')


    plt.clf()
    plot(rar, qrad/qplus, color='green',label='Qrad/Qplus')
    plot(rar, qadv/qplus, color='red',label='Qadv/Qplus')
    ylabel('$Q/Qplus$')
    xlabel('$r/rin$')   
 #   yscale('log')
    xscale('log')
    legend()
    savefig('Q.eps')



    plt.clf()
    plot(rar, p1ar, color='blue',label='P1')
    plot(rar, p3ar, color='green',label='P3')
    ylabel('$P$')
    xlabel('$r/rin$')
    ylim(0.,10.)
    legend()
 #   yscale('log')
#    xscale('log')
    savefig('P.eps')

    plt.clf()
    plot(rar, p1ar, color='blue',label='P1')
    plot(rar, p3ar, color='green',label='P3')
    ylabel('$P$')
    xlabel('$r/rin$')
    yscale('log')
    xscale('log')
    legend()
    savefig('Plog.eps')


    plt.clf()
    plot(bbeta, p1ar, color='blue')
    plot(bbeta, p3ar, color='green')
    ylabel('$P$')
    xlabel('$r/rin$')
    ylim(0.,10.)
    
 #   yscale('log')
#    xscale('log')
    savefig('P_beta.eps')

    plt.clf()
    plot(rar, taua, color='blue')
    plot(rar, taub, color='green')
    plot(rar, tauar,linestyle='dotted', color='red')
    plot([rin/rin],[tint], 'o')
    ylabel(r'$\tau$')
    xlabel('$r/rin$')
    yscale('log')
    xscale('log')
    savefig('tau.eps')

    plt.clf()
    plot(rar, oar, color='k')
    plot([rin/rin],[oint], 'o')
    ylabel(r'$\omega$')
    xlabel('$r/rin$')
    xscale('log')
    savefig('omega.eps')


    plt.clf()
    plot(rar, har,linestyle='dotted', color='red')
    plot(rar, hrssa/rar/rin, color='blue')
    plot(rar, hrssb/rar/rin, color='green')
    plot([rin/rin],[hint/rin], 'o')
    ylabel('$h/r$')
    xlabel('$r/rin$')
    xscale('log')
    yscale('log')
    savefig('har.eps')

    plt.clf()
    plot(rar, mdar, color='k')
    ylabel(r'$\dot{m}$')
    xlabel('$r/rin$')
    xscale('log')
    yscale('log')
    savefig('mdot.eps')

    if(isnan(omega)|(r>rin)):
        return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.),sqrt(-1.), sqrt(-1.) 
    h_in=b.fhin(rin, tau,mdotin)
#    print h_in
#    jjoij=raw_input()
    omega_in=b.oin(rin, hint,mdotin)
#    print omega_in
#    jjoij=raw_input()
    wrf_in=b.fwrfin(rin, hint,mdotin)
#    print wrf_in
#    jjoij=raw_input()
    
    print "rastr1: omega = "+str(omega)+"; h = "+str(h)+"; htormax = "+str(htormax)+"; mdotin = "+str(mdot)+"; wrf = "+str(wrf)+"\n"
    print "rastr1: omega_in = "+str(omega_in)+"; h_in = "+str(h_in)+" ; wrf_in = "+str(wrf)+"\n"

    return rar*rin, tauar

def rmesh_qeq(newmu, newmdot, newps):

    b.parset(newmu=newmu, newmdot=mdotglobal,newps=newps,neweta=0.,newalpha=0.1)
    parset(newmu=newmu, newmdot=mdotglobal,newps=newps,neweta=0.,newalpha=0.1)


    qeqmin=0.85
    qeqmax=0.95
    nr=30

    qeq=(qeqmax-qeqmin)*arange(nr)/double(nr-1)+qeqmin
  

    ocalc=np.zeros(nr, dtype=double)
    obc=np.zeros(nr, dtype=double)
    hcalc=np.zeros(nr, dtype=double)
    bmax=np.zeros(nr, dtype=double)
    hbc=np.zeros(nr, dtype=double)
    tcalc=np.zeros(nr, dtype=double)
    tbc=np.zeros(nr, dtype=double)
    so=np.zeros(nr, dtype=double)
    st=np.zeros(nr, dtype=double)
    sh=np.zeros(nr, dtype=double)
    smdot=np.zeros(nr, dtype=double)

    t=1
    k1=-0.1

    f=open('q'+str(t)+'.txt','w')


    for k in arange(nr):
        xi, tc=ordiv(newmu, newmdot, newps,0.0,qeq[k])
 #      tc=ordiv_tc(newmu, newmdot, newps,0.,qeq,xi)
        print "xi= "+str(xi)+'\n'
        print "qeq= "+str(qeq[k])+'\n'
        rin=xi*b.rafun()
        oo,hh,tt,hrmax,mdotin,ww=rastr(xi,qeq[k],tc)
        oin=b.oin(rin, hh,mdotin)
        hin=b.fhin(rin, tt, mdotin)
        wrfin=b.fwrfin(rin, hh,mdotin)
        ocalc[k]=oo
        hcalc[k]=hh
        tcalc[k]=ww
        obc[k]=oin
        hbc[k]=hin
        tbc[k]=wrfin
        so[k]=sign(oo-oin)
        st[k]=sign(tt-tt)
        sh[k]=sign(hin-hh)
        smdot[k]=mdotin

        
        f.write(str(qeq[k])+' '+str(xi)+' '+str(oo)+' '+str(oin)+' '+str(hh)+' '+str(hin)+' '+str(ww)+' '+str(wrfin)+' '+str(mdotin)+'\n')
           
    f.close()

    plt.clf()
    plot(qeq, smdot, color='k')
    ylabel('$\dot m$')
    xlabel('$qeq$')
    savefig('mdot'+str(t)+'.eps')

    plt.clf()
    plot(qeq, ocalc, color='k')
    plot(qeq, obc, color='r')
    ylabel('$\omega$')
    xlabel('$qeq$')
    ylim(0.1, 1.5)
    savefig('oos'+str(t)+'.eps')
 
    plt.clf()
    plot(qeq, hcalc, color='k')
    plot(qeq, hbc, color='r')
    ylabel('$h$')
    xlabel('$qeq$')
    yscale('log')
    savefig('h'+str(t)+'.eps')

    plt.clf()
    plot(qeq, tcalc, color='k')
    plot(qeq, tbc, color='r')
    ylabel('$wrf$')
    xlabel('$qeq$')
    yscale('log')
    savefig('wrf'+str(t)+'.eps')


def rmesh(newmu, newmdot, newps):

    b.parset(newmu=newmu, newmdot=mdotglobal,newp=newps,neweta=0.,newalpha=0.1)
    parset(newmu=newmu, newmdot=mdotglobal,newp=newps,neweta=0.,newalpha=0.1)


    qeqmin=0.1
    qeqmax=1.
    nr=10
    pr=10.
    ximin=0.1
    ximax=1.5
    tcmin=0.01
    tcmax=2.

#    qeq=(qeqmax-qeqmin)*arange(nr)/double(nr-1)+qeqmin
 #   xi=(ximax-ximin)*arange(nr)/double(nr-1)+ximin
#    tc=(tcmax-tcmin)*arange(nr)/double(nr-1)+tcmin
  

    ocalc=np.zeros([nr,nr,nr], dtype=double)
    obc=np.zeros([nr,nr,nr], dtype=double)
    hcalc=np.zeros([nr,nr,nr], dtype=double)
    bmax=np.zeros([nr,nr,nr], dtype=double)
    hbc=np.zeros([nr,nr,nr], dtype=double)
    tcalc=np.zeros([nr,nr,nr], dtype=double)
    tbc=np.zeros([nr,nr,nr], dtype=double)
    so=np.zeros([nr,nr,nr], dtype=double)
    st=np.zeros([nr,nr,nr], dtype=double)
    sh=np.zeros([nr,nr,nr], dtype=double)
    smdot=np.zeros([nr,nr,nr], dtype=double)

    t=1
    k1=-0.1

    f=open('q'+str(t)+'.txt','w')

    xi=ximin
    tc=tcmin
    qeq=qeqmin
    i=0
    
    for k in arange(nr):
        xi=xi+(ximax-ximin)/pr
        tc=tcmin
        for w in arange(nr):
            tc=tc+(tcmax-tcmin)/pr
            qeq=qeqmin
            for q in arange(nr):
                qeq=qeq+(qeqmax-qeqmin)/pr
                print 'qeq= '+str(qeq)+' tc= '+str(tc)+' xi= '+str(xi)+'\n'
                rin=xi*b.rafun()
                oo,hh,tt,hrmax,mdotin,ww=rastr(xi,qeq,tc)
                oin=b.oin(rin, hh,mdotin)
                hin=b.fhin(rin, tt, mdotin)
                wrfin=b.fwrfin(rin, hh,mdotin)
                ocalc[k,w,q]=oo
                hcalc[k,w,q]=hh
                tcalc[k,w,q]=ww
                obc[k,w,q]=oin
                hbc[k,w,q]=hin
                tbc[k,w,q]=wrfin
                so[k,w,q]=oo-oin
                st[k,w,q]=ww-wrfin
                sh[k,w,q]=-hin+hh
                smdot[k,w,q]=mdotin
                print i
                i+=1

        
                f.write(str(qeq)+' '+str(xi)+' '+str(tc)+' '+str(oo)+' '+str(oin)+' '+str(hh)+' '+str(hin)+' '+str(ww)+' '+str(wrfin)+' '+str(mdotin)+' '+str(oo-oin)+' ' +str(ww-wrfin)+' ' +str(hh-hin)+'\n')
           
    f.close()

def experimental(newmu,newps):
    global xiest, qeqest


    tstart=time.time()

    xiest=0.4
    xii=[]
    qeqq=[]
    md=[]
    newmdot1=10.
    newmdot2=10000.
    nmdot=16
    mdar=(newmdot2/newmdot1)**(arange(nmdot)/double(nmdot-1))*newmdot1
    refrad1 = 500.
    refrad2 = 1000.

    mumin=1. 
    mumax=100.
    nmu=15
    muar=(mumax/mumin)**(arange(nmu)/double(nmu-1))*mumin


    
    xx=xiest
    qq=qeqest
    f=open('experimental.txt','w')
    for ku in arange(nmu):
        xiest=xx
        qeqest=qq
        for i in arange(nmdot):
            newmu=muar[ku]
            newmdot=mdar[i]
            b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
            parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
            x=ordiv_smart(newmu, newmdot, newps)
            xiest=x[0]
            qeqest=x[1]
            if(i==0):
                xx=xiest
                qq=qeqest
            xii.append(xiest)
            qeqq.append(qeqest)
            md.append(newmdot)
            r, tau = rastr1(xiest, qeqest)
            print "inner radius "+str(r.min())
            taufun=interp1d(r, tau)
            tauref1=taufun(refrad1)
            tauref2=taufun(refrad2)
            f.write(str(newmu)+' '+str(newmdot)+' '+str(xiest)+' '+str(qeqest)+' '+str(tauref1)+' '+str(tauref2)+'\n')
            print str(newmu)+' '+str(newmdot)+' '+str(xiest)+' '+str(qeqest)+' '+str(tauref1)+' '+str(tauref2)+'\n'
           
    f.close()

def ordiv_smart(newmu, newmdot, newps):
    b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
    parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
#    print 'here'
#    print "ordiv_smart: xiest = "+str(xiest)+", qeqest = "+str(qeqest)
    tstart=time.time()
    co=scipy.optimize.root(vrapper,(xiest,qeqest),method='hybr',jac=None,tol=1e-4,callback=None,options=None)
    tend=time.time()
    print "co.x[0]= "+str(co.x[0])+' \n'
    print "co.x[1]= "+str(co.x[1])+' \n'
    print "calculation took "+str(tend-tstart)+"s = "+str((tend-tstart)/60.)+"min"

    return co.x

def vrapper(arg):
    xi=arg[0]
    qeq=arg[1]
    a,b=doffwrfin(xi,qeq)
#    b=doffo(xi,qeq,tc)
    print "vrap a= "+str(a)+" vrap b= "+str(b)+'\n'
    print "xi= "+str(xi)+" qeq= "+str(qeq)+'\n'

    return (a,b)

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

def searcon(a, b,qeq,tc):
    x1=a
    x2=b

    tol=1e-2
    
    while(abs(x1-x2)>tol):
        x=(x1+x2)/2.
        f=doffo(x,qeq,tc)
        if(isnan(f)):
            x2=x
        else:
            x1=x

    return x1

def doffo(xi,qeq,tc):
    ra=b.rafun()
    rin=ra*xi
    oin,hin,tauin,hrmax,mdotin,wrfin=rastr(xi,qeq,tc)
    oo=b.oin(rin, hin,mdotin)
    return oo-oin


def searchsign_old(x1, x2, f1,tol,qeq,tc):
    if(abs(x2-x1)<tol):
        print "too close"
        return sqrt(-1.)
    x=(x1+x2)/2.
    f=doffo(x,qeq,tc)
#    print f
    if((f*f1)<0.):
#        print "1"+str(x)
        return x
    else:
        xx1=searchsign_old(x1, x, f1,tol,qeq,tc)
        print "xx1= "+str(xx1)
        if(isfinite(xx1)):
#            print "2"+str(xx1)
            return xx1
        
        xx2=searchsign_old(x, x2, f1,tol,qeq,tc)
 #       print "3"+str(xx2)
        return xx2

def searconleft(xi1, xi2,qeq,tc):
    x1=xi1
    x2=xi2

    tol=1e-2
    
    while(abs(x1-x2)>tol):
        x=(x1+x2)/2.
        f=doffo(x,qeq,tc)
        if(isnan(f)):
            x1=x
        else:
            x2=x

    return x2


def rmmplot_old(fname):

    qeq,xi,tc,oo,oin,hh,hin,ww,wrfin,mdotin,difo,dwrf,dh=rk.rmmread(fname)

 #   mu2,md2,xi,o,htor,tau=rk.rmmread(fname)

    wuu=unique(qeq)
    wdd=unique(xi)
    wtc=unique(tc)
    nuu=size(wuu)
    ndd=size(wdd)
    ntc=size(wtc)

#    print qeq

  #  wuu=unique(mu2)
#    wdd=unique(md2)
#    nuu=size(wuu)
#    ndd=size(wdd)

    print nuu
    print ntc

    qeqar=reshape(asarray(qeq,dtype=double),[ntc,ndd,nuu])
    xiar=reshape(asarray(xi,dtype=double),[ntc,ndd,nuu])
    tcar=reshape(asarray(tc,dtype=double),[ntc,ndd,nuu])
    difoar=reshape(asarray(difo,dtype=double),[ntc,ndd,nuu])
    dwrfar=reshape(asarray(dwrf,dtype=double),[ntc,ndd,nuu])
    dhar=reshape(asarray(dh,dtype=double),[ntc,ndd,nuu])
    
#    muar=reshape(asarray(mu2,dtype=double),[nuu,ndd])
#    mdar=reshape(asarray(md2,dtype=double),[nuu,ndd])
#    xiar=reshape(asarray(xi,dtype=double),[nuu,ndd])
#    har=reshape(asarray(htor,dtype=double),[nuu,ndd])
    print qeqar[0,:,:].min()
    print qeqar[0,:,:].max()

    plt.clf()
    fig=figure()
    contourf(tcar[3,:,:], qeqar[3,:,:], difoar[3,:,:]) #,levels=rlev)
    title('difOmeg')
    c=colorbar()
    contour(tcar[0,:,:], qeqar[0,:,:], difoar[0,:,:],levels=[0.],linestyles='dashed', colors='w')
    contour(tcar[0,:,:], qeqar[0,:,:], dwrfar[0,:,:],levels=[0.],linestyles='dashed', colors='w')
    contour(tcar[0,:,:], qeqar[0,:,:], dhar[0,:,:],levels=[0.],linestyles='dashed', colors='w')

    contour(tcar[1,:,:], qeqar[1,:,:], difoar[1,:,:],levels=[0.], colors='g')
    contour(tcar[1,:,:], qeqar[1,:,:], dwrfar[1,:,:],levels=[0.], colors='g')
    contour(tcar[1,:,:], qeqar[1,:,:], dhar[1,:,:],levels=[0.], colors='g')

    contour(tcar[2,:,:], qeqar[2,:,:], difoar[2,:,:],levels=[0.], colors='b')
    contour(tcar[2,:,:], qeqar[2,:,:], dwrfar[2,:,:],levels=[0.], colors='b')
    contour(tcar[2,:,:], qeqar[2,:,:], dhar[2,:,:],levels=[0.], colors='b')

    contour(tcar[3,:,:], qeqar[3,:,:], difoar[3,:,:],levels=[0.], colors='yellow')
    contour(tcar[3,:,:], qeqar[3,:,:], dwrfar[3,:,:],levels=[0.], colors='yellow')
    contour(tcar[3,:,:], qeqar[3,:,:], dhar[3,:,:],levels=[0.], colors='yellow')

    contour(tcar[4,:,:], qeqar[4,:,:], difoar[4,:,:],levels=[0.], colors='magenta')
    contour(tcar[4,:,:], qeqar[4,:,:], dwrfar[4,:,:],levels=[0.], colors='magenta')
    contour(tcar[4,:,:], qeqar[4,:,:], dhar[4,:,:],levels=[0.], colors='magenta')

    contour(tcar[5,:,:], qeqar[5,:,:], difoar[5,:,:],levels=[0.], colors='black')
    contour(tcar[5,:,:], qeqar[5,:,:], dwrfar[5,:,:],levels=[0.], colors='black')
    contour(tcar[5,:,:], qeqar[5,:,:], dhar[5,:,:],levels=[0.], colors='black')

    contour(tcar[6,:,:], qeqar[6,:,:], difoar[6,:,:],levels=[0.],linestyles='dashed', colors='w')
    contour(tcar[6,:,:], qeqar[6,:,:], dwrfar[6,:,:],levels=[0.],linestyles='dashed', colors='w')
    contour(tcar[6,:,:], qeqar[6,:,:], dhar[6,:,:],levels=[0.],linestyles='dashed', colors='w')

    contour(tcar[7,:,:], qeqar[7,:,:], difoar[7,:,:],levels=[0.],linestyles='dashed', colors='g')
    contour(tcar[7,:,:], qeqar[7,:,:], dwrfar[7,:,:],levels=[0.],linestyles='dashed', colors='g')
    contour(tcar[7,:,:], qeqar[7,:,:], dhar[7,:,:],levels=[0.],linestyles='dashed', colors='g')

    contour(tcar[8,:,:], qeqar[8,:,:], difoar[8,:,:],levels=[0.],linestyles='dashed', colors='b')
    contour(tcar[8,:,:], qeqar[8,:,:], dwrfar[8,:,:],levels=[0.],linestyles='dashed', colors='b')
    contour(tcar[8,:,:], qeqar[8,:,:], dhar[8,:,:],levels=[0.],linestyles='dashed', colors='b')

    contour(tcar[9,:,:], qeqar[9,:,:], difoar[9,:,:],levels=[0.],linestyles='dashed', colors='magenta')
    contour(tcar[9,:,:], qeqar[9,:,:], dwrfar[9,:,:],levels=[0.],linestyles='dashed', colors='magenta')
    contour(tcar[9,:,:], qeqar[9,:,:], dhar[9,:,:],levels=[0.],linestyles='dashed', colors='magenta')


    
#    c.set_label(r'$\xi$',fontsize=16)
#    contour(qeqar[:,0,:], xiar[:,0,:], dwrfar[:,0,:], colors='w')
    ylabel('qeq',fontsize=18)
    xlabel('tc',fontsize=18)
#    xscale('log')
#    yscale('log')
#    xlim([muar.min(),muar.max()])
#    ylim([mdar.min(),mdar.max()])
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    fig.set_size_inches(6, 5)
    savefig(fname+'difOmeg.eps')


    plt.clf()
    fig=figure()
#    contourf(qeqar[:,0,:], tcar[:,0,:], dwrfar[:,0,:]) #,levels=rlev)
    contourf(tcar[3,:,:], qeqar[3,:,:], dwrfar[3,:,:])
    title('dwrf')
    c=colorbar()
#    c.set_label(r'$\xi$',fontsize=16)
    contour(tcar[3,:,:], qeqar[3,:,:], dwrfar[3,:,:],levels=[0.], colors='w')
    xlabel('qeq',fontsize=18)
    ylabel('tc',fontsize=18)
#    xscale('log')
#    yscale('log')
#    xlim([muar.min(),muar.max()])
#    ylim([mdar.min(),mdar.max()])
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    fig.set_size_inches(6, 5)
    savefig(fname+'dwrf.eps')


    plt.clf()
    fig=figure()
    contourf(tcar[3,:,:],qeqar[3,:,:], dhar[3,:,:]) #,levels=rlev)
    title('dh')
    c=colorbar()
    contour(tcar[3,:,:], qeqar[3,:,:], dhar[3,:,:],levels=[0.], colors='w')
#    c.set_label(r'$\xi$',fontsize=16)
#    contour(qeqar[:,0,:], xiar[:,0,:], dwrfar[:,0,:], colors='w')
    xlabel('qeq',fontsize=18)
    ylabel('tc',fontsize=18)
#    xscale('log')
#    yscale('log')
#    xlim([muar.min(),muar.max()])
#    ylim([mdar.min(),mdar.max()])
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    fig.set_size_inches(6, 5)
    savefig(fname+'dh.eps')


 #   plt.clf()
#    fig=figure()
#    contourf(muar, mdar, xiar) #,levels=rlev)
#    title(r'$\xi$')
#    c=colorbar()
#    c.set_label(r'$\xi$',fontsize=16)
#    contour(muar, mdar, har, levels=[0.03,0.1,1.0], colors='w')
#    xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=18)
#    ylabel('$\dot{m}$',fontsize=18)
#    xscale('log')
#    yscale('log')
#    xlim([muar.min(),muar.max()])
#    ylim([mdar.min(),mdar.max()])
#    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
#    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
#    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
#    fig.set_size_inches(6, 5)
#    savefig(fname+'_xi.eps')
#    savefig(fname+'_xi.jpg')


def rmm(newps, neweta):

    tstart=time.time()

    mumin=1. 
    mumax=100.
    mdmin=0.1
    mdmax=1000.
    nmu=15
    nd=16

    muar=(mumax/mumin)**(arange(nmu)/double(nmu-1))*mumin
    mdar=(mdmax/mdmin)**(arange(nd)/double(nd-1))*mdmin

    mu2=np.zeros([nmu,nd], dtype=double)+sqrt(-1.)
    md2=np.zeros([nmu,nd], dtype=double)+sqrt(-1.)
    xi=np.zeros([nmu,nd], dtype=double)+sqrt(-1.)
    qeq=np.zeros([nmu,nd], dtype=double)+sqrt(-1.)
    har=np.zeros([nmu,nd], dtype=double)+sqrt(-1.) # inner disc thickness (absolute)
    betar=np.zeros([nmu,nd], dtype=double)+sqrt(-1.) # maximal relative thickness
    oar=np.zeros([nmu,nd], dtype=double)+sqrt(-1.)
    tauar=np.zeros([nmu,nd], dtype=double)+sqrt(-1.)
    raam=np.zeros([nmu,nd], dtype=double)+sqrt(-1.)
    mddotint=np.zeros([nmu,nd], dtype=double)+sqrt(-1.)

    fname='rmm_op'+str(newps)+'_eta'+str(neweta)
    fout=open(fname+'.dat', 'w')
    
    for ku in arange(nmu):
        for kd in arange(nd):
            mu2[ku,kd]=muar[ku]
            md2[ku,kd]=mdar[kd]
            xi[ku,kd],qeq[ku,kd]=ordiv_smart(muar[ku], mdar[kd], -newps)
            print "xi (mu="+str(muar[ku])+", mdot="+str(mdar[kd])+") = "+str(xi[ku,kd])
            oin,hin,tin,bmax,mdotint,ll=rastr(xi[ku,kd], qeq[ku,kd])
            oar[ku,kd]=oin
            har[ku,kd]=hin
            tauar[ku,kd]=tin
            betar[ku,kd]=bmax
            mddotint[ku,kd]=mdotint
            raam[ku,kd]=raffugen(muar[ku], mdar[kd])
            fout.write(str(muar[ku])+' '+str(mdar[kd])+' '+str(xi[ku,kd])+' '+str(oin)+' '+str(hin/raam[ku,kd]/xi[ku,kd])+' '+str(tin)+' '+str(bmax)+' '+str(mddotint[ku,kd])+'\n')
            print str(muar[ku])+' '+str(mdar[kd])+' '+str(xi[ku,kd])+' '+str(oin)+' '+str(hin/raam[ku,kd]/xi[ku,kd])+' '+str(tin)+' '+str(bmax)

    fout.close()
    
    htor=har/raam/xi

    plt.clf()
    contourf(mu2, md2, xi) #,levels=rlev)
    colorbar()
    contour(mu2, md2, htor, levels=[1e-3,1e-2,1e-1,1.], colors='w')
    xlabel('$\mu$, $10^{30}$G\,cm$^{3}$')
    ylabel('$\dot{m}$')
    xscale('log')
    yscale('log')
    xlim([muar.min(),muar.max()])
    ylim([mdar.min(),mdar.max()])
    savefig(fname+'_xi.eps')
    savefig(fname+'_xi.jpg')

    plt.clf()
    plot(htor, xi, '.')
    xlabel('$H/R$')
    ylabel(r'$\xi$')
    xscale('log')
    yscale('log')
    savefig(fname+'_hh.eps')
    tend=time.time()

    print "RMM: the grid took "+str(tend-tstart)+"s = "+str((tend-tstart)/3600.)+"h"



def xiinf(hin, muu, mdott):
    return (alpha*lam*muu**2/mdott*hin)**(2./9.)/raffugen(muu, mdott)

def raffugen(muu, mdott):
    return (lam*muu**2/mdott)**(2./7.)*2.**(-1./7.)

def peqgen(muu, mdott):
    return pstar*(lam*muu**2/mdott)**(3./7.)*2.**(-3./14.)


def newname(fname):

    mu2,md2,xi,qeq=rk.rmmread(fname)
    print mu2
    nkj=raw_input()
    wuu=unique(mu2)
    wdd=unique(md2)
    nuu=size(wuu)
    ndd=size(wdd)
    

    print nuu
    print ndd

    muar=reshape(asarray(mu2,dtype=double),[nuu,ndd])
    mdar=reshape(asarray(md2,dtype=double),[nuu,ndd])
    xiar=reshape(asarray(xi,dtype=double),[nuu,ndd])
    qeqar=reshape(asarray(qeq,dtype=double),[nuu,ndd])
 
    plt.clf()
    fig=figure()
    contourf(muar, mdar, xiar) #,levels=rlev)
    title(r'$\xi$')
    c=colorbar()
    xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=18)
    ylabel('$\dot{m}$',fontsize=18)
    xscale('log')
    yscale('log')
    xlim([muar.min(),muar.max()])
    ylim([mdar.min(),mdar.max()])
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    fig.set_size_inches(6, 5)
    savefig(fname+'_xi.eps')

    plt.clf()
    fig=figure()
    contourf(muar, mdar, qeqar) #,levels=rlev)
    title(r'$qeq$')
    c=colorbar()
    xlabel('$qeq$, $10^{30}$G\,cm$^{3}$',fontsize=18)
    ylabel('$\dot{m}$',fontsize=18)
    xscale('log')
    yscale('log')
    xlim([muar.min(),muar.max()])
    ylim([mdar.min(),mdar.max()])
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    fig.set_size_inches(6, 5)
    savefig(fname+'_qeq.eps')



    close()


def difeta():

    fname='rmm_exp'
    fname1='rmm_exp_eta1'
    mu2,md2,xi,qeq=rk.rmmread(fname)
    print mu2
    nkj=raw_input()
    wuu=unique(mu2)
    wdd=unique(md2)
    nuu=size(wuu)
    ndd=size(wdd)
    

    print nuu
    print ndd

    muar=reshape(asarray(mu2,dtype=double),[nuu,ndd])
    mdar=reshape(asarray(md2,dtype=double),[nuu,ndd])
    xiar=reshape(asarray(xi,dtype=double),[nuu,ndd])
    qeqar=reshape(asarray(qeq,dtype=double),[nuu,ndd])
 

    mu21,md21,xi1,qeq1=rk.rmmread(fname1)

    wuu1=unique(mu21)
    wdd1=unique(md21)
    nuu1=size(wuu1)
    ndd1=size(wdd1)
    

    print nuu1
    print ndd1

    muar1=reshape(asarray(mu21,dtype=double),[nuu1,ndd1])
    mdar1=reshape(asarray(md21,dtype=double),[nuu1,ndd1])
    xiar1=reshape(asarray(xi1,dtype=double),[nuu1,ndd1])
    qeqar1=reshape(asarray(qeq1,dtype=double),[nuu1,ndd1])
 




    plt.clf()
    fig=figure()
    contourf(muar, mdar, xiar1-xiar) #,levels=rlev)
    title(r'$\Delta \xi$')
    c=colorbar()
    xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=18)
    ylabel('$\dot{m}$',fontsize=18)
    xscale('log')
    yscale('log')
    xlim([muar.min(),muar.max()])
    ylim([mdar.min(),mdar.max()])
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    fig.set_size_inches(6, 5)
    savefig(fname+'_xi_dif.eps')

    plt.clf()
    fig=figure()
    contourf(muar, mdar, qeqar-qeqar1) #,levels=rlev)
    title(r'$qeq$')
    c=colorbar()
    xlabel('$\Delta qeq$, $10^{30}$G\,cm$^{3}$',fontsize=18)
    ylabel('$\dot{m}$',fontsize=18)
    xscale('log')
    yscale('log')
    xlim([muar.min(),muar.max()])
    ylim([mdar.min(),mdar.max()])
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    fig.set_size_inches(6, 5)
    savefig(fname+'_qeq_dif.eps')



    close()








