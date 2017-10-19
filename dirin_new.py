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
import bocon_new as b
import ssdisk as ss
import readkit as rk

#physical parameters
mu=1.
ps=10.
mdotglobal=437.
#coefficients
alpha=0.1
eta=0.0
kt=0.5
epsilon=1.
psi=1.
#dimensionless coefficients
lam=3.9848e10
chi=8.8e-6
pstar=4.33e-5
#vertical structure
tvert=219./1024.
hvert=sqrt(5.)
#program coefficients
tol=1e-7
toff=1e-12

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
    ps=kwargs['newp']
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

def fefftemp(r,omega,h, tau, rin, omegain, wrfin,wrf):
    dwdr=fdomega(r,omega, tau, h, wrf)    
    return (1./(8.*pi)*wrf/sqrt(r)*abs(dwdr-1.5*omega/r))**0.25

def ABCfun(r, h, wrf, tau, omega, mdot):
    aa=2.*wrf*r 
    bb=16.*pi*mdot*r*b.ctemp(h,wrf,tau)**4/wrf/(1.+tvert*tau)
    cc=4.*mdot*omega/r**0.5
    kk=16.*pi*r**2.5*b.ctemp(h,wrf,tau)**4/(1.+tvert*tau)*omega*epsilon*(1.-psi)
    
    return -aa-bb+cc+kk+alpha*tau*(1.-omega**2)

def tausolve(r, omega, wrf, mdot, **kwargs):

    tau1=1.
    tau2=1.e10
    
    h1=hvert*sqrt(wrf*r/tau1/alpha)*r
    f1=ABCfun(r, h1, wrf, tau1, omega, mdot)
    h2=hvert*sqrt(wrf*r/tau2/alpha)*r
    f2=ABCfun(r, h2, wrf, tau2, omega, mdot)

    if((f1*f2)>=0.):
        tauest=kwargs['taurecent']
        tau1=tauest/2.
        tau2=tauest*2.
        h1=hvert*sqrt(wrf*r/tau1/alpha)*r
        f1=ABCfun(r, h1, wrf, tau1, omega, mdot)
        h2=hvert*sqrt(wrf*r/tau2/alpha)*r
        f2=ABCfun(r, h2, wrf, tau2, omega, mdot)

        if((f1*f2)>=0.):
            return sqrt(-1.), sqrt(-1.)
        
    while(abs((tau1-tau2)/(tau1+tau2))>tol):
        tau=sqrt(tau1*tau2)
        h=hvert*sqrt(wrf*r/tau/alpha)*r
        f=ABCfun(r, h, wrf, tau, omega, mdot)
        
        if((f1*f)<0.):
            tau2=tau
        else:
            tau1=tau
    return tau,h

def fdwrf(r,omega, wrf, h, tau, mdot):    
    return -2.*wrf/r+4.*mdot*omega/r**2.5-16.*pi*b.ctemp(h,wrf,tau)**4/(1.+tvert*tau)*(mdot/(r*wrf)-epsilon*sqrt(r)*omega*(1.-psi))

def fdomega(r,omega, tau, h, wrf):

    tt=b.ctemp(h,wrf, tau)
    return 1.5*omega/r-8.*pi*tt**4*sqrt(r)/wrf/(tvert*tau+1.)

def fdmdot(r,wrf,tau,h):
    return 8.*pi*epsilon*(sqrt(r)*b.ctemp(h,wrf,tau))**4/(1.+tvert*tau)

def vert_str(rrin,qeq):
    global oin
    b.parset(newmu=mu, newmdot=mdotglobal,newp=ps,neweta=0.,newalpha=0.1)
    parset(newmu=mu, newmdot=mdotglobal,newp=ps,neweta=0.,newalpha=0.1)

    ra=b.rafun()
    rin=ra*rrin
    print "rin= "+str(rin)
    rout=5.*rin
    ddr=-1.e-5
    mdot=mdotglobal
    print "mu= "+str(mu)+"dotm_global= "+str(mdotglobal)+"ps="+str(ps)
    if(mdotglobal<1.5):
        ddr*=sqrt(mdotglobal/1.5)
    ddr/=(1.+0.1*(b.peq()/ps))
    defac=0.99
    r=rout
    omega=1.0
    drout=1e-2
    rlast=rout
    # lists:
    rl=[]
    ol=[]
    tl=[]
    hl=[]
    mdot123=[]
    tau=1000.
    wrf=2.*mdot/r**2*omega*(sqrt(r)-qeq*sqrt(rin))

    while(r>=rin):
        dr=ddr*r
        r1=r+dr/2.

        tau,h=tausolve(r, omega, wrf, mdot, taurecent=tau)

        if(isnan(tau)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.)
        omega1=omega+fdomega(r, omega, tau, h, wrf)*dr/2.

        wrf1=wrf+fdwrf(r, omega, wrf,h,tau,mdot)*dr/2.

        if((h/r)>1.):
            mdot1=mdot+fdmdot(r, wrf, tau,h)*dr/2.
        else:
            mdot1=mdotglobal


        
        if(wrf1<=0.):
            print "negative stress! wrf = "+str(wrf1)
            print "dr="+str(dr)
            print "mdot= "+str(mdot)
            print "dmdot= "+str(fdmdot(r, wrf, tau,h))
            print "tcent= "+str(b.ctemp(h,wrf,tau))
            print "r="+str(r)
            return sqrt(-1.), sqrt(-1.), sqrt(-1.)
        print r
        tau1,h1=tausolve(r1, omega1, wrf1, mdot1, taurecent=tau)
        

        if(isnan(tau1)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.)
        omega=omega+fdomega(r1, omega1, tau1, h1, wrf1)*dr

        wrf=wrf+fdwrf(r1, omega1, wrf1,h1,tau1,mdot1)*dr
       
        if((h1/r1)>1.):
            mdot=mdot+fdmdot(r1, wrf1, tau1,h1)*dr
        else:
            mdot=mdotglobal
        
        r=(1.+ddr)*r
    print r
    oin=omega
    
    print "the end of the first part"

    mdotin=mdot
    wrfin=wrf

    mdot=mdotglobal
    
    if(mdotglobal<1.5):
        ddr*=sqrt(mdotglobal/1.5)

    ddr/=(1.+0.1*(b.peq()/ps))
    defac=0.99
    r=rout
    omega=1.0
    drout=1e-2
    rlast=rout
    # lists:
    rl=[]
    ol=[]
    tl=[]
    hl=[]
    mdot123=[]
    p4coef=[]
    p3coef=[]
    p1coef=[]
    md=[]
    tau=1000.
    wrf=2.*mdot/r**2*omega*(sqrt(r)-0.9055*sqrt(rin))

    f=open('rastr.txt','w')

    while(r>=rin):
        dr=ddr*r
        r1=r+dr/2.

        tau,h=tausolve(r, omega, wrf, mdot, taurecent=tau)

        if(isnan(tau)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.)
        omega1=omega+fdomega(r, omega, tau, h, wrf)*dr/2.

        wrf1=wrf+fdwrf(r, omega, wrf,h,tau,mdot)*dr/2.

       
        if ((h/r)<1):
            rsph=r

        if((h/r)>1.):
            mdot1=mdot+fdmdot(r, wrf, tau,h)*dr/2.
            mdotss=mdotglobal*r/rsph
        else:
            mdot1=mdotglobal
            mdotss=mdotglobal

        
        if(wrf1<=0.):
            print "negative stress! wrf = "+str(wrf1)
            print "dr="+str(dr)
            print "mdot= "+str(mdot)
            print "dmdot= "+str(fdmdot(r, wrf, tau,h))
            print "tcent= "+str(b.ctemp(h,wrf,tau))
            print "r="+str(r)
            return sqrt(-1.), sqrt(-1.), sqrt(-1.)
        print r
        tau1,h1=tausolve(r1, omega1, wrf1, mdot1, taurecent=tau)
        

        if(isnan(tau1)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.)
        omega=omega+fdomega(r1, omega1, tau1, h1, wrf1)*dr

        wrf=wrf+fdwrf(r1, omega1, wrf1,h1,tau1,mdot1)*dr

        if((h1/r1)>1.):
            mdot=mdot+fdmdot(r1, wrf1, tau1,h1)*dr
        else:
            mdot=mdotglobal

             
        r=(1.+ddr)*r
        if(r<(rlast/(1.+drout))):
            rl.append(r)
            ol.append(omega)
            tl.append(tau)
            hl.append(h)
            md.append(mdotss)
            mdot123.append(mdot1)
            fcorr=1.-sqrt(rin/r)
            fcorr=(fcorr+fabs(fcorr))/2.
            temp1=9.6e7*b.ctemp(h,wrf,tau)/(1.+tvert*tau)**0.25
            tempss=5.64e7*mdot**0.25*r**(-0.75)
            tausph=0.8*mdot/(alpha*sqrt(r))
            taussa=ss.ftaussa(r,mdot, alpha, rin)
            taussb=ss.ftaussb(r,mdot, alpha, rin)
 
            hssa=ss.hssa(mdot,r, rin)
            hssb=ss.hssb(r,mdot, alpha, rin)
            f.write(str(r)+' '+str(omega)+' '+str(h/r)+' '+str(hssa/r)+' '+str(hssb/r)+' '+str(tau)+' '+str(taussa)+' '+str(taussb)+' '+str(temp1)+' '+str(tempss)+' '+str(mdot)+' '+str(wrf)+'\n')
            rlast=r
    print r
    
    print "mdot= "+str(mdot)
    print "dmdot= "+str(fdmdot(r, wrf, tau,h))
    print "tcent= "+str(b.ctemp(h,wrf,tau))
    print "r="+str(r)
    print "tau= "+str(tau)
    print "h= "+str(h)

    f.close()
    rar=asarray(rl, dtype=double)
    ssmdot=asarray(md, dtype=double)
    oar=asarray(ol, dtype=double)
    tauar=asarray(tl, dtype=double)
    har=asarray(hl, dtype=double)
    mmdot=asarray(mdot123, dtype=double)

    tausph=0.8*mdot/(alpha*sqrt(rar))
    taua=ss.ftaussa(rar,mdot, alpha, rin)
    taub=ss.ftaussb(rar,mdot, alpha, rin)
    hrssa=ss.hssa(mdot, rar, rin)/rar
    hrssb=ss.hssb(rar,mdot, alpha, rin)/rar

    omin=oar.min()
    omax=oar.max()
    if(omax<1.):
        omax=1.

    oa=ss.ost(rar, rin, mdot, alpha, taua)
    ob=ss.ost(rar, rin, mdot, alpha, taub)
    oo=ss.ost(rar, rin, mdot, alpha, tauar)

    plt.clf()
    fig=figure()
    plot(rar/rin, oar,color='r')
    plot(rar/rin, oar*0.+oin, linestyle='dashed', color='k')
    plot(rar/rin, 1.-oa, color='k', linestyle='dotted')
    plot(rar/rin, 1.-ob, color='r', linestyle='dotted')
    ylabel(r'$\omega$',fontsize=18)
    xlabel(r'$R/R_{\rm in}$',fontsize=18)
    xscale('log')
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    fig.set_size_inches(6, 5)
    ylim(omin,omax)
    savefig('Romega.eps')
    close()



    plt.clf()
    fig=figure()
    plot(rar/rin, mmdot,color='r')
    plot(rar/rin, ssmdot,color='b')
    ylabel(r'$\dot m$',fontsize=18)
    xlabel(r'$R/R_{\rm in}$',fontsize=18)
 #   xscale('log')
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    fig.set_size_inches(6, 5)
#    ylim(omin,omax)
    savefig('Rdotm.eps')
    close()

    plt.clf()
    fig=figure()
#    plot(rar/rin, tausph, color='red',label='spherical')
    plot(rar/rin, tauar, color='r',label='our')
    plot(rar/rin, taua, color='green',label='SS a zone', linestyle='dashed')
    plot(rar/rin, taub, color='blue',label='SS b zone', linestyle='dotted')
#    plt.legend(loc='upper right')
    xscale('log')
    yscale('log')
    ylabel(r'$\tau $',fontsize=18)
    xlabel(r'$R/R_{\rm in}$',fontsize=18)
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    ylim(lgmin(tauar), lgmax(tauar*5.))
    xlim(1.,100.)
    fig.set_size_inches(6, 5)
    savefig('Rtau.eps')
    close()

    plt.clf()
    fig=figure()
    plt.plot(rar/rin, har/rar, color='r')
    plt.plot(rar/rin, hrssa*hvert, color='green',label='SS a zone', linestyle='dashed')
#    plt.plot(pr1torin, phrssa, color='green',linestyle='--',label='SS a zone')
    plt.plot(rar/rin, hrssb*hvert, color='blue',label='SS b zone', linestyle='dotted')
#    plt.legend(loc='upper right')
    xscale('log')
    ylabel('$H/R$',fontsize=18)
    xlabel(r'$R/R_{\rm in}$',fontsize=18)
    xlim(1.,100.)
    ymin=lgmin(har/rar)
    ymina=lgmin(hrssa*hvert)
    yminb=lgmin(hrssb*hvert)
    if(yminb>ymina):
        ymdisc=yminb
    else:
        ymdisc=ymina
        ymdisc=ymdisc/2.
    if(ymdisc<0.01):
        ymdisc=0.01
    if(ymin>ymdisc):
        ymin=ymdisc
    ymax=lgmax(har/rar)
    ymaxa=lgmax(hrssa*hvert)
    ymaxb=lgmax(hrssb*hvert)
    if(ymaxb<ymaxa):
        ymdisc=ymaxb
    else:
        ymdisc=ymaxa
    if(ymax<ymdisc):
        ymax=ymdisc
    ylim(ymin,ymax)
    yscale('log')
    tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
    fig.set_size_inches(6, 5)
    savefig('Rhr.eps')
    close()

# main routine calculating the radial structure
def rastr(rrin, fplot,qeq):
    ra=b.rafun()

    rin=ra*rrin

    rout=20.*rin
    ddr=-5.e-4
    mdot=mdotglobal
    if(mdot<1.5):
        ddr*=(mdot/1.5)
    ddr/=(1.+0.5*(b.peq()/ps))
    defac=0.99

    r=rout
    omega=1.
    drout=1e-2
    rlast=rout

    # lists:
    rl=[]
    ol=[]
    tl=[]
    hl=[]
    temp1=[]
    tss=[]
    tau=0.

    htormax=0.
    wrf=2.*mdot/r**2*omega*(sqrt(r)-qeq*sqrt(rin))

    while(rlast>=rin):
        dr=ddr*(r-rin*defac)
        r1=r+dr/2.
    
        tau,h=tausolve(r, omega, wrf, mdot, taurecent=tau)
        if(isnan(tau)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.)
        omega1=omega+fdomega(r, omega, tau, h, wrf)*dr/2.
        wrf1=wrf+fdwrf(r, omega, wrf,h,tau,mdot)*dr/2.
        if((h/r)>1.):
            mdot1=mdot+fdmdot(r, wrf, tau,h)*dr/2.
        else:
            mdot1=mdotglobal
 
        if(wrf1<=0.):
            print "negative stress! wrf = "+str(wrf1)
            return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.),sqrt(-1.), sqrt(-1.)  
        tau1,h1=tausolve(r1, omega1, wrf1, mdot1, taurecent=tau)
        if(isnan(tau1)):
            print "tausolve resulted in NaN"
            return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.)
        omega=omega+fdomega(r1, omega1, tau1, h1, wrf1)*dr
        wrf=wrf+fdwrf(r1, omega1, wrf1,h1,tau1,mdot1)*dr
        if((h1/r1)>1.):
            mdot=mdot+fdmdot(r1, wrf1, tau1,h1)*dr
        else:
            mdot=mdotglobal
       
        r=(1.+ddr)*r
        if(r<(rlast/(1.+drout))):
            rl.append(r)
            ol.append(omega)
            tl.append(tau)
            hl.append(h)
            temp=b.ctemp(h,wrf,tau)/(1.+tvert*tau)**0.25
            tempss=5.64/9.6*mdot**0.25*r**(-0.75)
            temp1.append(temp)
            tss.append(tempss)
            rlast=r
            if((h/r)>htormax):
                htormax=h/r
                rmax=r

    mdotin=mdot
    wrfin=wrf
    rar=asarray(rl, dtype=double)
    oar=asarray(ol, dtype=double)
    tauar=asarray(tl, dtype=double)
    har=asarray(hl, dtype=double)
    temperature=asarray(temp1, dtype=double)
    temperaturess=asarray(tss, dtype=double)

    if(isnan(oar.sum())|(size(rar)<2)):
        return sqrt(-1.), sqrt(-1.), sqrt(-1.), sqrt(-1.),sqrt(-1.), sqrt(-1.) 

    oint=interp1d(rar, oar)
    hint=interp1d(rar, har)
    tint=interp1d(rar, tauar)

    return oint(rin), hint(rin), tint(rin), htormax, mdotin,wrfin
 

    
def rmesh_qeq(newmu, newmdot, newps):

    b.parset(newmu=newmu, newmdot=mdotglobal,newp=newps,neweta=0.,newalpha=0.1)
    parset(newmu=newmu, newmdot=mdotglobal,newp=newps,neweta=0.,newalpha=0.1)


    qeqmin=0.
    qeqmax=1.0
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
        xi=ordiv(newmu, newmdot, newps,0.0,qeq[k])
        print "xi= "+str(xi)
        rin=xi*b.rafun()
        oo,hh,tt,hrmax,mdotin,ww=rastr(xi,False,qeq[k])
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




def picture_qeq(newmu, newmdot, newps):

    b.parset(newmu=newmu, newmdot=mdotglobal,newp=newps,neweta=0.,newalpha=0.1)
    parset(newmu=newmu, newmdot=mdotglobal,newp=newps,neweta=0.,newalpha=0.1)

    qeqmin=0.0
    ximin=0.9
    qeqmax=1.0
    ximax=1.1
    nr1=5
    nr2=6

    qeq=(qeqmax-qeqmin)*arange(nr2)/double(nr2-1)+qeqmin
    xi=(ximax-ximin)*arange(nr1)/double(nr1-1)+ximin

    do=zeros([nr1,nr2], dtype=double)
    dw=zeros([nr1,nr2], dtype=double)

    for k in arange(nr1):
        for p in arange(nr2):
            print "picture_qeq: xi = "+str(xi[k])+", qeq="+str(qeq[p])
            do[k,p]=doffo(xi[k],qeq[p])
            dw[k,p]=doffwrfin(xi[k],qeq[p])

    plt.clf()
    contour(qeq, xi, do,levels=[0.])
    contour(qeq, xi, dw,levels=[0.])
#    colorbar()
#    plot(qeq, smdot, color='k')
    xlabel('qeq')
    ylabel(r'$\xi$')
    savefig('xq.eps')

def doffo(xi,qeq):
    ra=b.rafun()
    rin=ra*xi
    oin,hin,tauin,hrmax,mdotin,wrfin=rastr(xi,False,qeq)
    oo=b.oin(rin, hin,mdotin)
    return oo-oin

def doffwrfin(xi,qeq):
    ra=b.rafun()
    rin=ra*xi
    oin,hin,tauin,hrmax,mdotin,wrfin=rastr(xi,False,qeq)
    return wrfin-b.fwrfin(rin, hin,mdotin)

def doffall(xiqeq):
    xi=xiqeq[0]
    qeq=xiqeq[1]
    print xiqeq
    return doffo(xi,qeq), doffwrfin(xi,qeq) 

def searcon(xi1, xi2,qeq):
    x1=xi1
    x2=xi2

    tol=1e-2
    
    while(abs(x1-x2)>tol):
        x=(x1+x2)/2.
        f=doffo(x,qeq)
        if(isnan(f)):
            x2=x
        else:
            x1=x

    return x1


def searconleft(xi1, xi2,qeq):
    x1=xi1
    x2=xi2

    tol=1e-2
    
    while(abs(x1-x2)>tol):
        x=(x1+x2)/2.
        f=doffo(x,qeq)
        if(isnan(f)):
            x1=x
        else:
            x2=x

    return x2

def searex(xi1, xi2,qeq):
    x1=xi1
    x2=xi2

    f1=doffo(x1,qeq)
    f2=doffo(x2,qeq)

    tol=1e-2

    while(abs(x1-x2)>tol):
        x=(x1+x2)/2.
        f=doffo(x,qeq)
        x0=par0([x1,x,x2], [f1,f,f2])
        if(x0<x):
            x2=x
        else:
            x1=x
        print str(x1)+" "+str(x2)
        print "("+str(f1)+" "+str(f2)+")"

    return x



def searchsign_old(x1, x2, f1,tol,qeq):
    if(abs(x2-x1)<tol):
        print "too close"
        return sqrt(-1.)
    x=(x1+x2)/2.
    f=doffo(x,qeq)
    print f
    if((f*f1)<0.):
        print "1"+str(x)
        return x
    else:
        xx1=searchsign_old(x1, x, f1,tol,qeq)
        print "xx1= "+str(xx1)
        if(isfinite(xx1)):
            print "2"+str(xx1)
            return xx1
        
        xx2=searchsign_old(x, x2, f1,tol,qeq)
        print "3"+str(xx2)
        return xx2

def searchsign(x1, x2, f1, tol,xi):
    if(abs(x2-x1)<tol):
        print "too close"
        return sqrt(-1.)
    x=(x1+x2)/2.
    f=doffo(xi,x)
    print f
    if((f*f1)<0.):
        print "1"+str(x)
        return x
    else:
        xx1=searchsign(x1, x, f1,tol,xi)
        print "xx1= "+str(xx1)
        if(isfinite(xx1)):
            print "2"+str(xx1)
            return xx1
        
        xx2=searchsign(x, x2, f1,tol,xi)
        print "3"+str(xx2)
        return xx2

def wrfincheck(xi, tol,qeq):
    xi1=xi*(1.-tol)
    xi2=xi*(1.+tol)

    wrf1=doffwrfin(xi1,qeq)
    wrf2=doffwrfin(xi2,qeq)
    print "wrf1= "+str(wrf1)+" wrf2= "+str(wrf2)+'\n'

    return ((wrf2*wrf1)<0.)

def ocheck(xi, tol,qeq):
    xi1=xi*(1.-tol)
    xi2=xi*(1.+tol)

    o1=doffo(xi1,qeq)
    o2=doffo(xi2,qeq)

    return ((o2*o1)<0.)

def qeqsearch(newmu, newmdot, newps,neweta):
    tstart=time.time()
    b.parset(newmu=newmu, newmdot=newmdot,newp=newps,neweta=neweta,newalpha=0.1)
    parset(newmu=newmu, newmdot=newmdot,newp=newps,neweta=neweta,newalpha=0.1)
    
    ra=b.rafun()

    sol=scipy.optimize.root(doffall, [1.,0.4], method='broyden1')

    tend=time.time()
    print "time is "+str(tend-tstart)

    print sol
    

    return 0.


 #   scipy.optimize.fsolve()

    


def qeqdiv(newmu, newmdot, newps,neweta):
    tstart=time.time()
    b.parset(newmu=newmu, newmdot=newmdot,newp=newps,neweta=neweta,newalpha=0.1)
    parset(newmu=newmu, newmdot=newmdot,newp=newps,neweta=neweta,newalpha=0.1)
    
    ra=b.rafun()

 
 #   nui=raw_input()

    qeq1=0.
    xi=ordiv(newmu, newmdot, newps,neweta,qeq1)
    print "xi= "+str(xi)
    f1=doffo(xi,qeq1)
    print "f1= "+str(f1) 

    qeq2=1.0
    xi=ordiv(newmu, newmdot, newps,neweta,qeq2)
    print "xi= "+str(xi)
    f2=doffo(xi,qeq2)
    print "f2= "+str(f2)

    if((f1*f2)>0.):
        qeq2=searchsign(qeq1, qeq2, f1,1e-2,xi)

    tol=1e-3

    while(abs(qeq2-qeq1)>tol):
        qeq=(qeq1+qeq2)/2.
        xi=ordiv(newmu, newmdot, newps,neweta,qeq)
   
        f=doffo(xi,qeq)
        if((f*f1)>0.):
            qeq1=qeq
            f1=f
            print "qeq = "+str(qeq1)+'..'+str(qeq2)
        else:
            qeq2=qeq
            f2=f
        print "(f*f1)>0."

    if((f1*f2)>0.):
        return sqrt(-1.)
    
    # linear interpolation:
    qeq=qeq1-f1/(f2-f1)*(qeq2-qeq1)

    tend=time.time()
    print "ordiv took "+str(tend-tstart)+"s"
    print "qeq= "+str(qeq)
    
    return qeq

def ordiv(newmu, newmdot, newps,neweta,qeq):
    tstart=time.time()
    b.parset(newmu=newmu, newmdot=newmdot,newp=newps,neweta=neweta,newalpha=0.1)
    parset(newmu=newmu, newmdot=newmdot,newp=newps,neweta=neweta,newalpha=0.1)
 #   qeq=0.5
    ra=b.rafun()

    xi1=0.1
    f1=doffwrfin(xi1,qeq)
    xi2=2.0
    f2=doffwrfin(xi2,qeq)
    if(isnan(f2)):
        xi2=searcon(0.5,xi2,qeq)
        f2=doffwrfin(xi2,qeq)
        print "rdiv: upgraded xi2 = "+str(xi2)

    if(isnan(f1)):
        print "NaN at left boundary "
        xi1=searconleft(xi1,0.5,qeq)
        f1=doffwrfin(xi1,qeq)
        print "rdiv: upgraded xi1 = "+str(xi1)

    if((f1*f2)>0.):
        xi2=searchsign_old(xi1, xi2, f1,1e-2,qeq)

    tol=1e-3

    while(abs(xi2-xi1)>tol):
        xi=(xi1+xi2)/2.
        f=doffwrfin(xi,1)
        if((f*f1)>0.):
            xi1=xi
            f1=f
        else:
            xi2=xi
            f2=f

    if((f1*f2)>0.):
        return sqrt(-1.)
    
    # linear interpolation:
    xi=xi1-f1/(f2-f1)*(xi2-xi1)

    tend=time.time()
    print "ordiv took "+str(tend-tstart)+"s"
    
    return xi
