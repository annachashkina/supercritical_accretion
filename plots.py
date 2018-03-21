from numpy import *
import matplotlib.pyplot as plt
import numpy.random
import time
import os
import bocon as b
import dirin as d
from matplotlib import *
from matplotlib.colors import BoundaryNorm

from parameters import *

# functions:
from physics import G, beta, comegaga, ctau, cwrf, ctemp, fh, ftau
from physics import Scal, Pcal, Qcal # SPQR
# constants
#from physics import  alpha, eta, kt, epsilon, psi, n, mmean, lam, chi, tvert, hvert
# more constants:
from physics import tol, toff, varold, varnew, qeqest, xiest, defac
from dirin import domega, dwrf, dtau, dtemp, dmdot
import ssdisk as ss
import readkit as rk

def lgmin(x):
    y=(log10(x)).min()
    y=floor(y)
    return 10.**y

def lgmax(x):
    y=(log10(x)).max()
    y=ceil(y)
    return 10.**y

# disc integration with a bunch of figures
def rastr1(rrin, qeq, newmu=d.mu, newmdot=mdotglobal, newps=ps):

    b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.0,newalpha=0.1)
    d.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.0,newalpha=0.1)
    
    ra=b.rafun()
    rin=ra*rrin
#    print rin
#    kok=raw_input()
    t=1
    rout=100.*rin
    ddr=-3.e-5
    mdot=d.mdotglobal ;   mu=d.mu
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

        
#        if((r*r)>(9.*tau1/(4.*64.*pi*tc**4.))):
        if(r<h):
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
        
 #       if((r*r)>(9.*tau1/(4.*64.*pi*tc**4.))):
        if(r<h1):
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
    ylabel('flux ratios')
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
    print "mdoyin= "+str(mdotin)+' \n'
    huh=raw_input()
    return rar*rin, tauar

# plots only xi and qeq
def rmmplot(fname):

    mu2,md2,xi,qeq=rk.rmmread(fname)
    print mu2
#    nkj=raw_input()
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

    levs=asarray([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4])
    norm = BoundaryNorm(boundaries=levs, ncolors=256)

    plt.clf()
    fig=plt.figure()
    plt.contourf(muar, mdar, xiar,levels=levs,cmap='jet') #,levels=rlev)
 #   title(r'$\xi$')
    c=plt.colorbar()
    plt.xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=18)
    plt.ylabel('$\dot{m}$',fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([muar.min(),muar.max()])
    plt.ylim([mdar.min(),mdar.max()])
   # tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
#    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
#    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
#    fig.set_size_inches(6, 5)
    plt.savefig(fname+'_xi.pdf')

    plt.clf()
    fig=plt.figure()
    plt.contourf(muar, mdar, qeqar,cmap='jet') #,levels=rlev)
    #title(r'$qeq$')
    c=plt.colorbar()
    plt.xlabel('$qeq$, $10^{30}$G\,cm$^{3}$',fontsize=18)
    plt.ylabel('$\dot{m}$',fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([muar.min(),muar.max()])
    plt.ylim([mdar.min(),mdar.max()])
   # tick_params(labelsize=16, length=8, width=1.0, which='major', pad=8)
#    tick_params(labelsize=16, length=4, width=1.0, which='minor', pad=8)
#    fig.tight_layout(pad=3.5,w_pad=0., h_pad=0.)
#    fig.set_size_inches(6, 5)
    plt.savefig(fname+'_qeq.pdf')

    plt.close()

# difference in xi for different eta values
def difeta():

    fname='rmm_eta_00_ps10peq_ready'
    fname1='rmm_eta0.1_p10eq_ready'

    mu2,md2,xi,qeq,oint,hint,htormax,mdotin,wint,madvterm=rk.rmmread_all(fname)

    wuu=unique(mu2)
    wdd=unique(md2)
    nuu=size(wuu)
    ndd=size(wdd)

    muar1=reshape(asarray(mu2,dtype=double),[nuu,ndd])
    mdar1=reshape(asarray(md2,dtype=double),[nuu,ndd])
    xiar1=reshape(asarray(xi,dtype=double),[nuu,ndd])

    muar_f1=muar1[::1,::50]
    mdar_f1=mdar1[::1,::50]
    xiar_f1=xiar1[::1,::50]

    mu22,md22,xi2,qeq2,oint2,hint2,htormax2,mdotin2,wint2,madvterm2=rk.rmmread_all(fname1)

    wuu=unique(mu22)
    wdd=unique(md22)
    nuu=size(wuu)
    ndd=size(wdd)

    muar12=reshape(asarray(mu22,dtype=double),[nuu,ndd])
    mdar12=reshape(asarray(md22,dtype=double),[nuu,ndd])
    xiar12=reshape(asarray(xi2,dtype=double),[nuu,ndd])

    muar_f2=muar12[::1,::50]
    mdar_f2=mdar12[::1,::50]
    xiar_f2=xiar12[::1,::50]

    deltaxi=xiar_f2-xiar_f1

    

    levs=np.arange(0.1,0.6,0.001)


    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar_f1, mdar_f1, deltaxi/xiar_f1, levels=levs ,cmap='jet') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal',ticks=[0.1,0.2,0.3,0.4,0.5,0.6,0.7])
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2 = ax.twinx()
    ax2.contourf(muar_f1, mdar_f1*2.*1.e38*0.1, deltaxi/xiar_f1, levels=levs,cmap='jet') #,levels=rlev)
    title(r'$\Delta \xi/\xi$',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
    ax.set_ylim(3.e3,1.e5)
 #   ax.set_ylim(1.e3,1.e5)
    ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim(3.e3*2.*1.e38*0.1,1.e5*2.*1.e38*0.1)
    fig.set_size_inches(10., 12.5)
    ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    plt.savefig(fname+'_deltaxi.eps')
    plt.close()

#
 
# plots an S-curve for mu=const, given RMMfile
def scurve(rmmfile, themu):

    mu2,md2,xi,qeq=rk.rmmread(rmmfile)

    wm=((mu2-themu)**2).argmin()
    wmu=where(mu2 == mu2[wm])
    
    nmu=size(wmu)
    mdar=asarray(md2[wmu]) ;   xiar=asarray(xi[wmu]) ;   qeqar=asarray(qeq[wmu])
    
    tinar=np.zeros(nmu, dtype=double)
    mdinar=np.zeros(nmu, dtype=double)
    mdoutar=mdar
    
    for k in arange(nmu):
        b.parset(newmu=mu2[wm], newmdot=mdar[k],newps=-10.,neweta=0.,newalpha=0.1)
        parset(newmu=mu2[wm], newmdot=mdar[k],newps=-10.,neweta=0.,newalpha=0.1)
        oint, hint, tint, htormax, mdotin, wint = rastr(xiar[k], qeqar[k])
        tinar[k]=tint
        mdinar[k]=mdotin

    clf()
    plot(tinar, mdoutar, color='r', label=r'$\dot{m}_{\rm out}$')
    plot(tinar, mdinar, color='k', label=r'$\dot{m}_{\rm in}$')
    xlabel(r'$\varkappa \Sigma$')
    ylabel(r'$\dot{m}$')
    legend()
    xscale('log')
    yscale('log')
    savefig('scurvein.eps')

# post-processing of a pre-calculated grid
def everything(rmmfile):
    
    mu2,md2,xi,qeq=rk.rmmread(rmmfile)

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
    

    fname=rmmfile+'fignya.txt'
    fout=open(fname, 'w')
    
    for ku in arange(nuu):
        for kd in arange(ndd):
            b.parset(newmu=muar[ku,kd], newmdot=mdar[ku,kd],newps=-10.,neweta=0.,newalpha=0.1)
            d.parset(newmu=muar[ku,kd], newmdot=mdar[ku,kd],newps=-10.,neweta=0.,newalpha=0.1)
            mu_e=muar[ku,kd]
            mdot_e=mdar[ku,kd]
            xi_e=xiar[ku,kd]
            qeq_e=qeqar[ku,kd]
            print mdotglobal
            print mu
            oint, hint, tint, htormax, mdotin, wint = d.rastr(xiar[ku,kd], qeqar[ku,kd])
            fout.write(str(mu_e)+' '+str(mdot_e)+' '+str(xi_e)+' '+str(qeq_e)+' '+str(oint)+' '+str(hint)+' '+str(htormax)+' '+str(mdotin)+' '+str(wint)+ '\n')
            fout.flush()
    fout.close()



def old_prep(rmmfile):
    
    mu2,md2,xi,qeq,oint,hint,htormax,mdotin,wint=rk.rmmread_all(rmmfile)

    wuu=unique(mu2)
    wdd=unique(md2)
    nuu=size(wuu)
    ndd=size(wdd)
    
    print nuu
    print ndd

    muar=reshape(asarray(mu2,dtype=double),[nuu,ndd])
    mdar=reshape(asarray(md2,dtype=double),[nuu,ndd])
    xiar=reshape(asarray(xi,dtype=double),[nuu,ndd])
    ointar=reshape(asarray(oint,dtype=double),[nuu,ndd])
    htormaxar=reshape(asarray(htormax,dtype=double),[nuu,ndd])

    fname='old_prep'
    fout=open(fname+'.txt', 'w')

    for ku in arange(nuu):
        for kd in arange(ndd):
            mu_e=muar[ku,kd]
            mdot_e=mdar[ku,kd]
            xi_e=xiar[ku,kd]
            o_e=ointar[ku,kd]
            h_e=htormaxar[ku,kd]
            fout.write(str(mu_e)+' '+str(mdot_e)+' '+str(xi_e)+' '+str(o_e)+' '+str(h_e)+' '+str(h_e)+ '\n')
            
    fout.close()


def rmmplot_all(fname):

    mu2,md2,xi,qeq,oint,hint,htormax,mdotin,wint,madvterm=rk.rmmread_all(fname)
    wuu=unique(mu2)
    wdd=unique(md2)
    nuu=size(wuu)
    ndd=size(wdd)
    print qeq.min()
    print qeq.max()
    print htormax.min()
    print htormax.max()
    print mdotin.min()
    print mdotin.max()
    print madvterm.min()
    print madvterm.max()
 #   lj=raw_input()
    
    print nuu
    print ndd

    muar1=reshape(asarray(mu2,dtype=double),[nuu,ndd])
    mdar1=reshape(asarray(md2,dtype=double),[nuu,ndd])
    xiar1=reshape(asarray(xi,dtype=double),[nuu,ndd])
    qeqar1=reshape(asarray(qeq,dtype=double),[nuu,ndd])

    htormaxar1=reshape(asarray(htormax,dtype=double),[nuu,ndd])
    mdotinar1=reshape(asarray(mdotin,dtype=double),[nuu,ndd])
    hintar1=reshape(asarray(hint,dtype=double),[nuu,ndd])

    madvtermar1=reshape(asarray(madvterm,dtype=double),[nuu,ndd])
    
    muar=muar1[::1,::50]
    mdar=mdar1[::1,::50]
    xiar=xiar1[::1,::50]
    qeqar=qeqar1[::1,::50]
    htormaxar=htormaxar1[::1,::50]
    mdotinar=mdotinar1[::1,::50]
    hintar=hintar1[::1,::50]
    madvtermar=madvtermar1[::1,::50]

    
    
    levs=np.arange(0.1,1.4,0.01)
    levs1=np.asarray([0.8,0.85,0.9,0.95,0.99], dtype=double)*100.
    levsxi=np.asarray([1.25,1.26,1.27], dtype=double)
 #   norm = colors.BoundaryNorm(boundaries=levs, ncolors=256)

    ksiold=2.**(-3./7.)*(1./3.*hintar)**(2./7.)
 
    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar, mdotinar*mdar, xiar,levels=levs,cmap='jet') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal',ticks=[0.1,0.3,0.5,0.7,0.9,1.1,1.3])
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2 = ax.twinx()
    ax2.contourf(muar, mdotinar*mdar*2.*1.e38*0.1, xiar,levels=levs,cmap='jet') #,levels=rlev)
    CS=ax2.contour(muar, mdotinar*mdar*2.*1.e38*0.1, xiar,levels=levsxi,colors='white')
    
    title(r'$\xi$',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
    ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    fig.set_size_inches(10., 12.5)
    ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    plt.savefig(fname+'_xi.eps')
    plt.close()

    xiarnew=xiar*mdotinar**(2./7.)
    difxi=xiar*(1.-mdotinar**(2./7.))
    
    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar, mdar*mdotinar, xiarnew,levels=levs,cmap='jet') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal',ticks=[0.1,0.3,0.5,0.7,0.9,1.1,1.3])
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2 = ax.twinx()
    ax2.contourf(muar, mdar*mdotinar*2.*1.e38*0.1, xiarnew,levels=levs,cmap='jet') #,levels=rlev)
    CS=ax2.contour(muar, mdar*mdotinar*2.*1.e38*0.1, xiarnew,levels=levsxi,colors='white')
    
    title(r'$\xi_{\rm eff}$',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
    ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    fig.set_size_inches(10., 12.5)
    ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    plt.savefig(fname+'_xi_modified.eps')
    plt.close()


    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar, mdar, difxi,cmap='jet') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal',ticks=[0.1,0.3,0.5,0.7,0.9,1.1,1.3])
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2 = ax.twinx()
    ax2.contourf(muar, mdar*2.*1.e38*0.1, difxi,cmap='jet') #,levels=rle
    
    title(r'$\Delta \xi$',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
    ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    fig.set_size_inches(10., 12.5)
    ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    plt.savefig(fname+'_xi_dif.eps')
    plt.close()



    levqeq=np.arange(0.8,1.42,0.01)


    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar, mdar, qeqar, levels=levqeq ,cmap='jet') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal')#,ticks=[0.1,0.3,0.5,0.7,0.9,1.1,1.3])
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2 = ax.twinx()
    ax2.contourf(muar, mdar*2.*1.e38*0.1, qeqar,levels=levqeq,cmap='jet') #,levels=rlev)
    title('$qeq$',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
    ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    fig.set_size_inches(10., 12.5)
    ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    plt.savefig(fname+'_qeq.eps')
    plt.close()


    levhr=np.arange(0.01,1.36,0.01)
 
    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar, mdar, htormaxar, levels=levhr , cmap='jet') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal')#,ticks=[0.1,0.3,0.5,0.7,0.9,1.1,1.3])
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2 = ax.twinx()
    ax2.contourf(muar, mdar*2.*1.e38*0.1, htormaxar, levels=levhr, cmap='jet') #,levels=rlev)
    title(r'$(H/R)_{\rm max}$',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
    ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    fig.set_size_inches(10., 12.5)
    ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    plt.savefig(fname+'_htormax.eps')
    plt.close()


    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar, mdar, hintar, cmap='jet') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal')#,ticks=[0.1,0.3,0.5,0.7,0.9,1.1,1.3])
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2 = ax.twinx()
    ax2.contourf(muar, mdar*2.*1.e38*0.1, hintar, cmap='jet') #,levels=rlev)
    title(r'$(H/R)_{\rm max}$',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
    ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    fig.set_size_inches(10., 12.5)
    ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    plt.savefig(fname+'_htor.eps')
    plt.close()

    

    levadv=np.arange(3.e-6,1.,0.01)

    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar, mdar, madvtermar, levels=levadv ,cmap='jet') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal')#,ticks=[0.1,0.3,0.5,0.7,0.9,1.1,1.3])
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2 = ax.twinx()
    ax2.contourf(muar, mdar*2.*1.e38*0.1, madvtermar, levels=levadv,cmap='jet') #,levels=rlev)
    title(r'$Q_{\rm adv}/Q_{\rm total}$',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
    ax.set_ylim(1.e3,1.e5)
    ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    fig.set_size_inches(10., 6.25)
    ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    plt.savefig(fname+'_adv.eps')
    plt.close()

    levmd=np.arange(0.7,0.99,0.005)
    levmd[-1]=1.1
    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar, mdar, mdotinar, levels=levmd,cmap='jet') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal')#,ticks=[0.1,0.3,0.5,0.7,0.9,1.1,1.3])
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2 = ax.twinx()
    ax2.contourf(muar, mdar*2.*1.e38*0.1, mdotinar, levels=levmd,cmap='jet') #,levels=rlev)
    CS=ax2.contour(muar, mdar*2.*1.e38*0.1, madvtermar*100.,levels=levs1,colors=('1', '0.75', '0.50','0.25','0'))  
    title(r'$\dot M_{\rm in}/\dot M_{\rm out}$',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
    ax.set_ylim(3.e3,1.e5)
    ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    fig.set_size_inches(10., 12.5)
    ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2.set_ylim(3.e3*2.*1.e38*0.1,1.e5*2.*1.e38*0.1)
    plt.savefig(fname+'_mdotin.eps')
    plt.close()

    plt.clf()
    plot(hintar,xiar*mdotinar**(2./7.),color='black')
    plot(hintar,xiar,color='blue')
    plot(hintar,1.113*hintar**(2./7.),color='red')
 #   plt.xscale('log')
    xlim(0.5,1.4)
    ylim(0.9,1.3)
    xlabel(r'$(H/R)_{\rm in}$')
    ylabel(r'$\xi$')
    plt.savefig(fname+'_htormax_vs_xiarpdf.eps')
    plt.close()



def rmmplot_crazy(fname):

    mu2,md2,xi,qeq=rk.rmmread(fname)
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
 #   muar2=muar[where(muar>0.9607*mdar**(1./2.))]
#    mdar2=mdar[where(muar>0.9607*mdar**(1./2.))]
#    nuu1=size(muar2)
#    ndd1=size(mdar2)
#    muarnew=reshape(asarray(muar2,dtype=double),[nuu1,ndd1])
#    mdarnew=reshape(asarray(mdar2,dtype=double),[nuu1,ndd1])


    dotp=6.31e-3*1.137**2.*mdar**(6./7.)*muar**(1./7.)
 #   dotpnew=1.38e-5*1.137**2.*mdarnew**(6./7.)*muarnew**(1./7.)
    
    pconst=pstar*lam**(3./7.)*muar**(6./7.)*xiar**(-0.25)*(1./3.)**(0.5)*mdar**(-3./7.)*2.**(3./7.)
    
    
    psc=muar-0.9607*mdar**(1./2.)
    levs=np.arange(0.1,1.4,0.01)
    levs1=np.asarray([0.8,0.85,0.9,0.95,0.99], dtype=double)*100.
    levsxi=np.asarray([1.25,1.26,1.27], dtype=double)
 #   norm = colors.BoundaryNorm(boundaries=levs, ncolors=256)

    #ksiold=2.**(-3./7.)*(1./3.*hintar)**(2./7.)

    plt.clf()
    plot(mdar,dotp)
    plt.savefig(fname+'_crazy_.eps')
    plt.close()
 
    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar, mdar, dotp, cmap='jet') #,levels=rlev)
#    cax1hat=ax.contour(muar, mdar, psc,levels=[0.],hatches=[None,'-'],color='none') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal')
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
 #   ax2 = ax.twinx()
#    ax2.contourf(muar, mdar*2.*1.e38*0.1, dotp,cmap='jet') #,levels=rlev)
 #   CS=ax2.contour(muar, mdar*2.*1.e38*0.1, dotp,colors='white')
    cax11hat=ax.contour(muar, mdar, psc,levels=[0., 100.],colors='white') #,levels=rlev)

    cax1hat=ax.contourf(muar, mdar, psc,levels=[-190.,0., 100.],colors='none',hatches=[None,'o'],extend='lower') #,levels=rlev)
    yt=ax.contour(muar, mdar, dotp,levels=[5.],colors='white',linestyles='dashed') #,levels=rlev)
    yt1=ax.contourf(muar, mdar, dotp,levels=[-100., 5.,100.],colors='none',hatches=['.',None],extend='lower') #,levels=rlev)
    yt2=ax.contour(muar, mdar, mdar,levels=[7500.],colors='black') #,levels=rlev)
    yt3=ax.contour(muar, mdar, mdar,levels=[5000.],colors='black') #,levels=rlev)
  
    title(r'$\dot p/(10^{-9}$\ s\ s$^{-1}$)',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
 #   ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
 #   ax2.set_yscale('log')
    fig.set_size_inches(10., 12.5)
 #   ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
 #   ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    plt.savefig(fname+'_crazy.eps')
    plt.close()

    plt.clf()
    fig=plt.figure()
    ax = fig.add_subplot(111)
    cax1=ax.contourf(muar, mdar, xiar, cmap='jet') #,levels=rlev)
    p0 = ax.get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], p0[3]+0.1, p0[2]-0.125, 0.025])
    cbar = plt.colorbar(cax1,cax=ax_cbar,orientation='horizontal')
    tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    ax2 = ax.twinx()
    ax2.contourf(muar, mdar*2.*1.e38*0.1, xiar,cmap='jet') #,levels=rlev)
 #   CS=ax2.contour(muar, mdar*2.*1.e38*0.1, dotp,colors='white')
    
    title(r'$\xi$',fontsize=20)
    ax.set_xlabel('$\mu$, $10^{30}$G\,cm$^{3}$',fontsize=20)
    ax.set_ylabel('$\dot{m}$',fontsize=20)
    ax2.set_ylabel('$L$, erg\,s$^{-1}$',fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_yscale('log')
    fig.set_size_inches(10., 12.5)
    ax2.tick_params(labelsize=18, length=8, width=1.0, which='major', pad=8)
    ax2.tick_params(labelsize=18, length=4, width=1.0, which='minor', pad=8)
    plt.savefig(fname+'_xi.eps')
    plt.close()




    

    
    
    
    

