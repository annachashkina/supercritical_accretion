import matplotlib
from matplotlib import rc
from matplotlib import axes
from numpy import *
import time

from parameters import *

import bocon as b
import dirin as d

# functions:
from physics import G, beta, comegaga, ctau, cwrf, ctemp, fh, ftau, xiinf
from physics import Scal, Pcal, Qcal # SPQR
# from physics import mu, ps, mdotglobal, alpha, eta, kt, epsilon, psi, n, mmean, lam, chi, tvert, hvert 
from dirin import domega, dwrf, dtau, dtemp, dmdot
import ssdisk as ss
import readkit as rk

def varpsi():
    newmu=1. ; newmdot=10. ; newps=-10.
    b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.0,newalpha=0.1)
    d.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.0,newalpha=0.1)
    d.XQset(0.55, 5.)
    psi1=1.45 ; psidiff=0.04 ; psi2=50.
    d.psiset(psi1)

    fout=open('psivar.dat', 'w+')
    xiar=[] ; psiar=[] ;  qeqar=[]
    converged=True
    while(converged&(psi1<psi2)):
        thp=d.ordiv_smart(newmu, newmdot, newps)
        xi=thp[0] ; qeq=thp[1]
        converged=thp[2]
        if(converged):
            d.XQset(xi, qeq)
            xiar.append(xi) ; psiar.append(psi1) ; qeqar.append(qeq)
            print "xi = "+str(xi)
            fout.write(str(psi1)+' '+str(xi)+' '+str(qeq)+'\n')
            fout.flush()
            psi1+=psidiff
            d.psiset(psi1)
    fout.close()
    xiar=asarray(xiar, dtype=double) ; qeqar=asarray(qeqar, dtype=double) ; psiar=asarray(psiar, dtype=double)
    clf()
    subplot(2,1,1)
    plot(psiar, xiar, '.k')
    xscale('log')
    #    xlabel(r'$\psi$')
    ylabel(r'$\xi$')
    subplot(2,1,2)
    plot(psiar, qeqar, '.k')
    xscale('log')
    xlabel(r'$\psi$')
    ylabel(r'$q$')
    savefig('psivar.eps')
    
    
# track the solution on the xi-qeq plane
def xiqeqplane():
    newmu=1. ; newmdot=10. ; newps=-10.
    b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.0,newalpha=0.1)
    d.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.0,newalpha=0.1)
    xi1=0.2 ; xi2=1.2 ; nxi=11
    xi=(xi2-xi1)*np.arange(nxi)/double(nxi-1)+xi1
    qeq1=0.2 ;  qeq2=1.2 ; nq=10
    qeq=(qeq2-qeq1)*np.arange(nq)/double(nq-1)+qeq1
    aar=np.zeros([nxi, nq], dtype=double)
    bar=np.zeros([nxi, nq], dtype=double)
    x2=np.zeros([nxi, nq], dtype=double)
    q2=np.zeros([nxi, nq], dtype=double)
    fout=open('xqout.txt', 'w')
    for kx in np.arange(nxi):
        for kq in np.arange(nq):
            at,bt=d.vrapper([xi[kx], qeq[kq]])
            aar[kx,kq]=at ; bar[kx,kq]=bt
            x2[kx,kq]=xi[kx] ; q2[kx,kq]=qeq[kq]
            fout.write(str(xi[kx])+' '+str(qeq[kq])+' '+str(at)+' '+str(bt)+'\n')
            fout.flush()
    fout.close()
    xi0,qeq0, conv = d.ordiv_smart(newmu, newmdot, newps)
    clf()
    contourf(x2, q2, np.log(aar**2+bar**2),nlevels=20)
    colorbar()
    contour(x2, q2, aar, levels=[0.], colors='w')
    contour(x2, q2, bar, levels=[0.], colors='k')
    if(conv):
        plot(xi0, qeq0, 'or')
    savefig('xiqeqmap.eps')
    
# varying rout test:
def varrout():
    newmu=1. ; newmdot=1e5
    d.XQset(1.34, 2.69)
    xi0=1.34 ; qeq0=2.69
    
    rmin=1e3 ; rmax=1.e4 ; nrads=10
    rr=(rmax/rmin)**(np.arange(nrads)/double(nrads-1))*rmin
    xiar=np.zeros(nrads, dtype=double)
    qeqar=np.zeros(nrads, dtype=double)
    qeqarest=(sqrt(rr)-sqrt(rr[0]))/2.+qeq0
    print qeqarest
    r=raw_input("ew")
    fout=open('varrout.txt', 'w')
    for k in np.arange(nrads):
        d.XQset(xi0, qeqarest[k])
        thp=d.ordiv_smart(newmu, newmdot, -10., routscale=rr[k])
        xiar[k]=thp[0] ; qeqar[k]=thp[1]
        if(thp[2]):
            #            d.XQset(xiar[k], qeqar[k])
            print "xi = "+str(xiar[k])
            fout.write(str(rr[k])+' '+str(xiar[k])+' '+str(qeqar[k])+'\n')
            fout.flush()
    fout.close()
    clf()
    subplot(121)
    plot(rr, xiar, '.k')
    xscale('log')
#    xlabel(r'$R_{\rm in}/R_{\rm out}$')
    ylabel(r'$\xi$')
    subplot(121)
    plot(rr, qeqar, '.k')
    xscale('log')
    xlabel(r'$R_{\rm in}/R_{\rm out}$')
    ylabel(r'$qeq$')
    savefig('xrouttest.eps')
    
# varying ps
def varps(newmu, newmdot):
    tstart=time.time()
    d.XQset(1.0, 0.9)
    ps1=b.peq()*10. ; psfac=0.95 ; newps=ps1
    converged=True

    psar=[] ; xar=[] ; qeqar=[]
    f=open('varps.txt','w')
    
    while(converged):
        rscale=newmdot
        if(newmdot<100.):
            rscale=100.
        b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
        d.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
        x0, x1, converged=d.ordiv_smart(newmu, newmdot, newps, routscale=rscale)
        d.XQset(x0, x1)
        print "varps: P = "+str(newps)+": "+str(x0)+", "+str(x1)
        print "(Peq = "+str(b.peq())+")"
        #        print str(xiest)+", "+str(qeqest)
      
        if(converged):
            psar.append(newps) ; xar.append(x0) ; qeqar.append(x1) 
            f.write(str(newps)+' '+str(x0)+' '+str(x1)+'\n')
            f.flush()
        newps*=psfac
        
    f.close()
    tend=time.time()
        
# a mu=const, ps=const mesh with variable mdot (former "experimental")
def varmdot(newmu,newps): 

    tstart=time.time()

    d.XQset(1.05, 0.965)
    newmdot1=1000.
    newmdot2=10000.
    nmdot=25
    mdar=(newmdot2/newmdot1)**(arange(nmdot)/double(nmdot-1))*newmdot1
   
    xx=xiest
    qq=qeqest
    f=open('varmdot.txt','w')
    for i in arange(nmdot):
        newmdot=mdar[i]
        b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
        d.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
        x0, x1, converged=d.ordiv_smart(newmu, newmdot, newps)
        d.XQset(x0, x1)
        f.write(str(newmu)+' '+str(newmdot)+' '+str(x0)+' '+str(x1)+'\n')
        print str(newmu)+' '+str(newmdot)+' '+str(x0)+' '+str(x1)+'\n'
        f.flush()
    f.close()
    tend=time.time()
    
# a mdot=const, ps=const bisection search for magnetic moment
def musearch(newmdot, newps): 

    tstart=time.time()

    d.XQset(1.0, 0.99)
    f=open('musearch.txt','w')
    newmu1=1. ; newmu2=100.

    # right boundary, probably no solution here
#    b.parset(newmu=newmu2, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
#    print "after bparset mu="+str(mu)
#    d.parset(newmu=newmu2, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
#    print "after parset mu="+str(mu)
    x0, x1, converged=d.ordiv_smart(newmu2, newmdot, newps)
    if(converged):
        print "musearch converged, unexpectedly;  increase your magnetic field"
        print "mu = "+str(newmu2)
        print "xi = "+str(x0)+", qeq = "+str(x1)
        ii=raw_input("?")

    xi=[] ; muu=[]
    while((newmu2/newmu1-1.)>0.1):
        newmu=sqrt(newmu1*newmu2)
        b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
        d.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
        x0, x1,converged=d.ordiv_smart(newmu, newmdot, newps)
        if(converged):
            d.XQset(x0, x1)
            xiest=d.xiest ; qeqest=d.qeqest
            newmu1=newmu
            f.write(str(newmu)+' '+str(xiest)+' '+str(qeqest)+'\n')
            print str(newmu)+' '+str(xiest)+' '+str(qeqest)+'\n'
            xi.append(xiest)  ; muu.append(newmu)
        else:
            newmu2=newmu
            print str(newmu)+' '+" not converged"
    f.close()
    tend=time.time()

    clf()
    plot(muu, xi, '.k')
    xlabel(r'$\mu$, $10^{30}$G')
    ylabel(r'$\xi$')
    savefig('musearch.eps')
    print "musearch took "+str(tend-tstart)
    return muu[-1]
    
# two-dimensional mesh computed on top of existing, with different ps and/or eta
def rmmclone(inspirefile, newps, neweta):

    tstart=time.time()
    #    muar=[] ; mdar=[] ; xi0ar=[] ; qeq0ar=[]; xiar=[] ; qeqar=[]
    
    muar, mdar, xi0ar, qeq0ar = rk.rmmread(inspirefile)
    nn=np.size(muar)
#    xiar=np.zeros(nn, dtype=double) ;   qeqar=np.zeros(nn, dtype=double)
    fout=open(inspirefile+'_ps'+str(newps)+'_eta'+str(neweta)+'.txt', 'w')
    for kk in np.arange(nn):
        d.XQset(xi0ar[kk], qeq0ar[kk])
        print "xi estimate = "+str(xi0ar[kk])
        rscale=mdar[kk] # outer radius is the maximum of mdot and 100
        if(rscale<100.):
            rscale=100.
        xiar,qeqar,conv =d.ordiv_smart(muar[kk], mdar[kk], -newps, neweta=neweta, routscale=rscale)
        fout.write(str(muar[kk])+' '+str(mdar[kk])+' '+str(xiar)+' '+str(qeqar)+'\n ')
        print str(muar[kk])+' '+str(mdar[kk])+' '+str(xiar)+' '+str(qeqar)+'\n '
        print "Delta xi = "+str(xiar-xi0ar[kk])
        fout.flush()
    fout.close()
    tend=time.time()
    print "RMMclone: the grid took "+str(tend-tstart)+"s = "+str((tend-tstart)/3600.)+"h"

# two-dimensional mesh for fixed Ps (or Ps/Peq, if newps<0) and eta
def rmm(newps, neweta):

    d.XQset(1., 1.1)

    tstart=time.time()

    mumin=1. ;   mumax=100.
    mdmin=1.e4 ;    mdmax=1.e5
    nmu=17 ;    nd=5

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
    conv=np.zeros([nmu, nd], dtype=bool)
    
    fname='rmm_op'+str(newps)+'_eta'+str(neweta)
    fout=open(fname+'.dat', 'w')
    
    for ku in arange(nmu):
        if(ku>0):
            d.XQset(xi[ku-1,0], qeq[ku-1,0]) # the closest value, calculated for the previous row, same column
        for kd in arange(nd):
            mu2[ku,kd]=muar[ku]
            md2[ku,kd]=mdar[kd]
            rscale=mdar[kd] # outer radius is the maximum of mdot and 100
            if(rscale<100.):
                rscale=100.
            xi[ku,kd],qeq[ku,kd],conv[ku,kd] =d.ordiv_smart(muar[ku], mdar[kd], -newps, neweta=neweta,  routscale=rscale)
            if(conv[ku,kd]):
                d.XQset(xi[ku,kd], qeq[ku,kd]) # 
            if(not(conv[ku,kd])):
                xi[ku,kd]=sqrt(-1)
                qeq[ku,kd]=sqrt(-1)
            print "xi (mu="+str(muar[ku])+", mdot="+str(mdar[kd])+") = "+str(xi[ku,kd])
            #            oin,hin,tin,bmax,mdotint,ll=d.rastr(xi[ku,kd], qeq[ku,kd])
#            oar[ku,kd]=oin
#            har[ku,kd]=hin
#            tauar[ku,kd]=tin
#            betar[ku,kd]=bmax
#            mddotint[ku,kd]=mdotint
#            raam[ku,kd]=raffugen(muar[ku], mdar[kd])
            fout.write(str(muar[ku])+' '+str(mdar[kd])+' '+str(xi[ku,kd])+' '+str(qeq[ku,kd])+'\n ')
                       #str(oin)+' '+str(hin/raam[ku,kd]/xi[ku,kd])+' '+str(tin)+' '+str(bmax)+' '+str(mddotint[ku,kd])+'\n')
            print str(muar[ku])+' '+str(mdar[kd])+' '+str(xi[ku,kd])+' '+str(qeq[ku,kd])+'\n '
            fout.flush()
    fout.close()
    tend=time.time()
    print "RMM: the grid took "+str(tend-tstart)+"s = "+str((tend-tstart)/3600.)+"h"

