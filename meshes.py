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

def xiqeqplane():
    newmu=1. ; newmdot=10. ; newps=-10.
    b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.0,newalpha=0.1)
    d.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.0,newalpha=0.1)
    xi1=0.5 ; xi2=1.5 ; nxi=17
    xi=(xi2-xi1)*np.arange(nxi)/double(nxi-1)+xi1
    qeq1=0.2 ;  qeq2=1.2 ; nq=21
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
    clf()
    contourf(x2, q2, np.log(aar**2+bar**2),nlevels=20)
    colorbar()
    contour(x2, q2, aar, levels=[0.], colors='w')
    contour(x2, q2, bar, levels=[0.], colors='k')
    savefig('xiqeqmap.eps')
    
# varying rout test:
def varrout():
    newmu=1. ; newmdot=10.
    d.XQset(0.5, 0.99)

    rmin=2. ; rmax=2000. ; nrads=5
    rr=(rmax/rmin)**(np.arange(nrads)/double(nrads-1))*rmin
    xiar=np.zeros(nrads, dtype=double)
    qeqar=np.zeros(nrads, dtype=double)
    fout=open('varrout.txt', 'w')
    for k in np.arange(nrads):
        thp=d.ordiv_smart(newmu, newmdot, -10., routscale=rr[k])
        xiar[k]=thp[0] ; qeqar[k]=thp[1]
        if(thp[2]):
            d.XQset(xiar[k], qeqar[k])
            print "xi = "+str(xiar[k])
            fout.write(str(rr[k])+' '+str(xiar[k])+' '+str(qeqar[k])+'\n')
            fout.flush()
    fout.close()
    clf()
    plot(rr, xiar, '.k')
    xscale('log')
    xlabel(r'$R_{\rm in}/R_{\rm out}$')
    ylabel(r'$\xi$')
    savefig('xrouttest.eps')
    
# varying ps
def varps(newmu, newmdot):
    tstart=time.time()
    d.XQset(0.5, 0.99)
    ps1=b.peq()*10. ; psfac=0.9 ; newps=ps1
    converged=True

    psar=[] ; xar=[] ; qeqar=[]
    f=open('varps.txt','w')
    
    while(converged):
        b.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
        d.parset(newmu=newmu, newmdot=newmdot,newps=newps,neweta=0.,newalpha=0.1)
        x0, x1, converged=d.ordiv_smart(newmu, newmdot, newps)
        d.XQset(x0, x1)
        print "varps: P = "+str(newps)+": "+str(x0)+", "+str(x1)
        print str(xiest)+", "+str(qeqest)
      
        if(converged):
            psar.append(newps) ; xar.append(x0) ; qeqar.append(x1) 
        f.write(str(newps)+' '+str(xiest)+' '+str(qeqest)+'\n')
        newps*=psfac
        
    f.close()
    tend=time.time()
        
# a mu=const, ps=const mesh with variable mdot (former "experimental")
def varmdot(newmu,newps): 

    tstart=time.time()

    d.XQset(0.4, 0.99)
    newmdot1=10.
    newmdot2=10000.
    nmdot=16
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
        f.write(str(newmu)+' '+str(newmdot)+' '+str(xiest)+' '+str(qeqest)+'\n')
        print str(newmu)+' '+str(newmdot)+' '+str(xiest)+' '+str(qeqest)+'\n'
        
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
    
def rmm(newps, neweta):

    d.XQset(1.0, 0.9)

    tstart=time.time()

    mumin=1. ;   mumax=100.
    mdmin=10. ;    mdmax=10000.
    nmu=15 ;    nd=16

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
            xi[ku,kd],qeq[ku,kd],conv[ku,kd] =d.ordiv_smart(muar[ku], mdar[kd], -newps)
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

