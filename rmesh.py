import numpy as np
import matplotlib.pyplot as plt
import numpy.random
import time
import os

import bocon as b
import dirin as d

# functions:
from physics import G, beta, comegaga, ctau, cwrf, ctemp, fh, ftau
from physics import Scal, Pcal, Qcal # SPQR
# constants
from physics import mu, ps, mdotglobal, alpha, eta, kt, epsilon, psi, n, mmean, lam, chi, tvert, hvert
# more constants:
from physics import tol, toff, varold, varnew, qeqest, xiest, defac

def rmesh_qeq(newmu, newmdot, newps):

    b.parset(newmu=newmu, newmdot=mdotglobal,newps=newps,neweta=0.,newalpha=0.1)
    d.parset(newmu=newmu, newmdot=mdotglobal,newps=newps,neweta=0.,newalpha=0.1)

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
        xi, tc=d.ordiv(newmu, newmdot, newps,0.0,qeq[k])
 #      tc=ordiv_tc(newmu, newmdot, newps,0.,qeq,xi)
        print "xi= "+str(xi)+'\n'
        print "qeq= "+str(qeq[k])+'\n'
        rin=xi*b.rafun()
        oo,hh,tt,hrmax,mdotin,ww=d.rastr(xi,qeq[k],tc)
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
    d.parset(newmu=newmu, newmdot=mdotglobal,newp=newps,neweta=0.,newalpha=0.1)


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
                oo,hh,tt,hrmax,mdotin,ww=d.rastr(xi,qeq,tc)
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
