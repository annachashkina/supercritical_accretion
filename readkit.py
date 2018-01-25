#READS DIFFERENT FILES
import matplotlib
from matplotlib import rc
from matplotlib import axes
from numpy import *
from pylab import *
from scipy.integrate import *
from scipy.interpolate import *
from scipy.optimize import curve_fit

from parameters import *

#Uncomment the following if you want to use LaTeX in figures 
rc('font',**{'family':'serif','serif':['Times']})
rc('mathtext',fontset='cm')
rc('mathtext',rm='stix')
rc('text', usetex=True)
# #add amsmath to the preamble
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amssymb,amsmath}"] 

# reading period-mesh file:
def perread(fname):
    
    p=[]
    xi=[]
    o=[]
    h=[]
    t=[]

    fmine=open(fname,'r')
    s=str.split(str.strip(fmine.readline()))
    while(s):
        p.append(s[0]) # in peq units
        xi.append(s[1])
        o.append(s[2])
        h.append(s[3])
        t.append(s[4])
        s=str.split(str.strip(fmine.readline()))
        
    fmine.close()

    p=asarray(p,dtype=double)
    xi=asarray(xi,dtype=double)
    h=asarray(h,dtype=double)
    o=asarray(o,dtype=double)
    t=asarray(t,dtype=double)

    

    return p,xi,o,h,t

# reading structure file r(GM/c^2) -- omega -- H/R -- tau

def readrastr(filename):

    ff=open(filename, 'r')

    r=[]
    o=[]
    htor=[]
    t=[]
    tc=[]
    wrf=[]
    ta=[]
    tb=[]
    q1=[]
    q2=[]
    q3=[]
    ha=[]
    hb=[]

    s=str.split(str.strip(ff.readline()))
    while(s):
        r.append(s[0])
        o.append(s[1])
        htor.append(s[2])
        t.append(s[3])
        tc.append(s[4])
        wrf.append(s[5])
        ta.append(s[6])
        tb.append(s[7])
        q1.append(s[8])
        q2.append(s[9])
        q3.append(s[10])
        ha.append(s[11])
        hb.append(s[12])

        
        s=str.split(str.strip(ff.readline()))

    r=asarray(r,dtype=double)
    o=asarray(o,dtype=double)
    htor=asarray(htor,dtype=double)
    t=asarray(t,dtype=double)
    tc=asarray(tc,dtype=double)
    wrf=asarray(wrf,dtype=double)
    ta=asarray(ta,dtype=double)
    tb=asarray(tb,dtype=double)
    qadv=asarray(q1,dtype=double)
    qrad=asarray(q2,dtype=double)
    qplus=asarray(q3,dtype=double)
    ha=asarray(ha,dtype=double)
    hb=asarray(hb,dtype=double)
    htora=ha/r
    htorb=hb/r
    

    ff.close()

    return r,o,htor,t, tc, wrf, ta, tb, qadv, qrad, qplus, htora, htorb




def rastr():

    fname1='anna_100.txt'
    fname2='rastr.txt'

    r,o,htor,t, tc, wrf, ta, tb, qadv, qrad, qplus, htora, htorb=readrastr(fname1)
    r_old,o_old,htor_old,t_old, tc_old, wrf_old, ta_old, tb_old, q1_old, q2_old, q3_old, htora_old, htorb_old=readrastr(fname2)
    rin=140.


    plt.clf()
    fig=figure()
    plt.subplot (3, 3, 1)
 #   plot(r*206746., ta, color='blue',label='rad')
#    plot(r*206746., tb, color='green',label='gas')
    plot(r*206746.*rin, t, color='red',label='calc-new')
    plot(r_old*206746.*rin, t_old, color='green',label='calc-old')
    ylabel(r'$\tau$')
    xlabel('$r$')
    yscale('log')
    xscale('log')
    legend()

    plt.subplot (3, 3, 2)
    plot(r*206746.*rin, htor, color='red',label='calc-new')
    plot(r_old*206746.*rin, htor_old, color='green',label='calc-old')
#    plot(r*206746., htora, color='blue',label='rad')
#    plot(r*206746., htorb, color='green',label='gas')
    ylabel('$h/r$')
    xlabel('$r$')
    xscale('log')
 #   yscale('log')
    legend()

    plt.subplot (3, 3, 3)
    plot(r*206746.*rin, wrf*2.56787e+21, color='red',label='calc-new')
    plot(r_old*206746.*rin, wrf_old*2.56787e+21, color='green',label='calc-old')
    ylabel(r'$W_{rf}$')
    xlabel('$r$')
    xscale('log')
    yscale('log')
    legend()


    plt.subplot (3, 3, 4)
    plot(r*206746.*rin, o, color='red',label='calc-new')
    plot(r_old*206746.*rin, o_old, color='green',label='calc-old')
    ylabel(r'$\Omega/\Omega_{\rm K}$')
    xlabel('$r$')
    xscale('log')
    legend()
#    yscale('log')

    plt.subplot (3, 3, 5)
    plot(r*206746.*rin, tc*9.6e7/(r*206746.)**(3./4.), color='red',label='calc-new')
    plot(r_old*206746.*rin, tc_old*9.6e7/(r_old*206746.)**(3./4.), color='green',label='calc-old')
    ylabel(r'$T_{\rm c}R^{-3/4}$')
    xlabel('$r$')
    xscale('log')
 #   yscale('log')
    legend()
  

    plt.subplot (3, 3, 6)
    plot(r*206746.*rin, fabs(qadv/qplus), color='red',label='Qadv-Qvis')
    plot(r*206746.*rin, fabs(qrad/qplus), color='green',label='Qrad-Qvis')
    ylim(0.0001,1.01)
    ylabel('$Q$')
    xlabel('$r$')
    yscale('log')
    xscale('log')
    legend()

    fig.set_size_inches(15, 15)
    savefig('comp.eps')
    



def stread(filename):

    ff=open(filename, 'r')

    r=[]
    o=[]
    htor=[]
    t=[]

    s=str.split(str.strip(ff.readline()))
    while(s):
        r.append(s[0])
        o.append(s[1])
        htor.append(s[2])
        t.append(s[3])
        s=str.split(str.strip(ff.readline()))

    r=asarray(r,dtype=double)
    o=asarray(o,dtype=double)
    htor=asarray(htor,dtype=double)
    t=asarray(t,dtype=double)

    ff.close()

    return r,o,htor,t

def rmmread(fname):
    fin=open(fname+'.txt', 'r')

    s=str.split(str.strip(fin.readline()))
    
    mu2=[]
    md2=[]
    xi=[]
    o=[]

    while(s):
        if(size(s)>=4):
            mu2.append(s[0])
            md2.append(s[1])
            xi.append(s[2])
            o.append(s[3])
        s=str.split(str.strip(fin.readline()))
    fin.close()
    print np.size(mu2)
    print np.size(md2)
    mmu2=asarray(mu2, dtype=double)
    mmd2=asarray(md2, dtype=double)
    mxi=asarray(xi, dtype=double)
    mo=asarray(o, dtype=double)
 
    return mmu2,mmd2,mxi,mo


def per():

    fname1='rastr.txt'
    fname2='rastr_new.txt'
    
    r1=[]
    omega1=[]
    hr1=[]
    tau1=[]
    tc1=[]
    wrf1=[]

    r2=[]
    omega2=[]
    hr2=[]
    tau2=[]
    tc2=[]
    wrf2=[]

    fmine=open(fname1,'r')
    s=str.split(str.strip(fmine.readline()))
    while(s):
        r1.append(s[0]) # in peq units
        omega1.append(s[1])
        hr1.append(s[2])
        tau1.append(s[3])
        tc1.append(s[4])
        wrf1.append(s[5])
        
        s=str.split(str.strip(fmine.readline()))
        
    fmine.close()

    rold=asarray(r1,dtype=double)
    omegaold=asarray(omega1,dtype=double)
    hrold=asarray(hr1,dtype=double)
    tauold=asarray(tau1,dtype=double)
    tcold=asarray(tc1,dtype=double)
    wrfold=asarray(wrf1,dtype=double)

    fmine=open(fname2,'r')
    s=str.split(str.strip(fmine.readline()))
    while(s):
        r2.append(s[0]) # in peq units
        omega2.append(s[1])
        hr2.append(s[2])
        tau2.append(s[3])
        tc2.append(s[4])
        wrf2.append(s[5])
        s=str.split(str.strip(fmine.readline()))
        
    fmine.close()

    rnew=asarray(r2,dtype=double)
    omeganew=asarray(omega2,dtype=double)
    hrnew=asarray(hr2,dtype=double)
    taunew=asarray(tau2,dtype=double)
    tcnew=asarray(tc2,dtype=double)
    wrfnew=asarray(wrf2,dtype=double)

    plt.clf()
    plot(rold, omegaold,label='old')
    plot(rnew, omeganew,label='new')
    ylabel(r'$\omega$')
    legend()
    xlabel('$r$')
 #   yscale('log')
    xscale('log')
    savefig('comp_omega.eps')

    plt.clf()
    plot(rold, hrold,label='old')
    plot(rnew, hrnew,label='new')
    ylabel('h/r')
    xlabel('$r$')
    legend()
    yscale('log')
    xscale('log')
    savefig('comp_hr.eps')

    plt.clf()
    plot(rold, tauold,label='old')
    plot(rnew, taunew,label='new')
    ylabel(r'$\tau$')
    xlabel('$r$')
    yscale('log')
    legend()
    xscale('log')
    savefig('comp_tau.eps')

    plt.clf()
    plot(rold, tcold,label='old')
    plot(rnew, tcnew,label='new')
    ylabel('tc')
    xlabel('$r$')
    legend()
    yscale('log')
    xscale('log')
    savefig('comp_tc.eps')


    plt.clf()
    plot(rold, wrfold,label='old')
    plot(rnew, wrfnew,label='new')
    ylabel('wrf')
    xlabel('$r$')
    legend()
    yscale('log')
    xscale('log')
    savefig('comp_wrf.eps')

    return 0.


def sigma_curve():

    fname1='exp500_1000.txt'

    
    mdot1=[]
    ksi1=[]
    qeq1=[]
    sigma1=[]
    sigma2=[]

    fmine=open(fname1,'r')
    s=str.split(str.strip(fmine.readline()))
    while(s):
        mdot1.append(s[0]) # in peq units
        ksi1.append(s[1])
        qeq1.append(s[2])
        sigma1.append(s[3])
        sigma2.append(s[4])
        
        
        
        s=str.split(str.strip(fmine.readline()))
        
    fmine.close()

    mdot=asarray(mdot1,dtype=double)
    sigma500=asarray(sigma1,dtype=double)
    sigma1000=asarray(sigma2,dtype=double)


    fname1='mdotexp_factor075.txt'

    
    mdot175=[]
    ksi175=[]
    qeq175=[]
    sigma175=[]
    sigma275=[]
    jj=[]

    fmine=open(fname1,'r')
    s=str.split(str.strip(fmine.readline()))
    while(s):
        mdot175.append(s[1]) # in peq units
        ksi175.append(s[2])
        qeq175.append(s[3])
        sigma175.append(s[4])
        sigma275.append(s[5])
        jj.append(s[0])
        
        s=str.split(str.strip(fmine.readline()))

    
    fmine.close()
    mdot75=asarray(mdot175,dtype=double)
    sigma50075=asarray(sigma175,dtype=double)
    sigma100075=asarray(sigma275,dtype=double)

    print sigma100075
    jij=raw_input()


 

    plt.clf()
    plot(sigma500,mdot,color='red', label='r/Rg = 1000')
    plot(sigma1000,mdot,color='blue',label='r/Rg = 500')
    plot(sigma50075,mdot75,color='black', label='r/Rg = 1000')
    plot(sigma100075,mdot75,color='green',label='r/Rg = 500')
    xlabel(r'$\Sigma \kappa$')
    legend()
    ylabel(r'$\dot m$')
    yscale('log')
    xscale('log')
    savefig('sigma_curve.eps')

    return 0.


def rmmread_all(fname):
    fin=open(fname+'.txt', 'r')

    s=str.split(str.strip(fin.readline()))
    
    mu2=[]
    md2=[]
    xi=[]
    qeq=[]
    oint=[]
    hint=[]
    htormax=[]
    mdotin=[]
    wint=[]

    while(s):
        if(size(s)>=4):
            mu2.append(s[0])
            md2.append(s[1])
            xi.append(s[2])
            qeq.append(s[3])
            oint.append(s[4])
            hint.append(s[5])
            htormax.append(s[6])
            mdotin.append(s[7])
            wint.append(s[8])
    
        s=str.split(str.strip(fin.readline()))
    fin.close()


    mmu2=asarray(mu2, dtype=double)
    mmd2=asarray(md2, dtype=double)
    mxi=asarray(xi, dtype=double)
    mqeq=asarray(qeq, dtype=double)
    moint=asarray(oint, dtype=double)
    mhint=asarray(hint, dtype=double)
    mhtormax=asarray(htormax, dtype=double)
    mmdotin=asarray(mdotin, dtype=double)
    mwint=asarray(wint, dtype=double)

    mrin=(lam*mmu2**2/mmd2)**(2./7.)*2.**(-1./7.)*mxi

    madvterm=3.*mmdotin*mhint/(5.*mrin**2.)

    return mmu2,mmd2,mxi,mqeq,moint,mhint,mhtormax,mmdotin/mmd2,mwint,madvterm


