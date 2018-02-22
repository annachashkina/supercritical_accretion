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


def rmmtransf(fname):

    fname_low='sample'
    mu_l,md_l,xi_l,qeq_l,oint_l,hint_l,htormax_l,mdotin_l,wint_l,adv_l=rmmread_all(fname_low)
    mu_h,md_h,xi_h,qeq_h,oint_h,hint_h,htormax_h,mdotin_h,wint_h,adv_h=rmmread_all(fname)

 #   print mdotin_l
#    cmlkmsl=raw_input()

    # mu for everybody
    muu=[]
    mu0=1.1
    while(mu0<100):
        muu.append(mu0)
        mu0=mu0+0.3
    mu=asarray(muu, dtype=double)
 
    print 'here1'
    # \dot m for everybody
    mdu=[]
    md0=1.1
    while(md0<1.e5):
        mdu.append(md0)
        md0=md0+10.
    md=asarray(mdu, dtype=double)
    print 'here2'

    mdu1=[]
    i=0
    md0=1.
    while(md0<1000):
        mdu1.append(md0)
        i+=1
        md0=md0+10.
    md_new_l=asarray(mdu1, dtype=double)
   
    print 'here3'

    wuu=unique(mu)
    nuu=size(wuu)
    wdd=unique(md)
    ndd=size(wdd)
    
    wuu_l=unique(mu_l)
    nuu_l=size(wuu_l)
    wdd_l=unique(md_l)
    ndd_l=size(wdd_l)

    wuu_h=unique(mu_h)
    nuu_h=size(wuu_h)
    wdd_h=unique(md_h)
    ndd_h=size(wdd_h)


    levs=asarray([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4])
   

    newxi_l=scipy.interpolate.interp2d(log10(wuu_l), log10(wdd_l), xi_l, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newqeq_l=scipy.interpolate.interp2d(log10(wuu_l), log10(wdd_l), qeq_l, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newoint_l=scipy.interpolate.interp2d(log10(wuu_l), log10(wdd_l), oint_l, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newhint_l=scipy.interpolate.interp2d(log10(wuu_l), log10(wdd_l), hint_l, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newhtormax_l=scipy.interpolate.interp2d(log10(wuu_l), log10(wdd_l), htormax_l, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newmdotin_l=scipy.interpolate.interp2d(log10(wuu_l), log10(wdd_l), mdotin_l, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newwint_l=scipy.interpolate.interp2d(log10(wuu_l), log10(wdd_l), wint_l, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newadv_l=scipy.interpolate.interp2d(log10(wuu_l), log10(wdd_l), adv_l, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
   
    newxi_h=scipy.interpolate.interp2d(log10(wuu_h), log10(wdd_h), xi_h, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newqeq_h=scipy.interpolate.interp2d(log10(wuu_h), log10(wdd_h), qeq_h, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newoint_h=scipy.interpolate.interp2d(log10(wuu_h), log10(wdd_h), oint_h, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newhint_h=scipy.interpolate.interp2d(log10(wuu_h), log10(wdd_h), hint_h, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newhtormax_h=scipy.interpolate.interp2d(log10(wuu_h), log10(wdd_h), htormax_h, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newmdotin_h=scipy.interpolate.interp2d(log10(wuu_h), log10(wdd_h), mdotin_h, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newwint_h=scipy.interpolate.interp2d(log10(wuu_h), log10(wdd_h), wint_h, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
    newadv_h=scipy.interpolate.interp2d(log10(wuu_h), log10(wdd_h), adv_h, kind='cubic', copy=True, bounds_error=False, fill_value=nan)
   



#    print "3162 "+str(newxi_h(1.,3162.))+'\n'
#    print "4216.96503429 "+str(newxi_h(1.,4216.96503429))+'\n'
#    print "4300 "+str(newxi_h(1.,4300))+'\n'
#    print "4600 "+str(newxi_h(1.,4600))+'\n'
#    print "5623.4132519"+str(newxi_h(1.,5623.4132519))+'\n'
#    print "7499"+str(newxi_h(1.,7499.))+'\n'
 #   rr=raw_input('f')



    fname='rmm_ready'
    fout=open(fname+'.txt', 'w')
    ku=0
    kd=0
    
    for ku in arange(nuu):
        
        for kd in arange(ndd):

            mu_e=mu[ku]
            mdot_e=md[kd]
            if(mdot_e<=1000.):
                xi_e=newxi_l(log10(mu_e),log10(mdot_e))
                qeq_e=newqeq_l(log10(mu_e),log10(mdot_e))
                oint_e=newoint_l(log10(mu_e),log10(mdot_e))
                hint_e=newhint_l(log10(mu_e),log10(mdot_e))
                htormax_e=newhtormax_l(log10(mu_e),log10(mdot_e))
                mdotin_e=newmdotin_l(log10(mu_e),log10(mdot_e))
                wint_e=newwint_l(log10(mu_e),log10(mdot_e))
                adv_e=newadv_l(log10(mu_e),log10(mdot_e))
                
            else:
                xi_e=newxi_h(log10(mu_e),log10(mdot_e))
                qeq_e=newqeq_h(log10(mu_e),log10(mdot_e))
                oint_e=newoint_h(log10(mu_e),log10(mdot_e))
                hint_e=newhint_h(log10(mu_e),log10(mdot_e))
                htormax_e=newhtormax_h(log10(mu_e),log10(mdot_e))
                mdotin_e=newmdotin_h(log10(mu_e),log10(mdot_e))
                wint_e=newwint_h(log10(mu_e),log10(mdot_e))
                adv_e=newadv_h(log10(mu_e),log10(mdot_e))

            fout.write(str(mu_e)+' '+str(mdot_e)+' '+str(xi_e[0])+' '+str(qeq_e[0])+ ' '+ str(oint_e[0])+' '+str(hint_e[0])+' '+str(htormax_e[0])+' '+str(mdotin_e[0])+ ' ' +str(wint_e[0]) +' '+str(adv_e[0])+'\n')
        
    fout.close()

    
    

def rmmcrazy(fname):

    mu_l,md_l,xi_l,qeq_l=rmmread(fname)
    xi_l=xi_l.ravel()
    xi_l=xi_l[xi_l>1.e-5]
    print xi_l

    
 
    muu=[]
    mu0=1.1
    while(mu0<100.):
        muu.append(mu0)
        mu0=mu0+0.3
    mu=asarray(muu, dtype=double)
 
    print 'here1'
    # \dot m for everybody
    mdu=[]
    md0=500.
    while(md0<1.e4):
        mdu.append(md0)
        md0=md0+10.
    md=asarray(mdu, dtype=double)
    print 'here2'

    
    print 'here3'

    wuu=unique(mu)
    nuu=size(wuu)
    wdd=unique(md)
    ndd=size(wdd)
    
    wuu_l=unique(mu_l)
    nuu_l=size(wuu_l)
    wdd_l=unique(md_l)
    ndd_l=size(wdd_l)

    f = SmoothBivariateSpline(log10(wuu_l),log10(wdd_l),xi_l,kx=1,ky=1)

    print f(log10(505.),log10(10.))

    nhu=raw_input()


    levs=asarray([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4])
   

 #   newxi_l=scipy.interpolate.interp2d(log10(wuu_l), log10(wdd_l), xi_l, kind='linear', copy=True, bounds_error=False, fill_value=nan)
    newqeq_l=scipy.interpolate.interp2d(log10(wuu_l), log10(wdd_l), qeq_l, kind='linear', copy=True, bounds_error=False, fill_value=nan)
 
    


    fname='rmm_ready_crazy'
    fout=open(fname+'.txt', 'w')
    ku=0
    kd=0
    
    for ku in arange(nuu):
        
        for kd in arange(ndd):

            mu_e=mu[ku]
            mdot_e=md[kd]
           
            xi_e=newxi_l(log10(mu_e),log10(mdot_e))
            qeq_e=newqeq_l(log10(mu_e),log10(mdot_e))
 
    

            fout.write(str(mu_e)+' '+str(mdot_e)+' '+str(xi_e[0])+' '+str(qeq_e[0])+'\n')
        
    fout.close()

  



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
    adv=[]

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
            adv.append(s[9])
    
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
    madvterm=asarray(adv, dtype=double)

    
    wmdotin1=where(mmdotin>1.)
    if(size(wmdotin1)>0):
        mmdotin[wmdotin1]=1.
    

    mrin=(lam*mmu2**2/mmd2)**(2./7.)*2.**(-1./7.)*mxi

    madvterm=3.*mmdotin*mhint/(5.*mrin**2.)
    madvterm=madvterm/(1.+madvterm)

    return mmu2,mmd2,mxi,mqeq,moint,mhint/mrin,mhtormax,mmdotin,mwint,madvterm


#    print np.arange(1, 10)
 #   print newxi_h(1., np.arange(k0, k1, n:101))


