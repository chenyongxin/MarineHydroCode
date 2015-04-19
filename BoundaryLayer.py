# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 10:29:23 2015

@author: yc11e14
"""





import numpy

def pohlF(eta): return 2*eta-2*eta**3+eta**4
def pohlG(eta): return eta/6*(1-eta)**3

import matplotlib.pyplot


def pohlPlot(lam):
    pyplot.xlabel(r'$u/u_e$', fontsize=16)
    pyplot.axis([-0.1,1.1,0,1])
    pyplot.ylabel(r'$y/\delta$', fontsize=16)
    eta = numpy.linspace(0.0,1.0,100)
    pyplot.plot(pohlF(eta),eta, lw=1, c='k', label=r'$P_F$')
    pyplot.plot(pohlF(eta)+lam*pohlG(eta),eta, lw=2, c='g', label=r'$P_F+\lambda P_G$')
    pyplot.legend(loc='upper left')
def disp_ratio(lam): return 3./10.-lam/120.
def mom_ratio(lam): return 37./315.-lam/945.-lam**2/9072.
def df_0(lam): return 2+lam/6.
    
#pyplot.xlabel(r'$\lambda$', fontsize=16)
#lam = numpy.linspace(-12,12,100)
#pyplot.plot(lam,disp_ratio(lam), lw=2, label=r'$\delta_1/\delta$')
#pyplot.plot(lam,mom_ratio(lam), lw=2, label=r'$\delta_2/\delta$')
#pyplot.plot(lam,df_0(lam)/10., lw=2, label=r'$c_f Re_\delta/20$')
#pyplot.legend(loc='upper right')

def g_1(lam): return df_0(lam)-lam*(disp_ratio(lam)+2*mom_ratio(lam))

from scipy.optimize import bisect
lam0 = bisect(g_1,-12,12)         # use bisect method to find root between -12...12
#print 'lambda_0 = ',lam0

def ddx_delta(Re_d,lam):
    if Re_d==0: return 0                     # Stagnation point condition
    return g_1(lam)/mom_ratio(lam)/Re_d      # delta'
    




#pyplot.xlabel(r'$\lambda$', fontsize=16)
#pyplot.ylabel(r'$g_1/F$', fontsize=16)
#pyplot.plot(lam,ddx_delta(1,lam), lw=2)
#pyplot.scatter(lam0,0, s=100, c='r')
#pyplot.text(lam0,3, r'$\lambda_0$',fontsize=15)





def heun(g,psi_i,i,dx,*args):
    g_i = g(psi_i,i,*args)                      # integrand at i
    tilde_psi = psi_i+g_i*dx                    # predicted estimate at i+1
    g_i_1 = g(tilde_psi,i+1,*args)              # integrand at i+1
    return psi_i+0.5*(g_i+g_i_1)*dx             # corrected estimate



def g_pohl(delta_i,i,u_e,du_e,nu):
    Re_d = delta_i*u_e[i]/nu            # compute local Reynolds number
    lam = delta_i**2*du_e[i]/nu         # compute local lambda 
    return ddx_delta(Re_d,lam)          # get derivative


#N = 20                                      # number of steps
#x = numpy.linspace(0,numpy.pi,N)            # set up x array from 0..pi
#psi = numpy.full_like(x,1.)                 # psi array with phi0=1
def g_test(psi,i): return psi               # define derivative function
#for i in range(N-1):                        # march!
#    psi[i+1] = heun(g_test,psi[i],i,(x[i+1]-x[i]))

def march(x,u_e,du_e,nu):
    delta0 = numpy.sqrt(lam0*nu/du_e[0])                # set delta0
    delta = numpy.full_like(x,delta0)                   # delta array
    lam = numpy.full_like(x,lam0)                       # lambda array
    for i in range(len(x)-1):                           # march!
        delta[i+1] = heun(g_pohl,delta[i],i,x[i+1]-x[i],    # integrate BL using...
                          u_e,du_e,nu)                          # additional arguments
        lam[i+1] = delta[i+1]**2*du_e[i+1]/nu               # compute lambda
        if abs(lam[i+1])>12: break                          # check stop condition
    return delta,lam,i                                  # return with separation index
    
    




#nu = 1e-4                                   # viscosity
#N = 32                                      # number of steps
#s = numpy.linspace(0,numpy.pi,N)            # distance goes from 0..pi
#u_e = 2.*numpy.sin(s)                       # velocity
#du_e = 2.*numpy.cos(s)                      # gradient
#delta,lam,iSep = march(s,u_e,du_e,nu)       # solve!
#
#
#
#
#
#pyplot.ylabel(r'$\delta/R$', fontsize=16)
#pyplot.xlabel(r'$s/R$', fontsize=16)
#pyplot.plot(s[:iSep+1],delta[:iSep+1],lw=2,label='Circle')
#pyplot.plot(s,s*5/numpy.sqrt(s/nu),lw=2,label='Flat plate')
#pyplot.legend(loc='upper left')
#pyplot.scatter(s[iSep],delta[iSep], s=100, c='r')
#pyplot.text(s[iSep]+0.1,delta[iSep],'separation between\n'
#            +'%.2f' % s[iSep]+'<s<'+'%.2f' % s[iSep+1],fontsize=12)
