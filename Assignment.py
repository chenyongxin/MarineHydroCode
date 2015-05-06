# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 23:28:00 2015

@author: ChenY
"""

import numpy
from matplotlib import pyplot
from VortexPanel import Panel,plot_flow, flow_velocity
from LiftBody import solve_gamma_kutta, make_jukowski
# 1. A circular cylinder
# 1.1 N = 32

def new_make_circle(N):
    # define the end-points of the panels
    x_ends = numpy.cos(numpy.linspace(numpy.pi, -numpy.pi, N+1))
    y_ends = numpy.sin(numpy.linspace(numpy.pi, -numpy.pi, N+1))
    # define the panels
    circle = numpy.empty(N, dtype=object)
    for i in xrange(N):
        circle[i] = Panel(x_ends[i], y_ends[i], x_ends[i+1], y_ends[i+1])
    return circle

#circle = new_make_circle(64)
#solve_gamma_kutta(circle)

# get the exact potential flow solution. C_P against theta
def potential_flow_solution():
    theta = numpy.linspace(0,numpy.pi,180)
    C_P = 1-4*numpy.sin(theta)**2
    pyplot.figure(figsize=(8,6))
    pyplot.ylabel("$C_P$",fontsize=16)
    pyplot.xlabel(r'$\theta$',fontsize=16)
    pyplot.plot(theta,C_P,lw=2, c='b', label='$C_P$')
    pyplot.legend(loc='lower right')

#first of all, define the coordinate of circle, then use these coordinates to find the velocity.
#directly use flow_velocity to find the velocity along the surface

def pressure(panels):
    x = numpy.array([p.xc for p in panels])
    y = numpy.array([p.yc for p in panels])
    u,v = flow_velocity(panels,x,y)    
    return 1-u**2-v**2 
#    return 2*(numpy.multiply([1-u**2-v**2],[p.S for p in panels]))
#print 'C_P = ', pressure(circle)

# define pressure coefficient plotting function
def CPplot(panels):
    theta = numpy.linspace(0,2*numpy.pi,len(panels))
    pyplot.figure(figsize=(8,6))
    pyplot.ylabel("$C_P$",fontsize=16)
    pyplot.xlabel(r'$\theta$',fontsize=16)
    pyplot.plot(theta[0:len(theta)/2],pressure(panels)[0:len(panels)/2],lw=2, c='g', label=r'$C_P$')
    pyplot.legend(loc='lower right')



#--------------------------
#|    This is 1.1         |
#--------------------------

# define 3 different resolutions of circular cylinders and plot in the same figure.
#circle32 = new_make_circle(32)
#solve_gamma_kutta(circle32)
#
#circle64 = new_make_circle(64)
#solve_gamma_kutta(circle64)
#
#circle128 = new_make_circle(128)
#solve_gamma_kutta(circle128)
#
#pyplot.figure(figsize=(8,6))
#pyplot.ylabel(r"$C_P$",fontsize=16)
#pyplot.xlabel(r'$\theta (rad)$',fontsize=16)
## Notice the theta value and number here should changed following CPplot function 
## because the range of theta in the figure is from 0 to pi
#pyplot.plot(numpy.linspace(0,2*numpy.pi,len(circle32))[0:len(circle32)/2],pressure(circle32)[0:len(circle32)/2],lw=2, c='c', label='$C_P--32$')
#pyplot.plot(numpy.linspace(0,2*numpy.pi,len(circle64))[0:len(circle64)/2],pressure(circle64)[0:len(circle64)/2],lw=2, c='b', label='$C_P--64$')
#pyplot.plot(numpy.linspace(0,2*numpy.pi,len(circle128))[0:len(circle128)/2],pressure(circle128)[0:len(circle128)/2],lw=2, c='r', label='$C_P--128$')
#pyplot.plot(numpy.linspace(0,numpy.pi,180),1-4*numpy.sin(numpy.linspace(0,numpy.pi,180))**2,lw=2, c='k', label='$C_P--exact$')
#pyplot.legend(loc='lower right')



#--------------------------
#|    This is 1.2         |
#--------------------------

#from BoundaryLayer import ddx_delta, heun, g_pohl, g_1, df_0, march
#lam0 = 7.05232310118

from BoundaryLayer import ddx_delta, heun, g_pohl, g_1, df_0, march
lam0 = 7.05232310118



#

# solve the boundary thickness delta




#
nu = 1e-5                                   # viscosity, based on Reynolds number, nu = UR/Re
N = 128                                     # number of steps
s = numpy.linspace(0,numpy.pi,N)            # distance goes from 0..pi
u_e = 2.*numpy.sin(s)                       # velocity
du_e = 2.*numpy.cos(s)                      # gradient
delta, lam, iSep = march(s,u_e,du_e,nu)     # solve!
#
##tau = nu*1000 *u_e/delta*df_0(lam)
#
## compute frictional drag. int tau ds, where ds is the lenght of vortex panel
#
##print 'drag coefficient = ', drag(tau,circle)/(1./2.*1000*1*numpy.pi)
#
#
#
#
##def c_f(lam, nu, delta, u_e):
##    Re_d = delta*u_e/nu
##    return 2*df_0(lam)/Re_d
#    
def half_c_f(lam, nu, delta, u_e):
    Re_d = delta*u_e/nu
    return df_0(lam)/Re_d
    

def tau_w(lam, nu, delta, u_e):
    if u_e == 0: return 0
    return half_c_f(lam, nu, delta, u_e)*u_e**2
    
#print tau_w(lam[0], nu, delta[0], u_e[0])
#tau = numpy.full_like(delta, 0)
#for i in range(iSep+1):
#    tau[i] = tau_w(lam[i], nu, delta[i], u_e[i])

tau = [tau_w(lam[i], nu, delta[i], u_e[i]) for i in range(N)]


sx = numpy.sin(s[0:iSep+1])

def drag(tau, sx, N):
#    return numpy.sum(tau[:iSep+1]*sx*numpy.pi/(N-1))
    return numpy.trapz(tau[:iSep+1]*sx, dx = numpy.pi/(N-1))


#print(sx[0:iSep])
# Notice: We just only compute half body, so in the end, we need to times 2 to get the whole cylinder's friction



#
#
#C_F_circle = 2*drag(tau, sx, N)/numpy.pi
#print ('Circle frictional coefficient = ' + '%0.2e' %C_F_circle)
#
#C_F_flat = 1.33 * numpy.sqrt(nu/numpy.pi)
#print("Flate plate: "+'%0.2e' %C_F_flat)
#
#s_x = numpy.sin(s)
#pyplot.plot(s, tau*s_x)
#pyplot.scatter(s[iSep], tau[iSep]*s_x[iSep])
#pyplot.xlabel('$s$',size=20)
#pyplot.ylabel(r'$\tau_w s_x$', size=20)
#















#######
# solve the boundary thickness delta
#nu = 1e-5                                   # viscosity, based on Reynolds number, nu = UR/Re
#N = 64                                     # number of steps
#s = numpy.linspace(0,numpy.pi,N)            # distance goes from 0..pi
#u_e = 2.*numpy.sin(s)                       # velocity
#du_e = 2.*numpy.cos(s)                      # gradient
#delta, lam, iSep = march(s,u_e,du_e,nu)     # solve!
#tau = numpy.full_like(delta, 0)
#for i in range(iSep+1):
#    tau[i] = tau_w(lam[i], nu, delta[i], u_e[i])
#
#sx = numpy.sin(s[0:iSep])
#print 'drag coefficient = ', 2*drag(tau, sx, N)/numpy.pi
#
#C_F_flat = 1.33*numpy.sqrt(nu/s[-1])
#print("Flate plate: "+ '%0.2e' %C_F_flat) 
########
#
#
##tau = nu*1000 *u_e/delta*df_0(lam)
#
## compute frictional drag. int tau ds, where ds is the lenght of vortex panel
#
##print 'drag coefficient = ', drag(tau,circle)/(1./2.*1000*1*numpy.pi)
#
#
#
#
##def c_f(lam, nu, delta, u_e):
##    Re_d = delta*u_e/nu
##    return 2*df_0(lam)/Re_d
#    
#def half_c_f(lam, nu, delta, u_e):
#    Re_d = delta*u_e/nu
#    return df_0(lam)/Re_d
#    
#
#def tau_w(lam, nu, delta, u_e):
#    if u_e == 0: return 0
#    return half_c_f(lam, nu, delta, u_e)*u_e**2
#    
##print tau_w(lam[0], nu, delta[0], u_e[0])
#
#
#
#def drag(tau_w, sx, N):
##    return numpy.sum(tau_w[:iSep]*sx*numpy.pi/N)
#    return numpy.trapz(tau_w[:iSep]*sx, dx = numpy.pi/N)


#print(sx[0:iSep])
# Notice: We just only compute half body, so in the end, we need to times 2 to get the whole cylinder's friction




#--------------------------
#|    This is 2.1         |
#--------------------------

#from LiftBody import make_jukowski, lift, jukowski_CL
#
#
#N1, N2, N3 = 32, 64, 128    # define different resolution
#foil1, foil2, foil3 = make_jukowski(N1，0.15), make_jukowski(N2，0.15), make_jukowski(N3，0.15)   # initialize foils with three resolutions
#n = 10                                # this number is the division of AoA
#alpha = numpy.linspace(0,10,n)         # set angle of attack which is a set of number from 0 to 10**0
## define a loop to calculate lift coefficient
#CL1, CL2, CL3 = numpy.empty(n), numpy.empty(n), numpy.empty(n) # initialize lift coefficient container
#for i in range(n):
#    solve_gamma_kutta(foil1,alpha[i]/180*numpy.pi)     # solve for gamma
#    CL1[i] = lift(foil1)
#
#    solve_gamma_kutta(foil2,alpha[i]/180*numpy.pi)     # solve for gamma
#    CL2[i] = lift(foil2)
#
#    solve_gamma_kutta(foil3,alpha[i]/180*numpy.pi)     # solve for gamma
#    CL3[i] = lift(foil3)    
#
#
#pyplot.figure(figsize=(8,6))
#pyplot.ylabel(r"$C_L$",fontsize=16)
#pyplot.xlabel(r'AoA $\alpha (degree)$',fontsize=16)
#pyplot.plot(alpha,CL1,lw=3, c='c', label='$C_L--32$')
#pyplot.plot(alpha,CL2,lw=3, c='m', label='$C_L--64$')
#pyplot.plot(alpha,CL3,lw=3, c='y', label='$C_L--128$')
#pyplot.plot(alpha,jukowski_CL(alpha/180*numpy.pi,0.15),lw=3, c='k', label='$C_L--exact$')
#pyplot.legend(loc='lower right')



#--------------------------
#|    This is 2.2         |
#--------------------------
from Separation import split_panels, Pohlhausen

def new_solve_plot_boundary_layers(panels,alpha=0,nu=1e-5):
    # split the panels
    top_panels,bottom_panels = split_panels(panels)
    # Set up and solve the top boundary layer
    top = Pohlhausen(top_panels,nu)
    top.march()
    # Set up and solve the bottom boundary layer
    bottom = Pohlhausen(bottom_panels,nu)
    bottom.march()
    return top,bottom

def new_predict_jukowski_separation(t_c,alpha=0,N=128):
    # set dx to gets the correct t/c
    foil = make_jukowski(N,dx=t_c-0.019)
    # find and print t/c
    x0 = foil[N/2].xc
    c = foil[0].xc-x0
    t = 2.*numpy.max([p.yc for p in foil])
    # solve potential flow and boundary layer evolution
    solve_gamma_kutta(foil,alpha)
    top, bottom = new_solve_plot_boundary_layers(foil,alpha)
    return (top.x_sep-x0)/c, (bottom.x_sep-x0)/c
#
#
#n = 10                                # this number is the division of AoA
#alpha = numpy.linspace(0,10,n) 
#sep = numpy.empty([n,2])                # Initialize a sep matrix to contain the separation points information.
#for i in range(n):
#    sep[i,:] = new_predict_jukowski_separation(0.15,alpha[i]*numpy.pi/180)
#pyplot.figure(figsize=(8,6))
#pyplot.ylabel(r'Separation position x/c from the leading edge',fontsize=16)
#pyplot.xlabel(r'AoA $\alpha (degree)$',fontsize=16)
#pyplot.plot(alpha,sep[:,0],lw=3, c='r', label='Top Separation')
#pyplot.scatter(alpha, sep[:,0])
#pyplot.plot(alpha,sep[:,1],lw=3, c='b', label='Bottom Separation')
#pyplot.scatter(alpha, sep[:,1])
#pyplot.legend(loc='upper left')
#
#for i in range(n):
#    print 'AoA: '+ '%.3f' %(alpha[i]*numpy.pi/180) + ' degree, Top Separation at x/c= ' + '%.3f'%(sep[i,0])+ ', Bottom Separation at x/c= ' + '%.3f'%(sep[i,1])



# This is going to estimate circular cylinder's separation.
# from LiftBody import make_circle


def new_predict_cylinder_separation(N=128,alpha=0,nu=1.e-5):
    # This is cylinder's separation information.
    s = numpy.linspace(0,2*numpy.pi,N)            # distance goes from 0..pi
    u_e = 2.*numpy.sin(s)                       # velocity
    du_e = 2.*numpy.cos(s) 
    delta, lam, iSep = march(s,u_e,du_e,nu)
    p = numpy.empty([N])                      
    circle = new_make_circle(N)
    solve_gamma_kutta(circle)
    p = pressure(circle)                    # This is pressure value 
    for i in range(N):
        if i > iSep:
            p[i] = p[i-1]
    return p[:N/2]
    
N = 32
p = new_predict_cylinder_separation(N,alpha=0,nu=1.e-5)

pyplot.figure(figsize=(8,6))
pyplot.ylabel(r'Pressure coefficient',fontsize=16)
pyplot.xlabel(r'AoA $\alpha (rad)$',fontsize=16)
pyplot.plot(numpy.linspace(0,numpy.pi,N/2),p,lw=3, c='r', label='$C_P$')
pyplot.legend(loc='upper right')


# Calculate drag coefficient with different resolutions.
N1, N2, N3 = 32, 64, 128
p1, p2, p3 = numpy.empty(N1), numpy.empty(N2), numpy.empty(N3)
p1 = new_predict_cylinder_separation(N1,alpha=0,nu=1.e-5)
p2 = new_predict_cylinder_separation(N2,alpha=0,nu=1.e-5)
p3 = new_predict_cylinder_separation(N3,alpha=0,nu=1.e-5)

def Cd(N,p):
    s = numpy.linspace(0,2*numpy.pi,N)
    sy = numpy.cos(s)
    ds = 2*numpy.pi/(N-1)               # This is interval
    return numpy.trapz(p*sy[:N/2], dx = 2*numpy.pi/(N-1))*2/numpy.pi
    

print 'Resolution N = ' +'%i' %(N1) + ', Drag coefficient=  ' + '%.3f' %(Cd(N1,p1))
print 'Resolution N = ' +'%i' %(N2) + ', Drag coefficient=  ' + '%.3f' %(Cd(N2,p2))
print 'Resolution N = ' +'%i' %(N3) + ', Drag coefficient=  ' + '%.3f' %(Cd(N3,p3))
    

    
    
    




#--------------------------
#|    This is 3           |
#--------------------------

#def naca_offset(x): 
#    return 5*0.15*(0.2969*numpy.sqrt(x)-0.1260*x-0.3516*x**2+0.2843*x**3-0.1036*x**4)
#
#def make_naca(N):
#    theta = numpy.linspace(0,-2.*numpy.pi,N+1)
#    x = (numpy.cos(theta)+1)/2.
#    y = naca_offset(x)*numpy.sign(numpy.sin(theta))
#    naca = numpy.empty(N, dtype=object)
#    for i in range(N):
#        naca[i] = Panel(x[i], y[i], x[i+1], y[i+1])
#    return naca
#    
#
#def predict_naca_separation(N,alpha=0):
#    # set dx to gets the correct t/c
#    naca = make_naca(N)
#
#    x0 = naca[N/2].xc
#    c = naca[0].xc-x0
#    # solve potential flow and boundary layer evolution
#    solve_gamma_kutta(naca,alpha)
#    top, bottom = new_solve_plot_boundary_layers(naca,alpha)
#    return (top.x_sep-x0)/c, (bottom.x_sep-x0)/c
#
#N = 128
#n = 10
#alpha = numpy.linspace(0,10,n) 
#sep = numpy.empty([n,2])                # Initialize a sep matrix to contain the separation points information.
#for i in range(n):
#    sep[i,:] = predict_naca_separation(N,alpha[i]*numpy.pi/180)
#
#pyplot.figure(figsize=(8,6))
#pyplot.ylabel(r'Separation position x/c from the leading edge',fontsize=16)
#pyplot.xlabel(r'AoA $\alpha (degree)$',fontsize=16)
#pyplot.plot(alpha,sep[:,0],lw=3, c='r', label='Top Separation')
#pyplot.scatter(alpha, sep[:,0])
#pyplot.plot(alpha,sep[:,1],lw=3, c='b', label='Bottom Separation')
#pyplot.scatter(alpha, sep[:,1])
#pyplot.legend(loc='upper left') 
#
#naca = make_naca(N)
#alpha = numpy.pi/16
#solve_gamma_kutta(naca,alpha)
#plot_flow(naca,alpha)
#
##Compare with jukouski foil
#n = 10                                # this number is the division of AoA
#alpha = numpy.linspace(0,10,n) 
#sep2 = numpy.empty([n,2])                # Initialize a sep matrix to contain the separation points information.
#for i in range(n):
#    sep2[i,:] = new_predict_jukowski_separation(0.15,alpha[i]*numpy.pi/180)
#
#
#
#pyplot.figure(figsize=(8,6))
#pyplot.ylabel(r'Separation position x/c from the leading edge',fontsize=16)
#pyplot.xlabel(r'AoA $\alpha (degree)$',fontsize=16)
#pyplot.plot(alpha,sep[:,0],lw=3, c='r', label='NACA Top Separation')
#pyplot.plot(alpha,sep2[:,0],lw=3, c='g', label='Jukowski Top Separation')
#pyplot.plot(alpha,sep[:,1],lw=3, c='b', label='NACA Bottom Separation')
#pyplot.plot(alpha,sep2[:,1],lw=3, c='y', label='Jukowski Bottom Separation')
#pyplot.legend(loc='bottom left') 