# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:28:45 2015

@author: yc11e14
"""

import numpy
from matplotlib import pyplot

from VortexPanel import Panel,plot_flow
from LiftBody import solve_gamma_kutta, make_jukowski,make_circle

# Define jukowski foil

alpha = numpy.pi/16
N = 64
foil = make_jukowski(N)
solve_gamma_kutta(foil,alpha)
plot_flow(foil,alpha)

# split panels into two sections based on the flow velocity
def split_panels(panels):
    # positive velocity defines `top` BL
    top = [p for p in panels if p.gamma<=0]      
    # negative defines the `bottom`
    bottom = [p for p in panels if p.gamma>=0]
    # reverse array so panel[0] is stagnation
    bottom = bottom[::-1]

    return top,bottom


def plot_segment(panels):
    pyplot.figure(figsize=(10,2))
    pyplot.axis([-1.2,1.2,-.3,.3])
    for i,p_i in enumerate(panels): 
        p_i.plot()
        if i%10 == 0:
            pyplot.scatter(p_i.xc,p_i.yc)
            pyplot.text(p_i.xc,p_i.yc+0.05, 
                'panel ['+'%i'%i+']',fontsize=12)
       



foil_top,foil_bottom = split_panels(foil)         
plot_segment(foil_top)






# Pohlhausen Boundary Layer class
class Pohlhausen:
    def __init__(self,panels,nu):
        self.u_e = [abs(p.gamma) for p in panels]   # tangential velocity
        self.s = numpy.empty_like(self.u_e)         # initialize distance array
        self.s[0] = panels[0].S
        for i in range(len(self.s)-1):              # fill distance array
            self.s[i+1] = self.s[i]+panels[i].S+panels[i+1].S           
        ds = numpy.gradient(self.s)     
        self.du_e = numpy.gradient(self.u_e,ds)     # compute velocity gradient

        self.nu = nu                                # kinematic viscosity
        self.xc = [p.xc for p in panels]            # x and ...
        self.yc = [p.yc for p in panels]            # y locations
        
    def march(self):
        # march down the boundary layer until separation
        from BoundaryLayer import march
        self.delta,self.lam,self.iSep = march(self.s,self.u_e,self.du_e,self.nu)
    
        # interpolate values at the separation point
        def sep_interp(y): return numpy.interp(    # interpolate function
            12,-self.lam[self.iSep:self.iSep+2],y[self.iSep:self.iSep+2])
        self.s_sep = sep_interp(self.s)
        self.u_e_sep = sep_interp(self.u_e)
        self.x_sep = sep_interp(self.xc)
        self.y_sep = sep_interp(self.yc)
        self.delta_sep = sep_interp(self.delta)

circle = make_circle(N)             # set-up circle
solve_gamma_kutta(circle)           # solve flow
top,bottom = split_panels(circle)   # split panels
nu = 1e-5                           # set viscosity
top = Pohlhausen(top,nu)            # get BL inputs
u_e = 2.*numpy.sin(top.s)           # analytic u_e
du_e = 2.*numpy.cos(top.s)          # analytic du_e

# compare the boundary layer inputs
pyplot.xlabel(r"$s$",fontsize=16)
pyplot.plot(top.s,top.u_e, lw=2, label=r'Panel $u_e$')
pyplot.plot(top.s,u_e, lw=2, label=r'Analytic $u_e$')
pyplot.plot(top.s,top.du_e, lw=2, label=r"Panel $u_e'$")
pyplot.plot(top.s,du_e, lw=2, label=r"Analytic $u_e'$")
pyplot.legend(loc='lower left')


top.march()         # solve the boundary layer flow
i = top.iSep+2      # last point to plot

# plot the boundary layer thicknes and separation point
pyplot.ylabel(r'$\delta$', fontsize=16)
pyplot.xlabel(r'$s$', fontsize=16)
pyplot.plot(top.s[:i],top.delta[:i],lw=2)
pyplot.scatter(top.s_sep,top.delta_sep, s=100, c='r')
pyplot.text(top.s_sep-0.6,top.delta_sep, 
    ' separation \n s='+'%.2f' % top.s_sep,fontsize=12)
    
    
    



def solve_plot_boundary_layers(panels,alpha=0,nu=1e-5):

    # split the panels
    top_panels,bottom_panels = split_panels(panels)
    
    # Set up and solve the top boundary layer
    top = Pohlhausen(top_panels,nu)
    top.march()

    # Set up and solve the bottom boundary layer
    bottom = Pohlhausen(bottom_panels,nu)
    bottom.march()
    
    # plot flow with separation points
    plot_flow(panels,alpha)
    pyplot.scatter(top.x_sep, top.y_sep, s=100, c='r')
    pyplot.scatter(bottom.x_sep, bottom.y_sep, s=100, c='g')
    
    return top,bottom






top,bottom = solve_plot_boundary_layers(circle)
