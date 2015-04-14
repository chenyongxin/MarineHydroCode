# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:02:02 2015

@author: yc11e14
"""

import numpy
from matplotlib import pyplot
from VortexPanel import Panel, solve_gamma, plot_flow




def make_circle(N):
    # define the end-points of the panels
    x_ends = numpy.cos(numpy.linspace(0, -2*numpy.pi, N+1))
    y_ends = numpy.sin(numpy.linspace(0, -2*numpy.pi, N+1))

    # define the panels
    circle = numpy.empty(N, dtype=object)
    for i in xrange(N):
        circle[i] = Panel(x_ends[i], y_ends[i], x_ends[i+1], y_ends[i+1])

    return circle

#N = 16
#circle = make_circle(N)  # set-up geom
#solve_gamma(circle)      # find gamma
#plot_flow(circle)        # plot the flow

def make_jukowski( N, dx = 0.2, dy = 0, dr = 0 ):
    # define the circle
    x_ends = numpy.cos(numpy.linspace(0, -2*numpy.pi, N+1))
    y_ends = numpy.sin(numpy.linspace(0, -2*numpy.pi, N+1))

    # shift circle
    r = numpy.sqrt((1+dx)**2+dy**2)+dr
    x2_ends = r*x_ends-dx
    y2_ends = r*y_ends-dy
    r2_ends = x2_ends**2+y2_ends**2

    # apply jukowski mapping
    x3_ends = x2_ends*(1+1./r2_ends)/2
    y3_ends = y2_ends*(1-1./r2_ends)/2
    
    # define the panels
    foil = numpy.empty(N, dtype=object)
    for i in xrange(N):
        foil[i] = Panel(x3_ends[i], y3_ends[i], x3_ends[i+1], y3_ends[i+1])
    
    return foil

#foil = make_jukowski(N)  # make foil
#solve_gamma(foil)        # solve for gamma
#plot_flow(foil)          # plot the flow
def lift(panels):
    c = panels[0].x[0]-panels[len(panels)/2].x[0]      # length scale
    return -4./c*numpy.sum([p.gamma*p.S for p in panels])
#print 'C_L =',lift(foil)

#alpha = numpy.pi/16         # set angle of attack
#solve_gamma(foil,alpha)     # solve for gamma
#plot_flow(foil,alpha)       # plot the flow
#print 'C_L =',lift(foil)    # print the lift
#
#alpha = 0                               # angle of attack
#dgamma = 0./(2*numpy.pi)                # vortex circulation
#N = 16                                  # number of panels
#
#circle = make_circle(N)                 # set-up geom
#solve_gamma(circle,alpha)               # find gamma
#for p in circle: p.gamma += dgamma      # increase gamma by dgamma
#plot_flow(circle,alpha)                 # plot the flow
#print 'C_L =',lift(circle)              # print the lift

# determine the vortex panel strength with Kutta Condition
def solve_gamma_kutta(panels,alpha=0):
    from VortexPanel import construct_A_b
    A,b = construct_A_b(panels,alpha)   # construct linear system
    A[:, 0] += 1                        # gamma[0]+ ...
    A[:,-1] += 1                        # gamma[N-1]=0
    gamma = numpy.linalg.solve(A, b)    # solve for gamma!
    for i,p_i in enumerate(panels):
        p_i.gamma = gamma[i]            # update panels
        
#N = 16
#alpha = numpy.pi/32
#foil = make_jukowski(N)             # make foil
#solve_gamma_kutta(foil,alpha)       # solve for gamma
#plot_flow(foil,alpha)               # plot the flow
#print 'C_L =',lift(foil)            # print the lift