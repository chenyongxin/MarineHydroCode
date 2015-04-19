# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 23:28:00 2015

@author: ChenY
"""

import numpy
import matplotlib.pyplot
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



# define 3 different resolutions of circular cylinders and plot in the same figure.
circle32 = new_make_circle(32)
solve_gamma_kutta(circle32)

circle64 = new_make_circle(64)
solve_gamma_kutta(circle64)

circle128 = new_make_circle(128)
solve_gamma_kutta(circle128)

pyplot.figure(figsize=(8,6))
pyplot.ylabel(r"$C_P$",fontsize=16)
pyplot.xlabel(r'$\theta$',fontsize=16)
# Notice the theta value and number here should changed following CPplot function 
# because the range of theta in the figure is from 0 to pi
pyplot.plot(numpy.linspace(0,2*numpy.pi,len(circle32))[0:len(circle32)/2],pressure(circle32)[0:len(circle32)/2],lw=2, c='c', label='$C_P--32$')
pyplot.plot(numpy.linspace(0,2*numpy.pi,len(circle64))[0:len(circle64)/2],pressure(circle64)[0:len(circle64)/2],lw=2, c='b', label='$C_P--64$')
pyplot.plot(numpy.linspace(0,2*numpy.pi,len(circle128))[0:len(circle128)/2],pressure(circle128)[0:len(circle128)/2],lw=2, c='r', label='$C_P--128$')
pyplot.plot(numpy.linspace(0,numpy.pi,180),1-4*numpy.sin(numpy.linspace(0,numpy.pi,180))**2,lw=2, c='k', label='$C_P--exact$')
pyplot.legend(loc='lower right')


