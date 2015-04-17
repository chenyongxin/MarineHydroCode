# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 23:28:00 2015

@author: ChenY
"""

import numpy
from VortexPanel import Panel,plot_flow
from LiftBody import solve_gamma_kutta, make_jukowski,make_circle
# 1. A circular cylinder
# 1.1 N = 32

circle = make_circle(128)
solve_gamma_kutta(circle)

#first of all, define the coordinate of circle, then use these coordinates to find the velocity.
#directly use flow_velocity to find the velocity along the surface

def pressure(panels):
    x = numpy.array([p.xc for p in panels])
    y = numpy.array([p.yc for p in panels])
    u,v = flow_velocity(panels,x,y)
    return 2*numpy.sum( numpy.multiply([1-x**2-y**2],[p.S for p in panels])) 
print 'C_P = ', pressure(circle)
