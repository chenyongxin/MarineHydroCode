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

circle = make_circle(32)

