#!/usr/bin/env python

import numpy as np


N = 10000
gaussian_parameters = np.array([ [250, 0.75, 0.15], #intensity, location, width
                                 [170, 1.0, 0.1],
                                 [100,  1.4, 0.15],
                                 [70,  1.7, 0.1] ])

samples = np.empty(0)

for p in gaussian_parameters:
  s = np.random.normal( p[1], p[2], int(N*p[0]/np.sum(gaussian_parameters[:,0])) )
  samples = np.append(s, samples)

np.savetxt('example.dat', samples)
