#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 10:52:21 2023

@author: ftucciar
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt


def ke_compute(u,v):
    ke = np.zeros_like(u)
    up2 = u**2
    vp2 = v**2
    ke[:, :, 1:-1, 1:-1] = ( up2[:, :, 1:-1, 0:-2] + up2[:, :, 1:-1, 1:-1] + \
                             vp2[:, :, 0:-2, 1:-1] + vp2[:, :, 1:-1, 1:-1] )
    ke = ke * 0.25
    return ke