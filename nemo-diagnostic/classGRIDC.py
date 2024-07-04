#!/usr/bin/env python3
#! -*- coding: utf-8 -*-


# Load modules
# ------------
import numpy as np

def ke_compute(u,v):
    ke = np.zeros_like(u)
    up2 = u**2
    vp2 = v**2
    ke[..., 1:-1, 1:-1] = ( up2[..., 1:-1, 0:-2] + up2[..., 1:-1, 1:-1] + \
                            vp2[..., 0:-2, 1:-1] + vp2[..., 1:-1, 1:-1] )
    return ke

