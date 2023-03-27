"""" Module to perform the Snapshot Proper Orthogonal Decomposition
     POD) procedure """

# Load modules
# ------------
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# Subroutines
# -----------
def RootMeanSquare(observable, weights, r1_weight_sum = None):

    if len(np.shape(weights)) == 2:
       # If weights are not present, compute them;
       if r1_weight_sum == None:
           r1_weight_sum = 1. / np.sum(weights, axis=(-2,-1))
       # If the shape of the observable is the same as of the 
       # weights, then it's a matter of scalar products
       if len(np.shape(observable)) == len(np.shape(weights)):
           RMS = np.sum(np.multiply(observable**2, weights), axis=(-2, -1))
           RMS = RMS * r1_weight_sum
       # Else, the product will be computed slice by slice and
       # the result is going to be an array
       else:
           RMS = np.sum(np.multiply(observable**2, weights[np.newaxis,:,:]), axis=(-2, -1))
           RMS = RMS * r1_weight_sum

    else:
       # If weights are not present, compute them;
       if r1_weight_sum == None:
           r1_weight_sum = 1. /  np.sum(weights, axis=(-3,-2,-1))
       RMS = np.sum(np.multiply(observable**2, weights), axis=(-3, -2, -1))
       RMS = RMS * r1_weight_sum

    RMS = np.sqrt(RMS)
    return RMS



def RootMeanSquareERROR(estimator, observable, weights, r1_weight_sum = None):
    error = estimator - observable
    if len(np.shape(weights)) == 2:
       # If weights are not present, compute them;
       if r1_weight_sum == None:
           r1_weight_sum = 1. / np.sum(weights, axis=(-2,-1))
       # If the shape of the observable is the same as of the 
       # weights, then it's a matter of scalar products
       if len(np.shape(error)) == len(np.shape(weights)):
           RMSE = np.sum(np.multiply(error**2, weights), axis=(-2, -1))
           RMSE = RMSE * r1_weight_sum
       # Else, the product will be computed slice by slice and
       # the result is going to be an array
       else:
           RMSE = np.sum(np.multiply(error**2, weights[np.newaxis,:,:]), axis=(-2, -1))
           RMSE = RMSE * r1_weight_sum

    else:
       # If weights are not present, compute them;
       if r1_weight_sum == None:
           r1_weight_sum =  1. / np.sum(weights, axis=(-3,-2,-1))
       RMSE = np.sum(np.multiply(error**2, weights), axis=(-3, -2, -1))
       RMSE = RMSE * r1_weight_sum

    RMSE = np.sqrt(RMSE)
    return RMSE
 


def PatternCorrelation(estimator, observable, weights):

    arg = np.multiply(estimator, observable)
    num = np.sum(np.multiply(arg, weights), axis=(-3, -2, -1))
    int1 = np.sum( np.multiply(estimator**2, weights), axis=(-3, -2, -1))
    int2 = np.sum( np.multiply(observable**2, weights), axis=(-3, -2, -1))
    den = int1 * int2
      
    PC = num / den
    
    return PC



def GaussianRelativeEntropy(est_avg, obs_avg, est_vnc, obs_vnc):

    tmp1 = np.divide((est_avg - obs_avg)**2, est_vnc)
    tmp2 = np.divide(obs_vnc, est_vnc)
    tmp3 = 1 + np.log(tmp2)
    GRE = 0.5 * (tmp1 + tmp2 - tmp3)
    
    return GRE



def AveragedGaussianRelativeEntropy(est_avg, obs_avg, est_vnc, obs_vnc, weights):

    GRE = GaussianRelativeEntropy(est_avg, obs_avg, est_vnc, obs_vnc)
    r1_weight_sum = 1. / np.sum(weights, axis=(-3,-2,-1))
    intGRE = np.sum(np.multiply(GRE, weights), axis=(-3, -2, -1))  
    wGRE = intGRE * r1_weight_sum
    
    return wGRE
    
