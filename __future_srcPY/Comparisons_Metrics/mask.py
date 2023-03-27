"""" Module to perform the Snapshot Proper Orthogonal Decomposition
     POD) procedure """

# Load modules
# ------------
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# Subroutines
# -----------
def createJetStreamMask(lat_in, lon_in):
    latct = 29
    lonct = -85
    latwd = 5
    lonwd = 20
    
    
    lat_out = lat_in
    lon_out = lon_in
    
    lon_out = np.where(lon_in > lonct - lonwd, lon_in, 361)
    lon_out = np.where(lon_out < lonct + lonwd, lon_out, 361)
    
    lat_out = np.where(lat_in > latct - latwd, lat_in, 91)
    lat_out = np.where(lat_out < latct + latwd, lat_out, 91)
    mask_out = np.where(lat_out != 91, np.ones_like(lat_out), 0)
    mask_out = np.multiply(np.where(lon_out != 361, np.ones_like(lat_out), 0), mask_out)

    border = 5
    mask_out[-border:, :] = 0
    mask_out[:, :border] = 0

    return lat_out, lon_out, mask_out


def findOverlap(HRlat, HRlon, LRlat, LRlon):
    # Get dimensions of input
    HRny = np.shape(HRlat)[0]
    HRnx = np.shape(HRlat)[1]
    LRny = np.shape(LRlat)[0]
    LRnx = np.shape(LRlat)[1]
    # Get ratio between inputs
    xR = int((HRnx - 2) / (LRnx - 2))
    yR = int((HRny - 2) / (LRny - 2))
    # Create replicating pattern in Y and X
    repY = yR*np.ones(LRny-2).astype(int)
    repX = xR*np.ones(LRnx-2).astype(int)
    # Replicate latitude and longitude in Y
    repLon = np.repeat(LRlon[1:-1, 1:-1], repeats=repY, axis=0)
    repLat = np.repeat(LRlat[1:-1, 1:-1], repeats=repY, axis=0)
    # Replicate latitude and longitude in Y
    repLon = np.repeat(repLon, repeats=repX, axis=1)
    repLat = np.repeat(repLat, repeats=repX, axis=1)
    # Find the overlaps traces (these are lines of equilatide and equilongitude)
    tr1 = repLon - HRlon[1:-1, 1:-1]
    tr2 = repLat - HRlat[1:-1, 1:-1]
    # Get the poins of intersections
    mask = tr1 == tr2

    
    cols = np.linspace(0, HRnx-1, HRnx)
    rows = np.linspace(0, HRny-1, HRny)
    cols = cols[1:-1]
    rows = rows[1:-1]

    cols = cols[np.sum(mask, axis=0) < 60]
    rows = rows[np.sum(mask, axis=1) < 1]
    
    return cols, rows