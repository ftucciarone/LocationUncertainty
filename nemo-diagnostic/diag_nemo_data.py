import numpy as np
from scipy import signal
from netCDF4 import Dataset


class DIAG_NEMO:
    """
    Diagnostics of NEMO data.
    """
    
    def __init__(self, param):

        # Read param.
        for kw in param:
            setattr(self, kw, param[kw]) 

        # Collect data
        self.read_data()
         
        # Create output 
        self.create_outfile()
   
        # Derived param.
        self.fs = 1/self.dt # sampling freq. [s^-1]
        self.mask_avg = np.zeros_like(self.sst[0]) # mask used for spatial average of PSDs
        self.mask_avg[self.nb:-self.nb,self.nb:-self.nb] = 1./np.prod(self.mask_avg[self.nb:-self.nb,self.nb:-self.nb].shape)

    def read_data(self):
        print(f'Reading data from {self.infiles_grid_U[0]}')
        nc = Dataset(self.infiles_grid_U[0],'r')
        # Read dimensions
        self.nt = nc.dimensions['time_counter'].size
        self.nz = nc.dimensions['depthu'].size
        self.ny = nc.dimensions['y'].size
        self.nx = nc.dimensions['x'].size
        # Read variables
        self.u = nc.variables['vozocrtx'][:,0].data.astype(self.dtype) # only surface  
        nc.close()
        
        print(f'Reading data from {self.infiles_grid_V[0]}')
        nc = Dataset(self.infiles_grid_V[0],'r')
        self.v = nc.variables['vomecrty'][:,0].data.astype(self.dtype) 
        nc.close()
        
        print(f'Reading data from {self.infiles_grid_W[0]}')
        nc = Dataset(self.infiles_grid_W[0],'r')
        self.w = nc.variables['vovecrtz'][:,0].data.astype(self.dtype) 
        nc.close()

        print(f'Reading data from {self.infiles_grid_T[0]}')
        nc = Dataset(self.infiles_grid_T[0],'r')
        self.ssh = nc.variables['sossheig'][...].data.astype(self.dtype) 
        self.sst = nc.variables['votemper'][:,0].data.astype(self.dtype)  
        nc.close()

        # Collect sequence of data
        for i in range(1,len(self.infiles_grid_U)):
            # U files
            nc = Dataset(self.infiles_grid_U[i],'r')
            self.nt += nc.dimensions['time_counter'].size
            self.u = np.append(self.u, nc.variables['vozocrtx'][:,0].data.astype(self.dtype), axis=0) 
            nc.close()
            # V files
            nc = Dataset(self.infiles_grid_V[i],'r')
            self.v = np.append(self.v, nc.variables['vomecrty'][:,0].data.astype(self.dtype), axis=0) 
            nc.close()
            # W files
            nc = Dataset(self.infiles_grid_W[i],'r')
            self.w = np.append(self.w, nc.variables['vovecrtz'][:,0].data.astype(self.dtype), axis=0) 
            nc.close()
            # T files
            nc = Dataset(self.infiles_grid_T[i],'r')
            self.ssh = np.append(self.ssh, nc.variables['sossheig'][...].data.astype(self.dtype), axis=0) 
            self.sst = np.append(self.sst, nc.variables['votemper'][:,0].data.astype(self.dtype), axis=0) 
            nc.close()
        

    def create_outfile(self):
        print(f'Initializing file: {self.outfile}')
        self.outnc = Dataset(self.outfile,'w',format='NETCDF4')
        # Create dimensions
        self.outnc.createDimension('t', self.nt)
        self.outnc.createDimension('y', self.ny)
        self.outnc.createDimension('x', self.nx)
        self.outnc.createDimension('f', self.nt//2+1)
        # Create variables
        dtype = 'f8' if self.dtype=='float64' else 'f4'
        self.outnc.createVariable('f', dtype, ('f',))
        self.outnc.createVariable('sst_mean', dtype, ('y','x',))
        self.outnc.createVariable('sst_std', dtype, ('y','x',))
        self.outnc.createVariable('sst_tpsd', dtype, ('f',))
        self.outnc.createVariable('ssh_mean', dtype, ('y','x',))
        self.outnc.createVariable('ssh_std', dtype, ('y','x',))
        self.outnc.createVariable('ssh_tpsd', dtype, ('f',))
        # Add attributes
        self.outnc.variables['f'].long_name = 'Frequency'
        self.outnc.variables['f'].units = 's^-1'
        self.outnc.variables['sst_mean'].long_name = 'Temporal mean of SST'
        self.outnc.variables['sst_mean'].units = 'K'
        self.outnc.variables['sst_std'].long_name = 'Temporal standard deviation of SST'
        self.outnc.variables['sst_std'].units = 'K'
        self.outnc.variables['sst_tpsd'].long_name = 'Temporal power spectral density of SST'
        self.outnc.variables['ssh_mean'].long_name = 'Temporal mean of SSH'
        self.outnc.variables['ssh_mean'].units = 'm'
        self.outnc.variables['ssh_std'].long_name = 'Temporal standard deviation of SSH'
        self.outnc.variables['ssh_std'].units = 'm'
        self.outnc.variables['ssh_tpsd'].long_name = 'Temporal power spectral density of SSH'


    def compute_temporal_PSD(self, q):
        if hasattr(self,'freq'):
            _, psd = signal.periodogram(q, self.fs, axis=0)
        else:  
            self.freq, psd = signal.periodogram(q, self.fs, axis=0)
        return np.sum(psd*self.mask_avg[None], axis=(-2,-1))


    def run(self):
        print('Run diagnostics')
        self.outnc.variables['sst_mean'][...] = np.mean(self.sst, axis=0)
        self.outnc.variables['sst_std'][...] = np.std(self.sst, axis=0)
        self.outnc.variables['sst_tpsd'][:] = self.compute_temporal_PSD(self.sst)
        self.outnc.variables['ssh_mean'][...] = np.mean(self.ssh, axis=0)
        self.outnc.variables['ssh_std'][...] = np.std(self.ssh, axis=0)
        self.outnc.variables['ssh_tpsd'][:] = self.compute_temporal_PSD(self.ssh)
        self.outnc.variables['f'][:] = self.freq
        self.outnc.close()


if __name__ == "__main__":
    
    dir_nemo = '/srv/storage/ithaca@storage2.rennes.grid5000.fr/ftucciarone/nemo-results/'
    
    # Path of R27 det. data
    dir_mode = 'data_R27'
    dir_subs = ['100-102y','102-104y','104-106y','106-108y','108-110y']
    filename = 'GYRE_5d_00010101_00021230'
    import os
    infiles_grid_U, infiles_grid_V, infiles_grid_W, infiles_grid_T = [], [], [], []
    for s in dir_subs:
        dir1 = os.path.join(dir_nemo, dir_mode, s)
        infiles_grid_U = np.append(infiles_grid_U, os.path.join(dir1, f'{filename}_grid_U.nc'))
        infiles_grid_V = np.append(infiles_grid_V, os.path.join(dir1, f'{filename}_grid_V.nc'))
        infiles_grid_W = np.append(infiles_grid_W, os.path.join(dir1, f'{filename}_grid_W.nc'))
        infiles_grid_T = np.append(infiles_grid_T, os.path.join(dir1, f'{filename}_grid_T.nc'))
    outfile = os.path.join(dir_nemo, dir_mode, 'diag.nc')

    # Set param
    param = {
        'dt': 5.*86400, # deta step [s]
        'nb': 4*3+1, # pts. near boundary (including boundary)
        'infiles_grid_U': infiles_grid_U,
        'infiles_grid_V': infiles_grid_V,
        'infiles_grid_W': infiles_grid_W,
        'infiles_grid_T': infiles_grid_T,
        'outfile': outfile,
        'dtype': 'float32', # 'float32' or 'float64'
    }   
    
    # Run diag
    diag = DIAG_NEMO(param)
    diag.run()

    # Check results
    nc = Dataset(outfile,'r')
    sst_mean = nc.variables['sst_mean'][...].data
    sst_std = nc.variables['sst_std'][...].data
    freq = nc.variables['f'][:].data
    sst_tpsd = nc.variables['sst_tpsd'][:].data
    nc.close()
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,2, constrained_layout=True, figsize=(8,4))
    ax[0].imshow(sst_mean, origin='lower', cmap='RdBu_r')
    ax[0].set_title('Temporal mean of SST')
    ax[1].imshow(sst_std, origin='lower', cmap='Reds')
    ax[1].set_title('Temporal std of SST')
    plt.show()
    
    fig, ax = plt.subplots(constrained_layout=True, figsize=(6,4))
    ax.loglog(freq, sst_tpsd)
    ax.set(xlabel='freq [1/s]', ylabel='PSD', title='Temporal spectra of SST')
    plt.grid(which='both', axis='both')
    plt.show()
