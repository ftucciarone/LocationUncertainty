import numpy as np
from scipy import signal
from netCDF4 import Dataset
import classArakawaC as Cgrid
import Spectrums as PSpD
import matplotlib.pyplot as plt


class DIAG_NEMO:
    """
    Diagnostics of NEMO data.
    """
    
    def __init__(self, param, klev=None):

        # Read param.
        for kw in param:
            setattr(self, kw, param[kw]) 
            
        # Safely overrides delay if not present 
        if 'dl' not in param.keys():
            self.dl = 0

        # Read dimensions from domain_cfg.nc
        self.read_dims()
        if klev==None:
            self.klev = [item for item in range(0, self.nz)]
        else:
            self.nz = len(self.klev)
        # Compute dimension for spatial spectrum
        self.sk2, self.sk1, efc = PSpD.init_2Dwvnb(self.nx, self.ny, self.dx, self.dx)
        
        self.read_data(self.klev)
        # Create output 
        self.create_outfile()
        # Compute means
        self.compute_mean()
        # Compute energy
        self.compute_energies()
        # Compute Hovmoller plots
        self.compute_hovmoller()
        # Compute spatial Power spectral density
        self.compute_spatial_PSD()
        # Derived param.
        self.fs = 1/self.dt # sampling freq. [s^-1]
        self.mask_avg = np.zeros_like(self.sst[0]) # mask used for spatial average of PSDs
        self.mask_avg[self.nb:-self.nb,self.nb:-self.nb] = 1./np.prod(self.mask_avg[self.nb:-self.nb,self.nb:-self.nb].shape)
        # Compute temporal Power spectral density
        self.u_tpsd = self.compute_temporal_PSD(self.u[self.dl:,...])
        self.v_tpsd = self.compute_temporal_PSD(self.v[self.dl:,...])
        self.tke_tpsd = (self.u_tpsd + self.v_tpsd)/2
        self.save_temporal_PSD('tke_tpsd', self.tke_tpsd)
        # Close
        self.outnc.close()


    def read_dims(self):
        print(f'Reading dimension from {self.infiles_grid}')
        nc = Dataset(self.infiles_grid,'r')
        # Read dimensions
        self.nz = nc.dimensions['nav_lev'].size
        self.ny = nc.dimensions['y'].size
        self.nx = nc.dimensions['x'].size
        self.nt = 0
        #
        self.dx = nc.variables['e1t'][0, 0, 0].data
        self.dy = nc.variables['e2t'][0, 0, 0].data
        # Close domain_cfg
        nc.close()

    def read_data(self, levs):
        # Collect 
        self.u = np.zeros([0, self.nz, self.ny, self.nx])
        self.v = np.zeros([0, self.nz, self.ny, self.nx])
        self.w = np.zeros([0, self.nz, self.ny, self.nx])
        self.t = np.zeros([0, self.nz, self.ny, self.nx])
        self.vort = np.zeros([0, self.nz, self.ny, self.nx])
        self.ssh = np.zeros([0, self.ny, self.nx])
        self.sst = np.zeros([0, self.ny, self.nx])
 
        # Collect data from sequence of U files
        for i in range(len(self.infiles_grid_U)):
            print(f'Reading data from {self.infiles_grid_U[i]}')
            # U files
            nc = Dataset(self.infiles_grid_U[i],'r')
            self.nt += nc.dimensions['time_counter'].size
            self.u = np.append(self.u, nc.variables['vozocrtx'][:,levs,...].data.astype(self.dtype), axis=0) 
            nc.close()
            
        # Collect data from sequence of V files
        for i in range(len(self.infiles_grid_V)):
            print(f'Reading data from {self.infiles_grid_V[i]}')
            # V files
            nc = Dataset(self.infiles_grid_V[i],'r')
            self.v = np.append(self.v, nc.variables['vomecrty'][:,levs,...].data.astype(self.dtype), axis=0) 
            nc.close()

        # Collect data from sequence of W files
        for i in range(len(self.infiles_grid_W)):
            print(f'Reading data from {self.infiles_grid_W[i]}')            
            # W files
            nc = Dataset(self.infiles_grid_W[i],'r')
            self.w = np.append(self.w, nc.variables['vovecrtz'][:,levs,...].data.astype(self.dtype), axis=0) 
            nc.close()

        # Collect data from sequence of T files
        for i in range(len(self.infiles_grid_T)):
            print(f'Reading data from {self.infiles_grid_T[i]}')              
            # T files
            nc = Dataset(self.infiles_grid_T[i],'r')
            self.ssh = np.append(self.ssh, nc.variables['sossheig'][...].data.astype(self.dtype), axis=0) 
            self.sst = np.append(self.sst, nc.variables['votemper'][:,0, ...].data.astype(self.dtype), axis=0) 
            self.t = np.append(self.t, nc.variables['votemper'][:,levs, ...].data.astype(self.dtype), axis=0) 
            
            if 'vorticity' in nc.variables.keys():
                self.vort = np.append(self.vort, nc.variables['vorticity'][:,levs, ...].data.astype(self.dtype), axis=0) 
            nc.close()

 
    def create_outfile(self):
        print(f'Initializing file: {self.outfile}')
        self.outnc = Dataset(self.outfile,'w',format='NETCDF4')
        # Create dimensions
        self.outnc.createDimension('t', self.nt)
        self.outnc.createDimension('z', self.nz)
        self.outnc.createDimension('y', self.ny)
        self.outnc.createDimension('x', self.nx)
        self.outnc.createDimension('f', (self.nt - self.dl)//2+1)
        self.outnc.createDimension('sk1', len(self.sk1)-1)
        
        
        # Create variables
        dtype = 'f8' if self.dtype=='float64' else 'f4'
        
        self.outnc.createVariable('sk1', dtype, ('sk1',))
        self.outnc.variables['sk1'][:] = self.sk1[1:]
        self.outnc.variables['sk1'].units = 'rad m^-1'
        
        self.outnc.createVariable('f', dtype, ('f',))
        # Create velocities
        self.outnc.createVariable('u_mean', dtype, ('z', 'y', 'x',))
        self.outnc.createVariable('v_mean', dtype, ('z', 'y', 'x',))
        self.outnc.createVariable('w_mean', dtype, ('z', 'y', 'x',))
        # Create vorticity 
        self.outnc.createVariable('vort_mean', dtype, ('z', 'y', 'x',))
        # Create 2d fields
        self.outnc.createVariable('sst_mean', dtype, ('y','x',))
        self.outnc.createVariable('ssh_mean', dtype, ('y','x',))
        self.outnc.createVariable('sst_std', dtype, ('y','x',))
        self.outnc.createVariable('ssh_std', dtype, ('y','x',))
        # Create kinetic energies: MEAN KE
        self.outnc.createVariable('hor_mke', dtype, ('z', 'y', 'x',))
        self.outnc.createVariable('ver_mke', dtype, ('z', 'y', 'x',))
        # Create kinetic energies: TURBULENT KE
        self.outnc.createVariable('hor_tke', dtype, ('z', 'y', 'x',))
        self.outnc.createVariable('ver_tke', dtype, ('z', 'y', 'x',))
        # Create kinetic energies: TOTAL KE
        self.outnc.createVariable('hor_ke', dtype, ('z', 'y', 'x',))
        self.outnc.createVariable('ver_ke', dtype, ('z', 'y', 'x',))
        # Create enstrophy: MEAN ENSTROPHY
        self.outnc.createVariable('menst', dtype, ('z', 'y', 'x',))
        # Create Hovmoller fields
        self.outnc.createVariable('t_hov', dtype, ('z', 't',))
        self.outnc.createVariable('u_hov', dtype, ('z', 't',))
        self.outnc.createVariable('v_hov', dtype, ('z', 't',))
        self.outnc.createVariable('w_hov', dtype, ('z', 't',))
        self.outnc.createVariable('tke_hov', dtype, ('z', 't',))
        # Create spatial spectrums
        self.outnc.createVariable('tke_spsd', dtype, ('z', 'sk1',))
        # Create temporal spectrums
        self.outnc.createVariable('tke_tpsd', dtype, ('z', 'f',))
        # Create spectrums
        self.outnc.createVariable('sst_tpsd', dtype, ('f',))
        self.outnc.createVariable('ssh_tpsd', dtype, ('f',))
        #
        #
        # Add attributes
        self.outnc.variables['f'].long_name = 'Frequency'
        self.outnc.variables['f'].units = 's^-1'
        # mean velocity
        self.outnc.variables['u_mean'].long_name = 'Temporal mean of U vel'
        self.outnc.variables['u_mean'].units = 'm/s'
        self.outnc.variables['v_mean'].long_name = 'Temporal mean of V vel'
        self.outnc.variables['v_mean'].units = 'm/s'
        self.outnc.variables['w_mean'].long_name = 'Temporal mean of W vel'
        self.outnc.variables['w_mean'].units = 'm/s'
        # mean vorticity
        self.outnc.variables['vort_mean'].long_name = 'Temporal mean of Vort'
        self.outnc.variables['vort_mean'].units = '1/s'
        # mean sea surface temperature
        self.outnc.variables['sst_mean'].long_name = 'Temporal mean of SST'
        self.outnc.variables['sst_mean'].units = 'K'
        # mean sea surface height 
        self.outnc.variables['ssh_mean'].long_name = 'Temporal mean of SSH'
        self.outnc.variables['ssh_mean'].units = 'm'
        # std of sea surface temperature
        self.outnc.variables['sst_std'].long_name = 'Temporal standard deviation of SST'
        self.outnc.variables['sst_std'].units = 'K'
        # std of sea surface height 
        self.outnc.variables['ssh_std'].long_name = 'Temporal standard deviation of SSH'
        self.outnc.variables['ssh_std'].units = 'm'
        # hovmoller plots
        self.outnc.variables['t_hov'].long_name = 'Hovmoller plot of temperature'
        self.outnc.variables['t_hov'].units = 'K'
        self.outnc.variables['u_hov'].long_name = 'Hovmoller plot of u tke'
        self.outnc.variables['u_hov'].units = 'm^2/s^2'
        self.outnc.variables['v_hov'].long_name = 'Hovmoller plot of v tke'
        self.outnc.variables['v_hov'].units = 'm^2/s^2'
        self.outnc.variables['w_hov'].long_name = 'Hovmoller plot of w tke'
        self.outnc.variables['w_hov'].units = 'm^2/s^2'
        self.outnc.variables['tke_hov'].long_name = 'Hovmoller plot of tke'
        self.outnc.variables['tke_hov'].units = 'm^2/s^2'
        # Spatial spectrums
        self.outnc.variables['tke_spsd'].long_name = 'Spatial TKE PSD'
        self.outnc.variables['tke_spsd'].units = 'rad^-1 m^3 s^-2'
        # Temporal spectrums
        self.outnc.variables['tke_tpsd'].long_name = 'Temporal TKE PSD'
        self.outnc.variables['tke_tpsd'].units = 'rad^-1 m^3 s^-2'


    def compute_temporal_PSD(self, q):
        if hasattr(self,'freq'):
            _, psd = signal.periodogram(q, self.fs, axis=0)
        else:  
            self.freq, psd = signal.periodogram(q, self.fs, axis=0)
        return np.sum(psd*self.mask_avg[None], axis=(-2,-1))

    def save_temporal_PSD(self, string, q):
        self.outnc.variables[string][...] = q[...].T
    
    
    def compute_spatial_PSD(self):
        # Energy in Fourier space is 
        # KE(k) = 0.5 \oint_{|\kappa|=k} ( |\hat{u}|^2 + |\hat{v}|^2 ) d\kappa
        #
        print('Compute and save spatial Power Spectral Density')
        self.fh = 0.5*(np.abs(np.fft.fft2(self.u[self.dl:,...], axes=(-2, -1)))**2 + \
                       np.abs(np.fft.fft2(self.v[self.dl:,...], axes=(-2, -1)))**2 )

        self.sk2, self.sk1, efc = PSpD.init_2Dwvnb(self.nx, self.ny, self.dx, self.dx)
        self.spsd = np.mean(PSpD.iso_spec(self.fh, self.sk2, self.sk1, efc), axis=(0))
        # Spatial spectrums
        self.outnc.variables['tke_spsd'][...] = self.spsd[...,1:]
        

    def compute_mean(self):
        print('Compute and save mean fields')
        # Mean velocities
        self.umean = np.mean(self.u[self.dl:,...], axis=0)
        self.vmean = np.mean(self.v[self.dl:,...], axis=0)
        self.wmean = np.mean(self.w[self.dl:,...], axis=0)
        # Mean kinetic energy
        self.hor_mke = Cgrid.ke_compute(self.umean, self.vmean)
        self.ver_mke = self.wmean**2
        # Mean enstrophy
        self.mvort = np.mean(self.vort[self.dl:,...], axis=0)
        self.menst = self.mvort**2
        # Save
        self.outnc.variables['u_mean'][...] = self.umean
        self.outnc.variables['v_mean'][...] = self.vmean
        self.outnc.variables['w_mean'][...] = self.wmean
        self.outnc.variables['sst_mean'][...] = np.mean(self.sst, axis=0)
        self.outnc.variables['ssh_mean'][...] = np.mean(self.ssh, axis=0)
        self.outnc.variables['vort_mean'][...] = self.mvort
        # Save Energies
        self.outnc.variables['hor_mke'][...] = self.hor_mke
        self.outnc.variables['ver_mke'][...] = self.ver_mke
        # Save Enstrophy
        self.outnc.variables['menst'][...] = self.menst
        
    def compute_std(self):
        print('Compute standard deviation fields')
        self.outnc.variables['sst_std'][...] = np.std(self.sst[self.dl:,..

.], axis=0)
        self.outnc.variables['ssh_std'][...] = np.std(self.ssh[self.dl:,...], axis=0)

    def compute_energies(self):
        print('Compute and save turbulent kinetic energy fields')
        # Turbulent kinetic energies
        self.utke = np.var(self.u[self.dl:,...], axis=0)
        self.vtke = np.var(self.v[self.dl:,...], axis=0)
        self.wtke = np.var(self.w[self.dl:,...], axis=0)
        # Horizontal and vertical Kinetic energies
        self.hor_tke = np.mean(Cgrid.ke_compute(self.utke, self.vtke), axis=0)
        self.ver_tke = np.mean(self.wtke, axis=0)
        # Total kinetic energy
        self.hor_ke = self.hor_mke + self.hor_tke
        self.ver_ke = self.ver_mke + self.ver_tke

        # Save Energies
        self.outnc.variables['hor_tke'][...] = self.hor_tke
        self.outnc.variables['ver_tke'][...] = self.ver_tke
        self.outnc.variables['hor_ke'][...] = self.hor_ke
        self.outnc.variables['ver_ke'][...] = self.ver_ke


    def compute_hovmoller(self):
        print('Compute and save hovmoller plots')
        # Turbulent kinetic energies
        utke = (self.u - self.u.mean(axis=0, keepdims=True))**2
        vtke = (self.v - self.v.mean(axis=0, keepdims=True))**2
        wtke = (self.w - self.w.mean(axis=0, keepdims=True))**2
        self.hov_t = np.mean(self.t, axis=(-2,-1))
        self.hov_utke = np.mean(utke, axis=(-2,-1))
        self.hov_vtke = np.mean(vtke, axis=(-2,-1))
        self.hov_wtke = np.mean(wtke, axis=(-2,-1))
        self.hov_htke = self.hov_utke + self.hov_vtke
        # Save Enstrophy
        self.outnc.variables['t_hov'][...] = self.hov_t[...,::-1].T
        self.outnc.variables['u_hov'][...] = self.hov_utke[...,::-1].T
        self.outnc.variables['v_hov'][...] = self.hov_vtke[...,::-1].T
        self.outnc.variables['w_hov'][...] = self.hov_wtke[...,::-1].T
        self.outnc.variables['tke_hov'][...] = self.hov_htke[...,::-1].T
        

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
