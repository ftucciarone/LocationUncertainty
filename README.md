# Transport under Location Uncertainty, a NEMO journey
```
src-config/
├── src-NOI/                   (contains Noise generation routines only)
│   ├── NOI/                     (routines for noise generation)
│   │   ├── tlu_dmd.F90            (Dynamical Mode decomposition)
│   │   ├── tlu_flt.F90            (Filtering routines)
│   │   ├── tlu_gss.F90            (Vertical profile parametrization)
│   │   ├── tlu_noi.F90            (General noise generation routines)
│   │   ├── tlu_pod.F90            (Proper Orthogonal decomposition)
│   │   ├── tlu_prj.F90            (Isopycnal projection)
│   │   ├── tlu_pso.F90            (Pseudo Observation POD)
│   │   ├── tlu_rnd.F90            (Random number generation)
│   │   └── tlu_wlt.F90            (Wavelet noise)
│   ├── USR/                     (usrdef routines to hijack NEMO)
│   │   └── tlu_usrdef_nam.F90     (domain modifications for wlt noise)
│   ├── nemogcm.F90              (Main file)
│   ├── step_oce.F90             (Dependencies)
│   ├── step.F90                 (Time stepping with call to noise generation only)
│   └── tlu.F90                  (Main initialization routine)
├── src-TLU/                   (contains Noise application routines only)
│   ├── DGN/                     (Diagnostics)
│   │   ├── tlu_dia.F90
│   │   └── tlu_nke.F90
│   ├── DYN/                     (Dynamics) 
│   │   ├── tlu_advDYN.F90         (Dynamical advection)
│   │   ├── tlu_sbcDYN.F90         (Dynamical surface BC)
│   │   └── tlu_wzvDYN.F90         (Vertical velocity)
│   ├── NENV/                    (Various)
│   │   └── tlu_MPI.F90
│   ├── STP/                     (Stochastic pressure modules)
│   │   └── tlu_STP.F90            (Stochastic pressure)
│   ├── TO_EXP/                  (Experiment files)
│   │   ├── namelist_cfg
│   │   └── namelist_ref
│   ├── TRA/                     (Tracer application)
│   │   ├── tlu_advTRA.F90.        (Tracer advection)
│   │   └── tlu_sbcTRA.F90.        (Tracer surface BC)
│   ├── nemogcm.F90              (Main file)
│   ├── step_oce.F90             (Dependencies)
│   ├── step.F90                 (Time stepping with generation and application)
│   └── tlu.F90                  (Main initialization routine)
└── xml-config
    ├── domain_def_nemo.xml
    ├── field_def_nemo.xml
    └── file_def_nemo.xml
```
## Docker


## Compilation of NEMO
 ./makenemo -m 'local' -r GYRE_PISCES -n 'Your_Config_Name' -j 32;
## Compilation of REBUILD_NEMO
First, replace `EXTERNAL` with `INTRINSIC` in `/tools/REBUILD_NEMO/src/rebuild_nemo.F90`, then launch 
 cd tools; ./maketools -m 'local' -n 'REBUILD_NEMO' -j 32;



 
