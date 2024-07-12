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

## Download of NEMO
`mkdir ~/nemo_dir` \
`cd ~/nemo_dir` \
`svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0/ .` \

## Compilation of NEMO
`cd ~/nemo_dir` \
`./makenemo -m 'local' -r GYRE_PISCES -n 'Your_Config_Name' -j 32;`\

## Compilation of REBUILD_NEMO
`cd ~/nemo_dir/tools` \
First, replace `EXTERNAL` with `INTRINSIC` in `/tools/REBUILD_NEMO/src/rebuild_nemo.F90`, then launch \
`./maketools -m 'local' -n 'REBUILD_NEMO' -j 32;`
The rebuild procedure needs a `nam_rebuild` file to be read, that has the following structure
```
&nam_rebuild
filebase='GYRE_det_5d_00010101_00151230_grid_T'
ndomain=15
/
```


## Download of this repository
`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`

```
curl -H 'Authorization: token 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' -H 'Accept: application/vnd.github.v3.raw' -O -L https://github.com/ftucciarone/LocationUncertainty
```

## Launch on POPOV
$\color{red}{\text{\textbf{WARNING:}}}$ POPOV does not like to be told the -host option, it messes up the tiling of NEMO for the subdomains 

Example of `./R3_launch.sh` (made executable with `chmod u+x R3_launch.sh`)
```
NPROCS=15

RUN_DIR=def_JAMES
RUN_CONFIG=R3_det
MODE=EXP00

# Choose one, those are different setups for different machines
export LD_LIBRARY_PATH=/home/${USER}/NEMO/installs/lib/
export LD_LIBRARY_PATH=/home/${USER}/nemo-base/lib/

HOME=/home/${USER}

RUN_DIR=${HOME}/${RUN_DIR}/cfgs/${RUN_CONFIG}/${MODE}

cd ${RUN_DIR}
date > t-tic
mpirun -np ${NPROCS} ./nemo
date > t-toc
```
and then launch with detached screen \
```
screen -d -m -S R3_det "./R3_launch.sh"
```
(utility to chech the running is `screen -list`, to kill is `screen -XS <numeric-id> quit`
