#
# Downloads the file directly from Github without needing git
#
RUN_DIR=def_JAMES
CFG_DIR=R3_det
EXP_DIR=EXP00

TOKEN="Authorization: token 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"
SOURCE=https://raw.githubusercontent.com/ftucciarone/LocationUncertainty/main/nemo-config

cd /home/${USER}/${RUN_DIR}/cfgs/${CFG_DIR}/MY_SRC
# Noise models
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-NOI/NOI/tlu_dmd.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-NOI/NOI/tlu_flt.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-NOI/NOI/tlu_gss.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-NOI/NOI/tlu_noi.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-NOI/NOI/tlu_pod.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-NOI/NOI/tlu_prj.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-NOI/NOI/tlu_pso.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-NOI/NOI/tlu_rnd.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-NOI/NOI/__future_tlu_wlt%20copie.F90
mv "__future_tlu_wlt%20copie.F90" tlu_wlt.F90

# Base configuration
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/tlu.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/step.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/step_oce.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/nemogcm.F90

# Diagnostic
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/DGN/tlu_dia.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/DGN/tlu_nke.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/DGN/tludgns.F90

# Dynamics
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/DYN/tlu_advDYN.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/DYN/tlu_sbcDYN.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/DYN/tlu_wzvDYN.F90

# New Environments
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/NENV/__future_tlu_MPI.F90
mv __future_tlu_MPI.F90 tlu_MPI.F90

# Stochastic pressure module
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/STP/tlu_STP.F90

# Tracer Advection
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/TRA/tlu_advTRA.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-config/src-TLU/TRA/tlu_sbcTRA.F90


cd /home/${USER}/${RUN_DIR}/cfgs/${CFG_DIR}/${EXP_DIR}
# Download xml files
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/xml-config/domain_def_nemo.xml
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/xml-config/field_def_nemo-oce.xml
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/xml-config/file_def_nemo.xml
