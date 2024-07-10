TOKEN="Authorization: token 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"
SOURCE=https://raw.githubusercontent.com/ftucciarone/LocationUncertainty/main/nemo-config/src-config/

# Noise models
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-NOI/NOI/tlu_dmd.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-NOI/NOI/tlu_flt.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-NOI/NOI/tlu_gss.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-NOI/NOI/tlu_noi.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-NOI/NOI/tlu_pod.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-NOI/NOI/tlu_prj.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-NOI/NOI/tlu_pso.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-NOI/NOI/tlu_rnd.F90

# Base configuration
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/tlu.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/step.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/step_oce.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/nemogcm.F90

# Diagnostic
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/DGN/tlu_dia.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/DGN/tlu_nke.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/DGN/tludgns.F90

# Dynamics
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/DYN/tlu_advDYN.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/DYN/tlu_sbcDYN.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/DYN/tlu_wzvDYN.F90

# New Environments
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/NENV/__future_tlu_MPI.F90
mv __future_tlu_MPI.F90 tlu_MPI.F90

# Stochastic pressure module
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/STP/tlu_STP.F90

# Tracer Advection
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/TRA/tlu_advTRA.F90
curl -H "${TOKEN}" -H 'Accept: application/vnd.github.v3.raw' -O -L ${SOURCE}/src-TLU/TRA/tlu_sbcTRA.F90
