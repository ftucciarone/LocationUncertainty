#!/bin/bash
EXEDIR=/srv/tempdd/${USER}/nemo2/release-4.0/tools/REBUILD_NEMO


cp nam_rebuild_restart nam_rebuild
$EXEDIR/rebuild_nemo.exe

sed -i s/restart/restart_trc/g nam_rebuild
$EXEDIR/rebuild_nemo.exe

cp nam_rebuild_state nam_rebuild
$EXEDIR/rebuild_nemo.exe

sed -i s/grid_U/grid_V/g nam_rebuild
$EXEDIR/rebuild_nemo.exe

sed -i s/grid_V/grid_W/g nam_rebuild
$EXEDIR/rebuild_nemo.exe

sed -i s/grid_W/grid_T/g nam_rebuild
$EXEDIR/rebuild_nemo.exe

sed -i s/grid_T/trendT/g nam_rebuild
$EXEDIR/rebuild_nemo.exe

sed -i s/trendT/trenda/g nam_rebuild
$EXEDIR/rebuild_nemo.exe

sed -i s/trenda/tlu_T_once/g nam_rebuild
$EXEDIR/rebuild_nemo.exe

sed -i s/tlu_T_once/grid_T/g nam_rebuild
sed -i s/5d/1y/g nam_rebuild
$EXEDIR/rebuild_nemo.exe

sed -i s/1y/1d/g nam_rebuild
$EXEDIR/rebuild_nemo.exe

FILE=mesh_mask_0000.nc
if [ -f "$FILE" ]; then
    $EXEDIR/rebuild_nemo.exe nam_rebuild_mask 
    sed -i s/mesh_mask/domain_cfg_out/g nam_rebuild_mask

    $EXEDIR/rebuild_nemo.exe nam_rebuild_mask
    sed -i s/domain_cfg_out/mesh_mask/g nam_rebuild_mask
fi
