"""Reads real files from storage. Run: python tests/test_datakit.py"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import datakit
from datakit import cmems, glorys04, glorys12

BOX = (500, 520, 300, 330)
g = glorys04.read("sossheig", "2015-01-01", "2015-01-03", box=BOX)
assert g.shape == (3, 20, 30), g.shape
assert len(glorys04.files("sossheig")) == 364
assert len(glorys04.files("vozocrtx", "2015-02")) == 731 - 31
assert len(glorys04.files("vozocrtx", "2015-02", "2015-02")) == 28
assert len(glorys04.files("vozocrtx", "2016")) == 366
assert glorys04.read("vozocrtx", "2015-01-01", "2015-01-01", box=BOX).shape == (1, 20, 30)
assert glorys04.read("vozocrtx", "2015-01-01", "2015-01-01", level=None, box=BOX).shape == (1, 75, 20, 30)
assert glorys04.read("vozocrtx", "2015-01-01", "2015-01-01", box=BOX, baroclinic=True).shape == (1, 20, 30)
assert glorys04.mesh("e1t").shape == (1021, 1442)
assert len(cmems.files("sossheig")) == 300
assert len(cmems.files("vozocrtx", "2015", "2015")) == 12
assert cmems.read("sossheig", "2015-01", "2015-02", box=BOX).shape == (2, 20, 30)
assert cmems.read("vomecrty", "2015-01", "2015-01", box=BOX).shape == (1, 20, 30)
assert cmems.mesh("e1t").shape == (1021, 1442)
assert glorys12.read("sossheig", "2015-01", "2015-01", box=BOX).shape == (1, 20, 30)
assert glorys12.read("vozocrtx", "1993-01", "1993-01", box=BOX).shape == (1, 20, 30)
assert len(glorys12.files("votemper", "1993", "1993")) == 12
assert glorys12.mesh("e1t").shape == (3059, 4322)
full = glorys04.read("sossheig", "2015-01-01", "2015-01-01")
print("SSH shape", full.shape, "land fraction", np.isnan(full).mean(), "zero fraction", (full == 0).mean(), "range", np.nanmin(full), np.nanmax(full))
u = cmems.read("vozocrtx", "2015-01", "2015-01")
print("CMEMS u land fraction", np.isnan(u).mean(), "zero fraction", (u == 0).mean())
got = list(glorys04.stream("vozocrtx", "vomecrty", start="2015-01-01", end="2015-01-03", box=BOX))
assert [g[0] for g in got] == ["2015-01-01", "2015-01-02", "2015-01-03"] and got[0][1].shape == (20, 30)
direct = glorys04.read("vomecrty", "2015-01-01", "2015-01-03", box=BOX)
assert np.array_equal(np.stack([g[2] for g in got]), direct, equal_nan=True)
assert len(list(cmems.stream("sossheig", start="2015-01", end="2015-02", box=BOX))) == 2
assert next(glorys12.stream("sossheig", box=BOX))[1].shape == (20, 30)
for name, levels in (("glorys04", 75), ("glorys12", 50), ("cmems", 75)):
    src = datakit.source(name)
    groups, missing = src.groups("u", "v", start="2015", end="2015") if name != "glorys12" else src.groups("u", "v", start="1993", end="1993")
    assert groups and not missing
    assert len(src.depth(groups[0][1])) == levels
    assert src.time_hours(groups[0][1]) > 0
    assert src.read("u", groups[0][0], groups[0][0], box=BOX).shape == (1, 20, 30)
    assert src.grid()["boundary"] == "t-fold"
# PP3 on the real grid: eta = H sin^2(lat) gives a purely eastward flow, -g H cos(lat) / (omega R)
from scipy.ndimage import binary_erosion

from orcakit import grid_c, orca
for name in ("cmems", "glorys12"):
    src = datakit.source(name)
    g = src.grid()
    tmask, umask, vmask, gphit, gphiu, gphiv = src.mesh("tmask", "umask", "vmask", "gphit", "gphiu", "gphiv")
    H, R, om = 1.0, 6371e3, 7.292116e-5
    Ug, Vg = grid_c.geostrophic_uv(H*np.sin(np.radians(gphit))**2, tmask[0].astype(bool), umask[0].astype(bool), vmask[0].astype(bool),
                                  g["e1t"], g["e2t"], g["e1v"], g["e2u"], gphiu, gphiv, orca.halo)
    # Open water only: the wet-only average is biased within two cells of the coast. South of 20N the grid follows latitude
    band = (gphiu < -10) & (gphiu > -40) & binary_erosion(tmask[0].astype(bool), iterations=2)
    exact = -9.80665*H*np.cos(np.radians(gphiu))/(om*R)
    err = np.nanmax(np.abs(Ug - exact)[band])/np.nanmax(np.abs(exact)[band])
    print(name, "geostrophic relative error on the real grid, open water 10-40S:", err)
    assert err < 1e-3

print("ALL OK")
