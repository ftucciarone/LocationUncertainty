"""Data loading and plots shared by the modredkit_* notebooks."""
import os
import sys

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.animation
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(".."))
from datakit import glorys04
from modredkit.modes import reconstruct as rebuild
from diagkit.fancy_plots import inria_dRdBu_r, seam_mask

COLOR = "#2a6f97"
DT = 1.0   # days
CMAP = inria_dRdBu_r.copy()
CMAP.set_bad((0, 0, 0, 0))   # transparent: masked seam cells span the whole map
GLOBAL_LON0 = -90
plt.rcParams["animation.html"] = "jshtml"
plt.rcParams["animation.embed_limit"] = 60


def gulf_box():
    lat, lon = glorys04.mesh("gphit", "glamt")
    jj, ii = np.where((lat > 30) & (lat < 45) & (lon > -75) & (lon < -45))
    return jj.min(), jj.max()+1, ii.min(), ii.max()+1


def _axes(fig, n, wet):
    """Map axes and pcolormesh coordinates for the global grid or the Gulf Stream box."""
    lat, lon = glorys04.mesh("gphit", "glamt")
    if wet.shape == lat.shape:
        proj = ccrs.EqualEarth(central_longitude=GLOBAL_LON0)
        seam = seam_mask(lon, lat, GLOBAL_LON0)
    else:
        j0, j1, i0, i1 = gulf_box()
        lat, lon = lat[j0:j1, i0:i1], lon[j0:j1, i0:i1]
        proj = ccrs.PlateCarree()
        seam = np.zeros(lat.shape, dtype=bool)
    axes = [fig.add_subplot(1, n, k+1, projection=proj) for k in range(n)]
    for ax in axes:
        ax.set_global() if wet.shape[0] > 500 else ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, facecolor="0.82", zorder=2)
        ax.coastlines(lw=0.4, color="0.35", zorder=3)
    # Project here, so that pcolormesh does not try to unwrap the ORCA grid itself
    pts = proj.transform_points(ccrs.PlateCarree(), lon, lat)
    return axes, pts[..., 0], pts[..., 1], seam


def _mesh(ax, x, y, field, lim, **kwargs):
    """Plain matplotlib pcolormesh: coordinates are already projected, and cartopy's wrapper would rewrap them."""
    ax.autoscale(False)
    return matplotlib.axes.Axes.pcolormesh(ax, x, y, field, cmap=CMAP, vmin=-lim, vmax=lim, shading="auto", **kwargs)


def _field(f, wet, seam):
    m = np.full(wet.shape, np.nan)
    m[wet] = f
    return np.ma.masked_invalid(np.where(seam, np.nan, m))


def _scale(fields):
    return np.percentile(np.abs(np.concatenate([np.ravel(f) for f in fields])), 99.5)


def load_ssh():
    """GLORYS04 SSH, 364 days. Returns anomalies (n, m), area weights and the wet mask."""
    ssh = glorys04.read("sossheig")
    wet = ~np.isnan(ssh).any(axis=0)
    X = ssh[:, wet].T.astype(np.float64)
    e1t, e2t = glorys04.mesh("e1t", "e2t")
    return X - X.mean(axis=1, keepdims=True), (e1t*e2t)[wet], wet


def load_gulfstream():
    """Surface U and V, 731 days, 30-45N, 75-45W. Returns anomalies U, V and their wet masks."""
    box = gulf_box()
    out = []
    for var in ("vozocrtx", "vomecrty"):
        a = glorys04.read(var, box=box)
        wet = ~np.isnan(a).any(axis=0)
        X = a[:, wet].T.astype(np.float64)
        out += [X - X.mean(axis=1, keepdims=True), wet]
    return out[0], out[2], out[1], out[3]


def relative_error(rec, X):
    return np.linalg.norm(rec - X) / np.linalg.norm(X)


def show_maps(fields, titles, wets, units=""):
    """Maps of 1D wet-point fields, side by side, on one diverging scale centred on zero."""
    fig = plt.figure(figsize=(6.2*len(fields), 3.4 if wets[0].shape[0] > 500 else 4), dpi=110)
    axes, lon, lat, seam = _axes(fig, len(fields), wets[0])
    lim = _scale(fields)
    for ax, f, title, wet in zip(axes, fields, titles, wets):
        im = _mesh(ax, lon, lat, _field(f, wet, seam), lim)
        ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=axes, shrink=0.8, label=units)
    plt.show()


def animate(panels, titles, wet, step=6, frames=60, units=""):
    """
    Animation of (n, m) snapshot matrices side by side, on one shared scale.
    Shows every `step`-th day, at most `frames` frames.
    """
    plt.ioff()
    fig = plt.figure(figsize=(5.6*len(panels), 3.2 if wet.shape[0] > 500 else 3.8), dpi=80)
    axes, lon, lat, seam = _axes(fig, len(panels), wet)
    days = np.arange(0, min(p.shape[1] for p in panels), step)[:frames]
    lim = _scale([p[:, days] for p in panels])
    meshes = [_mesh(ax, lon, lat, _field(p[:, 0], wet, seam), lim) for ax, p in zip(axes, panels)]
    fig.colorbar(meshes[0], ax=axes, shrink=0.8, label=units)

    def update(k):
        for ax, mesh, p, title in zip(axes, meshes, panels, titles):
            mesh.set_array(_field(p[:, days[k]], wet, seam).ravel())
            ax.set_title(f"{title}, day {days[k]}", fontsize=10)
        return meshes

    anim = matplotlib.animation.FuncAnimation(fig, update, frames=len(days), interval=150, blit=False)
    html = anim.to_jshtml()
    plt.close(fig)
    plt.ion()
    from IPython.display import HTML
    return HTML(html)


def show_variance(S, title):
    plt.figure(figsize=(7, 3))
    plt.semilogx(np.arange(1, S.size+1), np.cumsum(S**2)/np.sum(S**2), color=COLOR, lw=2)
    plt.xlabel("Number of modes"); plt.ylabel("Cumulative variance"); plt.title(title, fontsize=10)
    plt.grid(alpha=0.25); plt.tight_layout(); plt.show()


def show_spectrum(ct, b, title):
    keep = np.abs(ct.imag) > 1e-8
    plt.figure(figsize=(7, 3))
    plt.loglog(2*np.pi/np.abs(ct.imag[keep]), np.abs(b[keep]), "o", ms=5, color=COLOR)
    plt.xlabel("Period (days)"); plt.ylabel("|amplitude|"); plt.title(title, fontsize=10)
    plt.grid(alpha=0.25); plt.tight_layout(); plt.show()
