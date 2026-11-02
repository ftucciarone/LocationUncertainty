#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module fancy_plots

1. Project (INRIA) custom diverging colormaps -- same blue/white/red
   structure as matplotlib's RdBu/RdBu_r, but with project-specific endpoint
   colors instead of ColorBrewer's.

2. The AREAS OF INTEREST of the eORCA1 analysis, written in stone
   (ZONES, AGRIF_BOX, ZONE_COLOR, ZONE_BOX_EORCA1), and the region LOCATOR
   built on them (plot_zone_locator): a minimal Equal Earth land/ocean map
   with one area outlined, small enough to be dropped into a corner of any
   other figure.  cartopy is imported lazily, inside the functions that
   need it, so importing this module for the colormaps alone does not need
   cartopy.
"""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

white = np.array([1.0, 1.0, 1.0])
inria_ldRed = np.array([163, 42, 24]) / 255
inria_dBlue = np.array([20, 30, 55]) / 255
inria_lBlue = np.array([57, 66, 85]) / 255

# low=red, high=blue (matplotlib's "RdBu" convention)
inria_dRdBu = LinearSegmentedColormap.from_list("inria_dRdBu", [inria_ldRed, white, inria_dBlue], N=256)
inria_lRdBu = LinearSegmentedColormap.from_list("inria_lRdBu", [inria_ldRed, white, inria_lBlue], N=256)
# low=blue, high=red (matplotlib's "RdBu_r" convention)
inria_dRdBu_r = LinearSegmentedColormap.from_list("inria_dRdBu_r", [inria_dBlue, white, inria_ldRed], N=256)
inria_lRdBu_r = LinearSegmentedColormap.from_list("inria_lRdBu_r", [inria_lBlue, white, inria_ldRed], N=256)


# =============================================================================
# Areas of interest -- WRITTEN IN STONE
# =============================================================================
# The five areas used throughout the analysis (defined and checked in
# the earlier analysis notebook, section 5).  Do not "re-derive" them in a notebook:
# import them from here, so every figure and every reader uses the same boxes.
#
# Four lat/lon boxes -- {name: ((lat_min, lat_max), (lon_min, lon_max))}, in
# degrees.  Nino3's south bound is -8.0, not the -10 of the older notebooks:
# the AGRIF nest's north edge sits at -8.51 deg inside Nino3's longitude
# window, so the old -10 overlapped the nest; -8.0 is that edge plus a 0.5 deg
# margin, rounded up to the nearest 0.5 deg.  Only the south bound moved.
ZONES = {
    "ACC":        ((-72, -36), (-180, 180)),   # Southern Ocean circumpolar band
    "GulfStream": ((24, 49),   (-85, -20)),    # western North Atlantic -> mid-Atlantic
    "Kuroshio":   ((20, 45),   (118, 180)),    # western North Pacific -> dateline
    "Nino3":      ((-8.0, 10), (-150, -80)),   # eastern equatorial Pacific, clear of the AGRIF nest
}

# The fifth area, the AGRIF Pacific nest footprint, as a PARENT-grid index box
# (j0, j1, i0, i1), 0-based half-open: AGRIF_FixedGrids.in says
# "108 208 65 165" (1-based imin imax jmin jmax) for every nest built so far
# (rho = 3, 4 and 5), i.e. i[107:208], j[64:165].
AGRIF_BOX = (64, 165, 107, 208)

# One colour per area, used by every figure that draws them.
ZONE_COLOR = {"ACC": "C0",
              "GulfStream": "C1",
              "Kuroshio": "C2",
              "Nino3": "C3",
              "AGRIFZoom": "purple"}

# The index boxes (j0, j1, i0, i1) that latlon_box gives for ZONES on the
# eORCA1 parent grid (331 x 360), frozen as a tripwire: zone_boxes() asserts
# it still gets exactly these on that grid.
ZONE_BOX_EORCA1 = {
    "ACC":        (65, 133, 0, 360),
    "GulfStream": (226, 260, 202, 266),
    "Kuroshio":   (222, 256, 45, 107),
    "Nino3":      (165, 211, 137, 207),
    "AGRIFZoom":  AGRIF_BOX,
}


def latlon_box(nav_lat, nav_lon, lat_bnds, lon_bnds):
    r"""
    latlon_box: the index box of the grid points inside a latitude-longitude rectangle.

    On a curvilinear grid a region given in degrees is not an index range; this finds the smallest box
    of indices holding all the points inside it. The upper bounds are moved down so each side has an even
    length, which spectral routines need.

        Parameters:
            nav_lat  (float): (y, x) latitudes, in degrees
            nav_lon  (float): (y, x) longitudes, in degrees
            lat_bnds (tuple): (south, north), exclusive
            lon_bnds (tuple): (west, east), exclusive

        Returns:
            box (tuple): (j0, j1, i0, i1)

        Formula:
            (j_0, j_1, i_0, i_1) = \text{bounding box of } \{(j,i) : \varphi_{1} < \varphi_{j,i} < \varphi_{2},\ \lambda_{1} < \lambda_{j,i} < \lambda_{2}\}

        See also:
            zone_boxes
    """
    m = ((nav_lat > lat_bnds[0]) & (nav_lat < lat_bnds[1]) &
         (nav_lon > lon_bnds[0]) & (nav_lon < lon_bnds[1]))
    j, i = np.where(m)
    j0, i0 = int(j.min()), int(i.min())
    j1, i1 = int(j.max()) + 1, int(i.max()) + 1
    return j0, j1 - (j1 - j0) % 2, i0, i1 - (i1 - i0) % 2


def zone_boxes(nav_lat, nav_lon):
    r"""
    zone_boxes: the index boxes of the four study zones and of the AGRIF nest.

    On the eORCA1 parent grid the result is checked against frozen values, so a change to the zones or a
    different grid file cannot change the regions studied without notice.

        Parameters:
            nav_lat (float): (y, x) latitudes, in degrees
            nav_lon (float): (y, x) longitudes, in degrees

        Returns:
            boxes (dict): {name: (j0, j1, i0, i1)}, with 'AGRIFZoom' for AGRIF_BOX. Raises
                          AssertionError if they differ from the frozen eORCA1 boxes

        Formula:
            \{\text{zone}: \text{latlon\_box}(\text{zone bounds})\} \cup \{\text{AGRIFZoom}: \text{AGRIF\_BOX}\}

        See also:
            latlon_box, plot_zone_locator
    """
    boxes = {name: latlon_box(nav_lat, nav_lon, lat_b, lon_b) for name, (lat_b, lon_b) in ZONES.items()}
    boxes["AGRIFZoom"] = AGRIF_BOX
    if nav_lat.shape == (331, 360):
        assert boxes == ZONE_BOX_EORCA1, (f"zone boxes changed on the eORCA1 grid: {boxes} "
                                          f"!= frozen {ZONE_BOX_EORCA1}")
    return boxes


def seam_mask(nav_lon, nav_lat, central_longitude, factor=8):
    r"""
    seam_mask: the cells to hide so that pcolormesh does not stretch quads across the map.

    On eORCA1 the longitude wraps along a cut that runs diagonally through index space, and the north
    fold adds another. Plotted naively, cells across the cut become long bands over the whole map. Each
    point is projected, and cells whose jump to an index neighbour is far above the median spacing are
    masked, so it works for any central longitude.

        Parameters:
            nav_lon           (float): (y, x) longitudes, in degrees
            nav_lat           (float): (y, x) latitudes, in degrees
            central_longitude (float): Of ccrs.EqualEarth
            factor            (float): Jumps above factor times the median spacing count as a cut

        Returns:
            mask (bool): (y, x) True on cells touching a cut

        Formula:
            \text{mask}_{j,i} = \|\mathbf{p}_{j,i} - \mathbf{p}_{\text{nb}}\| > \text{factor} \cdot \mathrm{median}\|\Delta \mathbf{p}\| \text{ for any index neighbour}

        See also:
            plot_zone_locator
    """
    import cartopy.crs as ccrs
    proj = ccrs.EqualEarth(central_longitude=central_longitude)
    pts = proj.transform_points(ccrs.PlateCarree(), nav_lon, nav_lat)
    x, y = pts[..., 0], pts[..., 1]

    dist_i = np.hypot(np.diff(x, axis=1), np.diff(y, axis=1))  # (y, x-1): i -> i+1
    dist_j = np.hypot(np.diff(x, axis=0), np.diff(y, axis=0))  # (y-1, x): j -> j+1
    typical = np.median(np.concatenate([dist_i.ravel(), dist_j.ravel()]))

    bad = np.zeros(x.shape, dtype=bool)
    bad_i = dist_i > factor * typical
    bad[:, :-1] |= bad_i
    bad[:, 1:] |= bad_i
    bad_j = dist_j > factor * typical
    bad[:-1, :] |= bad_j
    bad[1:, :] |= bad_j
    return bad


def zone_outline_lonlat(zname, nav_lon, nav_lat, zone_box):
    r"""
    zone_outline_lonlat: polylines outlining a zone, along its grid-index perimeter.

    On a curvilinear grid the index perimeter of a zone and an ideal latitude-longitude rectangle differ
    visibly; the outline drawn is the region actually studied. The Antarctic Circumpolar Current zone is
    drawn as two circles at its latitude bounds, since its rows cross the index-space seam.

        Parameters:
            zname     (str): Key of zone_box
            nav_lon (float): (y, x) longitudes, in degrees
            nav_lat (float): (y, x) latitudes, in degrees
            zone_box (dict): From zone_boxes

        Returns:
            lines (list): (lons, lats) arrays, one per side

        Formula:
            \big(\lambda, \varphi\big) \text{ along rows } j_0, j_1 - 1 \text{ and columns } i_0, i_1 - 1

        See also:
            plot_zone_locator, zone_boxes
    """
    j0, j1, i0, i1 = zone_box[zname]
    if zname == "ACC":
        lat_b, _ = ZONES[zname]
        lons_ring = np.linspace(-180, 180, 361)
        return [(lons_ring, np.full_like(lons_ring, lat_b[0])),
                (lons_ring, np.full_like(lons_ring, lat_b[1]))]
    south = (nav_lon[j0, i0:i1], nav_lat[j0, i0:i1])
    north = (nav_lon[j1 - 1, i0:i1], nav_lat[j1 - 1, i0:i1])
    west = (nav_lon[j0:j1, i0], nav_lat[j0:j1, i0])
    east = (nav_lon[j0:j1, i1 - 1], nav_lat[j0:j1, i1 - 1])
    return [south, north, west, east]


# =============================================================================
# Region locator
# =============================================================================

def plot_zone_locator(zname, nav_lon, nav_lat, land, zone_box, ax=None, central_longitude=-90,
                      box_color="red", land_color="0.15", show_agrif=False, agrif_color="purple",
                      linewidth=1.5):
    r"""
    plot_zone_locator: a small map showing where a zone is: the land silhouette with the zone outlined.

    Meant as an inset beside a regional figure, so it has no title, ticks, gridlines or colorbar.

        Parameters:
            zname    (str): 'ACC', 'GulfStream', 'Kuroshio', 'Nino3' or 'AGRIFZoom'
            nav_lon (float): (y, x) longitudes of the parent grid, in degrees
            nav_lat (float): (y, x) latitudes of the parent grid, in degrees
            land     (bool): (y, x) True on land, e.g. bathy <= 0
            zone_box (dict): From zone_boxes
            ax    (GeoAxes): Axes with projection ccrs.EqualEarth(central_longitude=...), or None
                             for a new small figure
            central_longitude (float): Of the projection, when ax is None
            box_color    (str): Colour of the zone outline
            land_color   (str): Colour of land
            show_agrif  (bool): Also outline the AGRIF nest, unless zname is 'AGRIFZoom'
            agrif_color  (str): Colour of the AGRIF outline
            linewidth  (float): Width of the outlines

        Returns:
            ax (GeoAxes): The axes drawn on

        Formula:
            \text{land}(\lambda, \varphi) \text{ on Equal Earth} + \text{zone\_outline\_lonlat}

        See also:
            zone_outline_lonlat, seam_mask
    """
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    own_fig = ax is None
    if own_fig:
        fig = plt.figure(figsize=(2.4, 1.3))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.EqualEarth(central_longitude=central_longitude))

    # Plain land/ocean silhouette, not bathymetry: a locator only has to be
    # recognisable as "the world map".  Seam-masked like every other map.
    seam = seam_mask(nav_lon, nav_lat, central_longitude)
    silhouette = np.ma.masked_where(seam, np.asarray(land).astype(float))   # 0 ocean / 1 land
    ax.pcolormesh(nav_lon, nav_lat, silhouette, transform=ccrs.PlateCarree(),
                  cmap=ListedColormap(["white", land_color]), vmin=0, vmax=1, shading="auto")
    ax.set_global()

    if show_agrif and zname != "AGRIFZoom":
        for lons_e, lats_e in zone_outline_lonlat("AGRIFZoom", nav_lon, nav_lat, zone_box):
            ax.plot(lons_e, lats_e, transform=ccrs.PlateCarree(), color=agrif_color, linewidth=linewidth)
    for lons_e, lats_e in zone_outline_lonlat(zname, nav_lon, nav_lat, zone_box):
        ax.plot(lons_e, lats_e, transform=ccrs.PlateCarree(), color=box_color, linewidth=linewidth)

    for spine in ax.spines.values():
        spine.set_linewidth(0.6)

    if own_fig:
        plt.tight_layout(pad=0.1)
        plt.show()
    return ax

