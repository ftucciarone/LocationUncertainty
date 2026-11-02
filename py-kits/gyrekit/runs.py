#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module runs

Find and read GYRE experiments from config/paths_R3.yaml and paths_R9.yaml.
Replaces OpenYAML.yaml2dict of the earlier analysis code. The configs point to the
result folders set in paths.toml (gyre_root).
"""
import os

import netCDF4
import numpy as np
import yaml

from orcakit import paths

CONFIG = os.path.join(os.path.dirname(__file__), "config")
KINDS = ("curlfile", "tempfile", "uvelfile", "vvelfile", "wvelfile")


def experiments(resolution):
    r"""
    experiments: the GYRE experiments of a resolution, with their files.

    Paths come from paths_<resolution>.yaml, with the high-resolution runs included; they are not checked.

        Parameters:
            resolution (str): 'R3' or 'R9'

        Returns:
            exps (dict): {name: {'label', 'color', 'files': {kind: [paths]}}}

        Formula:
            \{\text{name}: (\text{label}, \text{color}, \{\text{kind}: [\text{paths}]\})\}

        See also:
            gridfile, read_field
    """
    with open(os.path.join(CONFIG, f"paths_{resolution}.yaml")) as f:
        cfg = yaml.safe_load(f)
    root = paths.get("gyre_root") + cfg["data_dir"]
    out = {}
    for exp in {**cfg["HRexps"],
                **cfg["LRexps"]}.values():
        base = root + exp["base_dir"]
        files = {k: [f"{base}{sub}/{exp[k]}{cfg['ext']}".replace("//", "/") for sub in exp["subs_dir"]] for k in KINDS}
        out[exp["name"]] = {"label": exp["label"],
                            "color": exp["color"],
                            "files": files}
    return out


def gridfile(resolution):
    r"""
    gridfile: the path of the domain_cfg file of a GYRE resolution.

        Parameters:
            resolution (str): 'R3', 'R9' or 'R27'

        Returns:
            path (str): domain_cfg_<resolution>.nc in the configured folder

        Formula:
            \text{gyre\_root}/\text{dom\_dir}/\text{domain\_cfg\_}R.nc

        See also:
            read_spacing
    """
    with open(os.path.join(CONFIG, "paths_R3.yaml")) as f:
        cfg = yaml.safe_load(f)
    return f"{paths.get('gyre_root')}{cfg['dom_dir']}/domain_cfg_{resolution}.nc"


def read_field(paths, var, levels=None, dtype="float32"):
    r"""
    read_field: read a variable from several files, joined along time.

        Parameters:
            paths  (list): Files, in time order
            var     (str): Variable
            levels  (int): Vertical levels of 4D fields, index or slice; None for all
            dtype   (str): Type of the result

        Returns:
            field (float): Variable concatenated along its first axis

        Formula:
            X = \big[\, f_1[v];\ f_2[v];\ \dots \,\big] \text{ along time}

        See also:
            experiments
    """
    out = []
    for p in paths:
        with netCDF4.Dataset(p) as nc:
            v = nc[var]
            data = v[:, levels] if (levels is not None and v.ndim == 4) else v[:]
            out.append(np.ma.getdata(data).astype(dtype))
    return np.concatenate(out, axis=0)


def read_spacing(gridfile):
    r"""
    read_spacing: the uniform grid spacing of a GYRE domain_cfg file.

    GYRE grids are regular, so the first cell's widths hold everywhere.

        Parameters:
            gridfile (str): domain_cfg file

        Returns:
            dx (float): e1t at the first point, in metres
            dy (float): e2t at the first point, in metres

        Formula:
            (\Delta x, \Delta y) = (e_{1t}, e_{2t})_{0,0}

        See also:
            gridfile
    """
    with netCDF4.Dataset(gridfile) as nc:
        return float(nc["e1t"][0, 0, 0]), float(nc["e2t"][0, 0, 0])
