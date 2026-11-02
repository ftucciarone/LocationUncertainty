"""Synthetic-data checks. Run: python tests/test_kits.py"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from modredkit import dmd, modes, pod
from noisekit import correlation, noise, running, structure, welford  # noqa: F401
from orcakit import orca, grid_c, netcdf, filters  # noqa: F401
from diagkit import spectra, conservation, fancy_plots, metrics, flow, stats  # noqa: F401
from gyrekit import grids, runs, diag

DT = 1.0
OMEGA = np.array([0.3, 0.7])


def field(nt=80, n=200):
    x = np.linspace(0, 2*np.pi, n)
    t = np.arange(nt) * DT
    # Two travelling waves: rank 4, so two conjugate eigenvalue pairs
    return (np.outer(np.sin(x), np.cos(OMEGA[0]*t)) + np.outer(np.cos(x), np.sin(OMEGA[0]*t))
            + 0.5*(np.outer(np.sin(2*x), np.cos(OMEGA[1]*t)) + np.outer(np.cos(2*x), np.sin(OMEGA[1]*t))))


def test_pod_reconstructs():
    X = field() + 3.0
    m, S, c, mean = pod.pod(X, r=4)
    assert np.allclose(pod.reconstruct(m, c, mean), X, atol=1e-9)
    assert np.allclose(m.T @ m, np.eye(4), atol=1e-9)


def test_pod_weights():
    X = field()
    w = np.linspace(1, 3, X.shape[0])
    m, S, c, mean = pod.pod(X, r=4, weights=w)
    assert np.allclose(m.T @ (w[:, None]*m), np.eye(4), atol=1e-9)
    assert np.allclose(pod.reconstruct(m, c, mean), X, atol=1e-9)


def test_rank_from_energy():
    assert pod.rank_from_energy(np.array([3., 1., .1]), .85) == 1
    assert pod.rank_from_energy(np.array([3., 1., .1]), .995) == 2


def _freqs(ct):
    big = ct[np.argsort(-np.abs(ct.imag))]
    return np.sort(np.unique(np.round(np.abs(big.imag), 6)))


def test_exact_dmd_frequencies():
    X = field()
    mds, dt_eig, ns = dmd.exactDMD(X[:, :-1], X[:, 1:], 4)[:3]
    ct = dmd.DMD_continuous_eig(dt_eig, DT)
    assert np.allclose(_freqs(ct), OMEGA, atol=1e-6), _freqs(ct)


def test_kernel_dmd_frequencies():
    X = field()
    out = dmd.kernelDMD(X.T @ X, X[:, :-1].T, DT, 4)
    assert np.allclose(_freqs(out[1]), OMEGA, atol=1e-6), _freqs(out[1])
    assert out[0].shape == (X.shape[0], 4)


def test_optimal_amplitudes_reconstruct():
    X = field()
    mds, dt_eig, ns = dmd.exactDMD(X[:, :-1], X[:, 1:], 4)[:3]
    b = dmd.DMD_amplitudes_optimal(mds, dt_eig, X)
    rec = (mds * b) @ (dt_eig[:, None] ** np.arange(X.shape[1])[None, :])
    assert np.linalg.norm(rec.real - X) / np.linalg.norm(X) < 1e-6


def test_pod_multi_matches_concatenated():
    X = field(); Y = 2*field()[::-1] + np.outer(np.ones(200), np.sin(0.2*np.arange(80)))
    w1, w2 = np.linspace(1, 3, 200), np.linspace(2, 4, 200)
    ms, S, c, mu = pod.pod_multi([X, Y], weights=[w1, w2])
    _, S0, _, _ = pod.pod(np.vstack((X, Y)), weights=np.r_[w1, w2], r=S.size)
    assert np.allclose(S, S0, rtol=1e-6)
    G = ms[0].T @ (w1[:, None]*ms[0]) + ms[1].T @ (w2[:, None]*ms[1])
    assert np.allclose(G, np.eye(S.size), atol=1e-6)
    assert np.allclose(pod.reconstruct(ms[0], c, mu[0]), X, atol=1e-6)
    assert np.allclose(pod.reconstruct(ms[1], c, mu[1]), Y, atol=1e-6)


def test_pod_multi_full_weight_matrix():
    X = field(); w = np.linspace(1, 3, 200)
    a = pod.pod_multi([X], weights=[w])[1]
    b = pod.pod_multi([X], weights=[np.diag(w)])[1]
    assert np.allclose(a, b)


def test_pod_multi_drops_null_modes():
    assert pod.pod_multi([field()])[1].size == 4


def test_mrdmd():
    X = field(nt=160)
    mds, dtE, ct, amp, idx = dmd.mrDMD(X, 4, DT, 3)
    assert idx.shape == (ct.size, 3) and mds.shape[1] == ct.size
    assert sorted(set(map(tuple, idx[:, :2].tolist()))) == [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2), (2, 3)]
    assert np.all(np.abs(dtE) <= 1 + 1e-6)
    assert np.allclose(_freqs(ct[idx[:, 0] == 2]), OMEGA, atol=1e-5), _freqs(ct[idx[:, 0] == 2])


def test_mrdmd_keeps_damped_modes():
    X = field(nt=160) * np.exp(-0.02*np.arange(160))
    ct = dmd.mrDMD(X, 4, DT, 1)[2]
    assert ct.size == 4 and np.allclose(ct.real, -0.02, atol=1e-5), ct


def test_svd_rank():
    X = field()
    assert dmd.SVD(X, 4)[1].size == 4
    assert dmd.SVD(X, None)[1].size == X.shape[1]
    full = np.random.default_rng(0).normal(size=(50, 30))
    assert dmd.SVD(full, 100.0)[1].size == 30   # a percentage, not a rank
    assert dmd.SVD(X, 0.5)[1].size <= 4


def test_modes():
    ct = np.array([1j*OMEGA[0], -1j*OMEGA[0], 1j*OMEGA[1], -1j*OMEGA[1], -0.5+0j])
    eig = np.exp(ct*DT)
    assert modes.mask_one_per_pair(ct).sum() == 3
    assert modes.mask_unit_circle(eig).sum() == 4
    assert modes.mask_not_growing(np.array([1.5, 0.5, 1.0])).tolist() == [False, True, True]
    assert np.allclose(np.abs(ct[modes.mask_period(ct, 15., 30.)].imag), OMEGA[0])
    mds = np.eye(5, dtype=complex)
    got, a = modes.apply_mask(modes.mask_top_energy(mds, np.arange(5.)+1, 2), mds, np.arange(5.))
    assert got.shape == (5, 2) and a.tolist() == [3., 4.]
    assert np.allclose(modes.time_dynamics(eig, 5)[:, 3], np.exp(ct*3))
    assert np.allclose(np.abs(modes.time_dynamics(eig, 5, "undamped")), 1)
    X = field()
    m, e, _ = dmd.exactDMD(X[:, :-1], X[:, 1:], 4)
    b = dmd.DMD_amplitudes_optimal(m, e, X)
    assert np.linalg.norm(modes.reconstruct(m, b, e, X.shape[1]) - X) / np.linalg.norm(X) < 1e-6
    assert np.allclose(modes.gram(X), X.T @ X)
    assert np.isclose(modes.gram(X, "poly", 2, 1.0)[0, 0], (1 + (X.T @ X)[0, 0])**2)


def test_pack_unpack():
    rng = np.random.default_rng(1)
    u, v = rng.normal(size=(7, 4, 5)), rng.normal(size=(7, 4, 5))
    wu, wv = rng.random((4, 5)) > 0.3, rng.random((4, 5)) > 0.3
    X = modes.pack([u, v], [wu, wv])
    assert X.shape == (wu.sum() + wv.sum(), 7)
    a, b = modes.unpack(X, [wu, wv])
    assert np.array_equal(a[:, wu], u[:, wu]) and np.isnan(a[:, ~wu]).all() and np.array_equal(b[:, wv], v[:, wv])


def test_rolling_stats():
    snaps = [np.full((2, 2), float(k)) for k in range(5)]
    got = list(running.rolling_stats(snaps, 3))
    assert [g[0] for g in got] == [2, 3, 4]
    assert np.allclose(got[0][1], 1) and np.allclose(got[0][2], np.std([0, 1, 2]))


def test_downsample_matches_loop():
    rng = np.random.default_rng(2)
    f = rng.normal(size=(3, 20, 30)) + 1j*rng.normal(size=(3, 20, 30))
    mask = rng.random((4, 5)) > 0.4
    ilat, ilon = rng.integers(0, 20, (4, 5)), rng.integers(0, 30, (4, 5))
    out = orca.downsample_field(f, ilat, ilon, mask)
    loop = np.zeros((3, 4, 5), dtype=complex)
    for j in range(4):
        for i in range(5):
            if mask[j, i]:
                loop[:, j, i] = f[:, ilat[j, i], ilon[j, i]]
    assert np.array_equal(out, loop) and orca.downsample_field(f[0].real, ilat, ilon, mask).shape == (4, 5)


def test_coastal_taper_and_modes_file():
    import tempfile
    import netCDF4
    m = np.ones((3, 20, 20)); m[:, :5] = 0
    t = orca.coastal_taper(m, 1.5, 4)
    assert t[:, :5].max() == 0 and 0 < t[1, 6, 10] < 1 and t[1, 18, 10] > t[1, 6, 10]
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "m.nc")
        netcdf.write_modes(p, {"uavg": np.ones((2, 3, 4))}, {"umode": (1+2j)*np.ones((5, 2, 3, 4))}, np.arange(5.))
        with netCDF4.Dataset(p) as nc:
            assert nc["umode_imag_r"].shape == (5, 2, 3, 4) and np.allclose(nc["umode_imag_r"][:], 2)
            assert np.allclose(nc["omega"][:], np.arange(5))


def test_geostrophy_shape_and_sign():
    ny, nx = 12, 14
    one = np.ones((ny, nx))
    ssh = np.tile(np.arange(ny, dtype=float), (nx, 1)).T
    ug, vg = noise.uv_from_ssh_dA(ssh, one, one, 1e-4*one, one, one)
    assert np.allclose(ug, -9.806/1e-4)
    assert np.allclose(vg, 0)


def test_metrics():
    w = np.ones((4, 5))
    a = np.arange(20.).reshape(4, 5) + 1
    assert np.isclose(metrics.pattern_correlation(a, 3*a, w), 1)
    assert np.isclose(metrics.pattern_correlation(a, -a, w), -1)
    assert np.isclose(metrics.rmse(a+2, a, w), 2)
    assert np.isclose(metrics.rms(2*np.ones((4, 5)), w), 2)
    assert metrics.rms(np.ones((7, 4, 5)), w).shape == (7,)
    o = np.ones((4, 5))
    assert np.allclose(metrics.gaussian_relative_entropy(o, o, o, o), 0)
    assert np.isclose(metrics.averaged_gre(o, o, o, o, w), 0)


def _staggered(ny=20, nx=24, h=2.0):
    y, x = np.mgrid[0:ny, 0:nx] * h
    return x, y, h


def test_okubo_weiss_analytic():
    x, y, h = _staggered()
    # Solid-body rotation: u at (x+h/2, y), v at (x, y+h/2)
    u, v = -y, x
    inner = (slice(2, -2), slice(2, -2))
    assert np.allclose(flow.vorticity(u, v, h, h)[inner], 2)
    assert np.allclose(flow.strain_normal(u, v, h, h)[inner], 0)
    assert np.allclose(flow.strain_shear(u, v, h, h)[inner], 0)
    assert np.allclose(flow.okubo_weiss(u, v, h, h)[inner], -4)
    # Pure normal strain
    assert np.allclose(flow.okubo_weiss(x, -y, h, h)[inner], 4)
    # Pure shear strain
    assert np.allclose(flow.okubo_weiss(y, x, h, h)[inner], 4)


def test_stats():
    u = np.ones((6, 8, 10))
    assert np.allclose(stats.kinetic_energy(u, u)[:, 1:-1, 1:-1], 1)
    q = np.arange(6.)[:, None, None]*np.ones((6, 8, 10))
    assert np.allclose(stats.fluctuations(q, 2)[:2].mean(axis=0), 0)
    assert np.allclose(stats.fluctuations(q, 2)[2:].mean(axis=0), 0)
    assert stats.hovmoller(q).shape == (6,)


def test_ke_spectrum_peak():
    n = 64
    x = np.arange(n)
    u = np.tile(np.sin(2*np.pi*5*x/n), (n, 1))
    k, spec = stats.ke_spectrum(u, 0*u, 1.0, 1.0)
    assert np.argmax(spec) == 4
    assert spectra.init_3Dwvnb(4, 4, 4, 1., 1., 1.)[1].size > 0


def test_gyrekit():
    assert len(runs.experiments("R3")) >= 2 and len(runs.experiments("R9")) >= 1
    f = next(iter(runs.experiments("R3").values()))["files"]["uvelfile"]
    assert f[0].endswith(".nc")
    lat, lon = np.meshgrid(np.linspace(15, 45, 40), np.linspace(-100, -70, 50), indexing="ij")
    m = grids.jet_mask(lat, lon, 100., 100.)
    assert m.shape == lat.shape and m.sum() > 0 and m[-1].sum() == 0
    p = np.ones((26, 26))
    assert grids.coarse_grain_P(p, 3).shape[0] > 0
    u = np.random.default_rng(0).normal(size=(30, 2, 16, 16))
    out = diag.velocity_summary(u, u, 1000., 1000., 86400., delay=5, nb=2)
    assert out["mke"].shape == (2, 16, 16) and out["hov_eke"].shape == (30, 2)
    assert out["psd_eke"].shape[0] == out["freq"].size


def test_netcdf_and_map():
    import tempfile
    import netCDF4
    with tempfile.TemporaryDirectory() as d:
        grid = os.path.join(d, "grid.nc")
        with netCDF4.Dataset(grid, "w") as nc:
            nc.createDimension("z", 2); nc.createDimension("y", 3); nc.createDimension("x", 4)
        out = os.path.join(d, "out.nc")
        netcdf.write_fields(grid, out, {"a": np.full((3, 4), np.nan),
                                        "b": np.ones((2, 3, 4))})
        with netCDF4.Dataset(out) as nc:
            assert nc["a"].dimensions == ("y", "x") and nc["b"].dimensions == ("z", "y", "x")
            assert nc["a"][:].mask.all() and np.allclose(nc["b"][:], 1)
        lat_c, lon_c = np.meshgrid(np.array([0., 10.]), np.array([0., 20.]), indexing="ij")
        lat_f, lon_f = np.meshgrid(np.arange(-5., 16.), np.arange(-5., 26.), indexing="ij")
        path = os.path.join(d, "map.nc")
        first = orca.load_or_make_map(path, lat_c, lon_c, lat_f, lon_f)
        again = orca.load_or_make_map(path, lat_c, lon_c, lat_f, lon_f)
        assert np.array_equal(first[0], again[0]) and np.array_equal(first[1], again[1])
        assert lat_f[first[0][1, 1], first[1][1, 1]] == 10 and lon_f[first[0][1, 1], first[1][1, 1]] == 20


def test_w_conserves_mass():
    from noisekit.pipelines.w_modes import snapshot_w
    ny, nx, nz = 12, 16, 4
    y, x = np.mgrid[0:ny, 0:nx].astype(float)
    # A constant du/dx = 1e-5 with uniform metrics, west of the wrap
    u = np.tile(1e-5*1000.*x, (nz, 1, 1))
    v = np.zeros_like(u)
    one = np.ones((ny, nx))
    m = {"e1t": 1000.*one,
         "e2t": 1000.*one,
         "e1v": 1000.*one,
         "e2u": 1000.*one,
         "tmask": np.ones((nz, ny, nx)),
         "e3t": 10.*np.ones((nz, ny, nx))}
    w = snapshot_w(u, v, m)
    # dw/dz = -div with k increasing downward: w at the surface is -div times the depth above the bottom cell
    assert np.allclose(w[0, 5, 5], -1e-5*10.*(nz-1)) and np.allclose(w[-1], 0)
    v2 = np.tile(1e-5*1000.*y, (nz, 1, 1))
    assert np.allclose(snapshot_w(0*v2, v2, m)[0, 5, 5], -1e-5*10.*(nz-1))


def test_correlation():
    rho = np.exp(-np.arange(1, 21)/5.0)[:, None]
    L, c = correlation.efolding_lag(rho, 20)
    assert abs(L[0]-5) < 0.02 and not c[0]
    assert correlation.efolding_lag(np.ones((5, 1)), 5)[1][0]
    rng = np.random.default_rng(3)
    # AR(1) in time with a known e-folding time
    a = 0.9
    x = np.zeros((20000, 1))
    for t in range(1, 20000):
        x[t] = a*x[t-1] + rng.normal()
    T = correlation.efolding_lag(correlation.temporal_correlation(x, 40), 40)[0][0]
    assert abs(T - (-1/np.log(a))) < 1.5, T
    from scipy.ndimage import gaussian_filter
    f = np.stack([gaussian_filter(rng.normal(size=(60, 200)), (0, 6), mode="wrap") for _ in range(100)])
    rx, ry = correlation.spatial_correlation([f, f], 20)
    Lx = np.nanmean(correlation.efolding_lag(rx[:, 20:40, 20:150], 20)[0])
    Ly = np.nanmean(correlation.efolding_lag(ry[:, 20:40, 20:150], 20)[0])
    assert Lx > 5*Ly, (Lx, Ly)


def test_structure():
    n = 60
    e = 1000.*np.ones((n, n))
    wet = np.ones((n, n), dtype=bool); wet[:10, :10] = False
    Dxx, Dyy, Dxy, k2 = structure.tensor_field(8000., 8000., 0., wet.shape)
    A, M = structure.assemble_A(Dxx, Dyy, Dxy, k2, e, e, wet)
    assert abs(A - A.T).max() < 1e-8
    rng = np.random.default_rng(4)
    b = rng.normal(size=n*n)*wet.ravel()
    import scipy.sparse.linalg as spla
    direct = spla.splu(A.tocsc()).solve(b)
    x, _ = structure.pcg_jacobi(A, b, 1/A.diagonal(), tol=1e-10)
    assert np.linalg.norm(x-direct)/np.linalg.norm(direct) < 1e-6
    lu = structure.build_filter(8000., 8000., 0., e, e, wet)
    s = structure.unit_variance_scale(lu, e, e, wet, rng, 300)
    v = np.mean([(structure.draw(lu, e, e, wet, rng)/s)[20:40, 20:40]**2 for _ in range(300)])
    assert abs(v-1) < 0.1, v
    f = structure.draw(lu, e, e, wet, rng)
    assert np.all(f[:10, :10] == 0)


def test_clean():
    from datakit import clean
    a = np.array([1.0, np.nan, np.inf, 11.0, -3.0])
    assert np.allclose(clean(a, fill=0.0), [1, 0, 0, 0, -3])
    assert np.isnan(clean(a)[3])


if __name__ == "__main__":
    fails = 0
    for k, f in list(globals().items()):
        if k.startswith("test_"):
            try:
                f(); print("PASS", k)
            except Exception as e:
                fails += 1; print("FAIL", k, type(e).__name__, str(e)[:300])
    sys.exit(fails)
