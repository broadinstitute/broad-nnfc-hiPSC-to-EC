"""
Shift detection between two heatmaps H1, H2 with axes:
  x = genomic coordinate (columns or rows depending on your convention)
  k = kernel size (the other axis)

This file provides TWO methods:

A) Sliding-window normalized local cross-correlation (NCC) across x-shifts
   - best shift Δx(x0)
   - improvement over zero-shift: S(x0) = NCC_max - NCC_zero
   - optional permutation Z-score
   - plotting helpers

B) Optical-flow-like local displacement field (Δx(x0)) using OpenCV Farneback
   - treats heatmaps like images
   - estimates dense flow (dx, dy) between H1 and H2
   - aggregates to a 1D Δx(x0) track and a confidence track

Assumptions:
- H1 and H2 are numpy arrays with shape (K, X) by default:
    rows = kernel sizes (K), cols = genomic positions (X)
  If yours are (X, K), set axis_x=0 in method A and transpose before method B.

Dependencies:
- numpy, scipy, matplotlib
- Optional:
    - opencv-python for optical flow (method B)
"""

from __future__ import annotations
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
import matplotlib.pyplot as plt

# Optional (for permutation p-values / z-score convenience)
try:
    from scipy.stats import norm
except Exception:
    norm = None


# ---------------------------
# Utilities
# ---------------------------

def _zscore(a: np.ndarray, axis=None, eps: float = 1e-8) -> np.ndarray:
    m = np.mean(a, axis=axis, keepdims=True)
    s = np.std(a, axis=axis, keepdims=True)
    return (a - m) / (s + eps)

def _validate_same_shape(H1: np.ndarray, H2: np.ndarray):
    if H1.shape != H2.shape:
        raise ValueError(f"H1 and H2 must have the same shape. Got {H1.shape} vs {H2.shape}")


# ============================================================
# A) Sliding-window NCC across x shifts
# ============================================================

def local_shift_ncc(
    H1: np.ndarray,
    H2: np.ndarray,
    *,
    axis_x: int = 1,
    window_bp: int = 401,
    max_shift_bp: int = 200,
    step_bp: int = 1,
    normalize: str = "window",  # "window" or "per_k"
    return_full: bool = False,
):
    """
    Compute local shifts along genomic axis using normalized cross-correlation (NCC).

    Parameters
    ----------
    H1, H2 : np.ndarray
        Heatmaps. Default expected shape (K, X), where X is genomic coordinate.
    axis_x : int
        Axis corresponding to genomic coordinate. Default 1 (columns).
    window_bp : int
        Sliding window size in x (must be odd recommended).
    max_shift_bp : int
        Shifts tested in [-max_shift_bp, +max_shift_bp].
    step_bp : int
        Step for window centers along x.
    normalize : {"window", "per_k"}
        - "window": z-score whole patch (over both k and x in the window).
        - "per_k": z-score each kernel row within the window (over x only).
    return_full : bool
        If True, also returns the full correlation curves per x0 (can be large).

    Returns
    -------
    x_centers : np.ndarray, shape (N,)
        Window center indices along genomic axis.
    best_shift : np.ndarray, shape (N,)
        Best shift Δx(x0) in index units (bp if 1 index = 1 bp bin).
        Positive means H2 needs to be shifted RIGHT to match H1 (see note below).
    ncc_max : np.ndarray, shape (N,)
        Maximum NCC over shifts.
    ncc_zero : np.ndarray, shape (N,)
        NCC at shift 0.
    improvement : np.ndarray, shape (N,)
        ncc_max - ncc_zero (your S(x)).
    corr_curves (optional) : np.ndarray, shape (N, 2*max_shift_bp+1)
        NCC vs shift for each center.

    Note on sign:
        We compute corr(W1, W2 shifted by +dx). If dx > 0 maximizes correlation,
        then shifting H2 to the RIGHT aligns it to H1. That implies peaks in H2
        are LEFT of H1 by dx. Flip sign if you prefer the opposite convention.
    """
    _validate_same_shape(H1, H2)

    # Move x axis to last for easier slicing: (other_dims..., X)
    H1m = np.moveaxis(H1, axis_x, -1)
    H2m = np.moveaxis(H2, axis_x, -1)

    *other, X = H1m.shape
    if len(other) != 1:
        # We expect one "k" axis (K) and one "x" axis. If more dims, flatten them.
        H1m = H1m.reshape((-1, X))
        H2m = H2m.reshape((-1, X))
    K = H1m.shape[0]

    if window_bp > X:
        raise ValueError(f"window_bp ({window_bp}) cannot exceed X ({X}).")
    if window_bp % 2 == 0:
        # not required, but typical for symmetric windows
        window_bp += 1

    half = window_bp // 2
    shifts = np.arange(-max_shift_bp, max_shift_bp + 1, dtype=int)

    # valid centers where both window and shifted windows fit
    min_center = half + max_shift_bp
    max_center = X - half - 1 - max_shift_bp
    if max_center < min_center:
        raise ValueError(
            "Window+shift too large for X. "
            f"Need X >= window_bp + 2*max_shift_bp. Got X={X}, window_bp={window_bp}, max_shift_bp={max_shift_bp}."
        )

    x_centers = np.arange(min_center, max_center + 1, step_bp, dtype=int)
    N = len(x_centers)

    best_shift = np.zeros(N, dtype=int)
    ncc_max = np.full(N, np.nan, dtype=float)
    ncc_zero = np.full(N, np.nan, dtype=float)
    improvement = np.full(N, np.nan, dtype=float)

    corr_curves = np.full((N, shifts.size), np.nan, dtype=float) if return_full else None

    # Precompute patch indices relative to center
    rel = np.arange(-half, half + 1, dtype=int)  # length window_bp

    for i, c in enumerate(x_centers):
        # W1 patch
        idx1 = c + rel
        W1 = H1m[:, idx1]  # (K, window_bp)

        # Normalize W1 once
        if normalize == "window":
            W1n = _zscore(W1, axis=None)
        elif normalize == "per_k":
            W1n = _zscore(W1, axis=1)
        else:
            raise ValueError("normalize must be 'window' or 'per_k'")

        # Flatten for dot products
        v1 = W1n.ravel()
        denom1 = np.sqrt(np.dot(v1, v1)) + 1e-12

        # For each shift, compute NCC
        for j, dx in enumerate(shifts):
            idx2 = (c + dx) + rel
            W2 = H2m[:, idx2]

            if normalize == "window":
                W2n = _zscore(W2, axis=None)
            else:  # per_k
                W2n = _zscore(W2, axis=1)

            v2 = W2n.ravel()
            denom2 = np.sqrt(np.dot(v2, v2)) + 1e-12
            ncc = float(np.dot(v1, v2) / (denom1 * denom2))

            if return_full:
                corr_curves[i, j] = ncc

            if dx == 0:
                ncc_zero[i] = ncc

        if return_full:
            jmax = int(np.nanargmax(corr_curves[i]))
            ncc_max[i] = corr_curves[i, jmax]
            best_shift[i] = int(shifts[jmax])
        else:
            # recompute quickly without storing (still fine)
            # We'll just re-run to find max; simpler and still not too slow for moderate N.
            cvals = []
            for dx in shifts:
                idx2 = (c + dx) + rel
                W2 = H2m[:, idx2]
                W2n = _zscore(W2, axis=None) if normalize == "window" else _zscore(W2, axis=1)
                v2 = W2n.ravel()
                denom2 = np.sqrt(np.dot(v2, v2)) + 1e-12
                cvals.append(float(np.dot(v1, v2) / (denom1 * denom2)))
            cvals = np.array(cvals, dtype=float)
            jmax = int(np.argmax(cvals))
            ncc_max[i] = cvals[jmax]
            best_shift[i] = int(shifts[jmax])
            # ncc_zero already set above

        improvement[i] = ncc_max[i] - ncc_zero[i]

    if return_full:
        return x_centers, best_shift, ncc_max, ncc_zero, improvement, corr_curves
    return x_centers, best_shift, ncc_max, ncc_zero, improvement


def local_shift_permutation_zscore(
    H1: np.ndarray,
    H2: np.ndarray,
    *,
    axis_x: int = 1,
    window_bp: int = 401,
    max_shift_bp: int = 200,
    step_bp: int = 1,
    normalize: str = "window",
    n_perm: int = 200,
    perm_seed: int = 0,
):
    """
    Compute permutation-based Z-score for improvement S(x)=NCC_max-NCC_zero.

    Null: for each window center, we circularly shift H2 within the local window
    by a random offset (different per permutation), destroying alignment while
    preserving local marginal structure.

    Returns
    -------
    x_centers, best_shift, improvement, zscore, pvalue (if scipy available else None)
    """
    rng = np.random.default_rng(perm_seed)

    # First compute observed
    x_centers, best_shift, _, _, improvement = local_shift_ncc(
        H1, H2, axis_x=axis_x, window_bp=window_bp, max_shift_bp=max_shift_bp,
        step_bp=step_bp, normalize=normalize, return_full=False
    )

    # Prepare axis arrangement
    H1m = np.moveaxis(H1, axis_x, -1)
    H2m = np.moveaxis(H2, axis_x, -1)
    *other, X = H1m.shape
    if len(other) != 1:
        H1m = H1m.reshape((-1, X))
        H2m = H2m.reshape((-1, X))

    half = (window_bp if window_bp % 2 == 1 else window_bp + 1) // 2
    shifts_test = np.arange(-max_shift_bp, max_shift_bp + 1, dtype=int)
    rel = np.arange(-half, half + 1, dtype=int)

    # Compute null improvements
    null_impr = np.zeros((len(x_centers), n_perm), dtype=float)

    for p in range(n_perm):
        # Randomly phase-scramble H2 locally by circular shifts within each window
        # For each center c, we take the same window but circularly roll W2 by a random amount along x.
        for i, c in enumerate(x_centers):
            idx = c + rel
            W1 = H1m[:, idx]
            W2 = H2m[:, idx].copy()

            # roll each kernel row by a random amount within the window
            # (you can also roll the whole patch by one offset; per-row tends to be more stringent)
            rolls = rng.integers(0, W2.shape[1], size=W2.shape[0])
            for r in range(W2.shape[0]):
                W2[r] = np.roll(W2[r], int(rolls[r]))

            # Compute S_null = max NCC over shifts - NCC at 0, but with W2 "scrambled" baseline.
            # Here, since W2 is only defined on the window, we emulate shifts by rolling within window.
            # This yields a conservative null for local alignment.
            if normalize == "window":
                W1n = _zscore(W1, axis=None)
            else:
                W1n = _zscore(W1, axis=1)
            v1 = W1n.ravel()
            denom1 = np.sqrt(np.dot(v1, v1)) + 1e-12

            cvals = []
            for dx in shifts_test:
                W2s = np.roll(W2, shift=dx, axis=1)
                W2n = _zscore(W2s, axis=None) if normalize == "window" else _zscore(W2s, axis=1)
                v2 = W2n.ravel()
                denom2 = np.sqrt(np.dot(v2, v2)) + 1e-12
                cvals.append(float(np.dot(v1, v2) / (denom1 * denom2)))
            cvals = np.array(cvals, dtype=float)
            Snull = float(np.max(cvals) - cvals[dx_to_index(0, max_shift_bp)])
            null_impr[i, p] = Snull

    mu = null_impr.mean(axis=1)
    sd = null_impr.std(axis=1) + 1e-8
    z = (improvement - mu) / sd

    pval = None
    if norm is not None:
        # one-sided: large z indicates stronger-than-null shift evidence
        pval = 1.0 - norm.cdf(z)

    return x_centers, best_shift, improvement, z, pval


def dx_to_index(dx: int, max_shift: int) -> int:
    """Map dx in [-max_shift, max_shift] to an index in the shift array."""
    return int(dx + max_shift)


def plot_shift_tracks(x_centers, best_shift, improvement, z=None, title_prefix=""):
    """
    Plot:
      - best shift Δx(x)
      - improvement S(x) = NCC_max - NCC_0
      - optional Z-score
    """
    fig = plt.figure(figsize=(12, 7))
    ax1 = plt.subplot(3 if z is not None else 2, 1, 1)
    ax1.plot(x_centers, best_shift)
    ax1.set_ylabel("best shift Δx (bins)")
    ax1.set_title(f"{title_prefix}Best local shift vs x")

    ax2 = plt.subplot(3 if z is not None else 2, 1, 2, sharex=ax1)
    ax2.plot(x_centers, improvement)
    ax2.set_ylabel("S(x)=NCC_max - NCC_0")
    ax2.set_title(f"{title_prefix}Shift-driven similarity improvement")

    if z is not None:
        ax3 = plt.subplot(3, 1, 3, sharex=ax1)
        ax3.plot(x_centers, z)
        ax3.set_ylabel("Z-score")
        ax3.set_xlabel("genomic position (index)")
        ax3.set_title(f"{title_prefix}Permutation Z-score for S(x)")
    else:
        ax2.set_xlabel("genomic position (index)")

    plt.tight_layout()
    return fig


# ============================================================
# B) Optical-flow-like displacement field (image-based)
# ============================================================

def optical_flow_displacement(
    H1: np.ndarray,
    H2: np.ndarray,
    *,
    axis_x: int = 1,
    pyr_scale: float = 0.5,
    levels: int = 3,
    winsize: int = 25,
    iterations: int = 3,
    poly_n: int = 5,
    poly_sigma: float = 1.2,
    aggregate: str = "weighted_mean",  # how to reduce 2D flow to 1D track over x
    weight_power: float = 1.0,         # used if weighted_mean
):
    """
    Estimate a dense displacement field between heatmaps using Farneback optical flow.
    Returns a 1D Δx(x) track and a confidence track.

    Requires: opencv-python

    Parameters
    ----------
    H1, H2 : np.ndarray
        Heatmaps. Default expected shape (K, X).
    axis_x : int
        Genomic coordinate axis.
    aggregate : {"mean", "median", "weighted_mean"}
        How to reduce flow dx(K,X) -> dx_track(X)
    weight_power : float
        weights = (|grad| + eps)^weight_power for weighted aggregation.

    Returns
    -------
    dx_track : np.ndarray, shape (X,)
        Estimated displacement along x at each genomic position.
        Positive dx means pixels in H2 are displaced to the right relative to H1
        (OpenCV convention: flow gives where a pixel in H1 moves to in H2).
    conf_track : np.ndarray, shape (X,)
        A simple confidence proxy (higher = more texture / signal).
    flow : np.ndarray, shape (K, X, 2)
        Dense flow field (dx, dy).
    """
    _validate_same_shape(H1, H2)

    try:
        import cv2
    except ImportError as e:
        raise ImportError("Method B requires opencv-python. Install with: pip install opencv-python") from e

    # Move to (K, X) with x as columns
    if axis_x != 1:
        H1m = np.moveaxis(H1, axis_x, 1)
        H2m = np.moveaxis(H2, axis_x, 1)
    else:
        H1m, H2m = H1, H2

    # Convert to float32 and normalize to 0..1 (improves flow stability)
    A = H1m.astype(np.float32)
    B = H2m.astype(np.float32)

    # Contrast normalize
    A = _zscore(A, axis=None).astype(np.float32)
    B = _zscore(B, axis=None).astype(np.float32)

    # OpenCV expects single-channel 2D images; ours are 2D (K x X), OK.
    flow = cv2.calcOpticalFlowFarneback(
        prev=A,
        next=B,
        flow=None,
        pyr_scale=pyr_scale,
        levels=levels,
        winsize=winsize,
        iterations=iterations,
        poly_n=poly_n,
        poly_sigma=poly_sigma,
        flags=0,
    )
    # flow[...,0]=dx, flow[...,1]=dy in pixel units
    dx = flow[..., 0]  # (K, X)
    dy = flow[..., 1]  # (K, X)

    # Confidence proxy: gradient magnitude of A (places with structure provide better flow)
    # (You can swap to |A| or local variance depending on your data.)
    gy, gx = np.gradient(A)
    gradmag = np.sqrt(gx * gx + gy * gy) + 1e-6

    if aggregate == "mean":
        dx_track = dx.mean(axis=0)
        conf_track = gradmag.mean(axis=0)
    elif aggregate == "median":
        dx_track = np.median(dx, axis=0)
        conf_track = np.median(gradmag, axis=0)
    elif aggregate == "weighted_mean":
        w = gradmag ** float(weight_power)
        dx_track = (dx * w).sum(axis=0) / (w.sum(axis=0) + 1e-12)
        conf_track = w.mean(axis=0)
    else:
        raise ValueError("aggregate must be one of: mean, median, weighted_mean")

    return dx_track, conf_track, flow


def plot_optical_flow_tracks(dx_track, conf_track=None, title_prefix=""):
    fig = plt.figure(figsize=(12, 5))
    ax1 = plt.subplot(2 if conf_track is not None else 1, 1, 1)
    ax1.plot(dx_track)
    ax1.set_ylabel("Δx (pixels/bins)")
    ax1.set_title(f"{title_prefix}Optical-flow Δx(x)")

    if conf_track is not None:
        ax2 = plt.subplot(2, 1, 2, sharex=ax1)
        ax2.plot(conf_track)
        ax2.set_ylabel("confidence proxy")
        ax2.set_xlabel("genomic position (index)")
        ax2.set_title(f"{title_prefix}Optical-flow confidence")
    else:
        ax1.set_xlabel("genomic position (index)")

    plt.tight_layout()
    return fig


# ============================================================
# Example usage (edit to your arrays)
# ============================================================
if __name__ == "__main__":
    # Example: H1, H2 are numpy arrays of shape (K, X)
    # Replace these with your real arrays.
    K, X = 40, 5000
    rng = np.random.default_rng(0)
    H1 = rng.normal(size=(K, X)) * 0.2
    H2 = H1.copy()

    # Inject a shifted peak region for demo
    for k in range(K):
        H1[k, 2000:2050] += np.exp(-((np.arange(50) - 25) ** 2) / (2 * (5 + k/10) ** 2))
        H2[k, 2020:2070] += np.exp(-((np.arange(50) - 25) ** 2) / (2 * (5 + k/10) ** 2))

    # ---- A) NCC sliding window
    x_centers, best_dx, ncc_max, ncc_0, S = local_shift_ncc(
        H1, H2,
        axis_x=1,
        window_bp=401,
        max_shift_bp=80,
        step_bp=5,
        normalize="per_k",
    )
    plot_shift_tracks(x_centers, best_dx, S, title_prefix="NCC ")

    # Optional Z-score (can be slow; start small n_perm)
    # x_centers, best_dx, S, z, p = local_shift_permutation_zscore(
    #     H1, H2,
    #     axis_x=1, window_bp=401, max_shift_bp=80, step_bp=10,
    #     normalize="per_k", n_perm=100, perm_seed=0
    # )
    # plot_shift_tracks(x_centers, best_dx, S, z=z, title_prefix="NCC+Perm ")

    # ---- B) Optical flow
    try:
        dx_track, conf_track, flow = optical_flow_displacement(
            H1, H2, axis_x=1, winsize=25, iterations=3, aggregate="weighted_mean", weight_power=1.0
        )
        plot_optical_flow_tracks(dx_track, conf_track, title_prefix="OptFlow ")
    except ImportError:
        print("opencv-python not installed; skipping optical flow demo.")

    plt.show()
