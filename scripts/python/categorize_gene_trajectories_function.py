from future import annotations

import re from typing import Optional import numpy as np import pandas
as pd import matplotlib.pyplot as plt from
matplotlib.backends.backend_pdf import PdfPages

def _parse_day_from_col(col: str) -> Optional[int]: m =
re.search(r”day()“, str(col)) return int(m.group(1)) if m else None

def _anova_f_oneway(groups): clean_groups = [] for g in groups: g =
np.asarray(g, dtype=float) g = g[~np.isnan(g)] if g.size == 0: return
np.nan clean_groups.append(g)

    k = len(clean_groups)
    ns = np.array([g.size for g in clean_groups], dtype=float)
    N = ns.sum()
    if k < 2 or N <= k:
        return np.nan

    means = np.array([g.mean() for g in clean_groups], dtype=float)
    grand = np.sum(means * ns) / N

    ss_between = np.sum(ns * (means - grand) ** 2)
    ss_within = np.sum([np.sum((g - m) ** 2) for g, m in zip(clean_groups, means)])

    df_between = k - 1
    df_within = N - k
    if df_within <= 0:
        return np.nan

    ms_between = ss_between / df_between
    ms_within = ss_within / df_within
    if ms_within <= 0:
        return np.inf if ms_between > 0 else np.nan

    return float(ms_between / ms_within)

def _zscore_rows(x): mu = np.nanmean(x, axis=1, keepdims=True) sd =
np.nanstd(x, axis=1, keepdims=True) sd = np.where(sd == 0, 1.0, sd)
return (x - mu) / sd

def _corr_distance_matrix(X): C = np.corrcoef(X) C = np.clip(C, -1.0,
1.0) D = 1.0 - C np.fill_diagonal(D, 0.0) return D

def _cluster_pearson(X, n_clusters): from scipy.spatial.distance import
squareform from scipy.cluster.hierarchy import linkage, fcluster

    D = _corr_distance_matrix(X)
    dcond = squareform(D, checks=False)
    Z = linkage(dcond, method="average")
    labels = fcluster(Z, t=n_clusters, criterion="maxclust")
    return labels.astype(int)

def categorize_gene_trajectories( expr: pd.DataFrame, n_clusters: int,
f_threshold: float = 10.0, scale_for_clustering: bool = True,
scale_for_plots: bool = True, pdf_path: str = “trajectory_patterns.pdf”,
):

    col_days = {c: _parse_day_from_col(c) for c in expr.columns}
    days_sorted = sorted(set(col_days.values()))
    day_to_cols = {d: [c for c in expr.columns if col_days[c] == d] for d in days_sorted}

    genes = expr.index.astype(str).to_list()
    F = np.full(len(genes), np.nan, dtype=float)
    expr_values = expr.to_numpy(dtype=float)

    day_to_idx = {d: [expr.columns.get_loc(c) for c in cols] for d, cols in day_to_cols.items()}

    for i in range(expr_values.shape[0]):
        groups = [expr_values[i, day_to_idx[d]] for d in days_sorted]
        F[i] = _anova_f_oneway(groups)

    is_constant = np.isfinite(F) & (F <= float(f_threshold))

    traj = np.zeros((len(genes), len(days_sorted)), dtype=float)
    for j, d in enumerate(days_sorted):
        cols = day_to_cols[d]
        traj[:, j] = expr[cols].mean(axis=1).to_numpy(dtype=float)

    labels = np.full(len(genes), 0, dtype=int)
    idx_nc = np.where(~is_constant)[0]
    if idx_nc.size > 0:
        X = traj[idx_nc, :].copy()
        if scale_for_clustering:
            X = _zscore_rows(X)
        labs = _cluster_pearson(X, n_clusters=n_clusters)
        labels[idx_nc] = labs

    group = np.array(["constant"] * len(genes), dtype=object)
    for k in range(1, n_clusters + 1):
        group[(labels == k)] = f"C{k}"

    grouping_df = pd.DataFrame(
        {"gene": genes, "group": group, "is_constant": is_constant, "F_stat": F}
    )

    traj_plot = _zscore_rows(traj.copy()) if scale_for_plots else traj.copy()

    with PdfPages(pdf_path) as pdf:
        for gname in sorted(grouping_df["group"].unique()):
            idx = np.where(grouping_df["group"].to_numpy() == gname)[0]
            Xg = traj_plot[idx, :]
            mean_g = np.nanmean(Xg, axis=0)

            fig = plt.figure(figsize=(11, 8.5))
            ax = fig.add_subplot(111)

            for r in range(Xg.shape[0]):
                ax.plot(days_sorted, Xg[r, :], linewidth=1)

            ax.plot(days_sorted, mean_g, linewidth=3)
            ax.set_title(f"{gname} (n={idx.size})")
            ax.set_xticks(days_sorted)
            ax.set_ylim(min(ax.get_ylim()[0], 0), ax.get_ylim()[1])

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    return grouping_df, pdf_path
