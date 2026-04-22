from __future__ import annotations

import re
from typing import Iterable, Optional, Tuple

import numpy as np
import pandas as pd


def _parse_day_from_col(col: str) -> Optional[int]:
    """Extract an integer day from a column name.

    Expected pattern: contains the substring "day<INTEGER>" (e.g. "day0", "day12").
    Returns None when no day can be parsed.
    """

    match = re.search(r"day(\d+)", str(col))
    return int(match.group(1)) if match else None


def _anova_f_oneway(groups: Iterable[np.ndarray]) -> float:
    """Compute a one-way ANOVA F-statistic (no p-value).

    Notes:
    - This is a lightweight implementation to avoid depending on scipy.stats.
    - NaNs are dropped within each group.
    - Returns NaN if there isn't enough data to compute a statistic.
    """

    clean_groups: list[np.ndarray] = []
    for group in groups:
        arr = np.asarray(group, dtype=float)
        arr = arr[~np.isnan(arr)]
        if arr.size == 0:
            return float("nan")
        clean_groups.append(arr)

    k = len(clean_groups)
    group_sizes = np.array([g.size for g in clean_groups], dtype=float)
    n_total = float(group_sizes.sum())
    if k < 2 or n_total <= k:
        return float("nan")

    group_means = np.array([g.mean() for g in clean_groups], dtype=float)
    grand_mean = float(np.sum(group_means * group_sizes) / n_total)

    sum_squares_between = float(
        np.sum(group_sizes * (group_means - grand_mean) ** 2)
    )
    sum_squares_within = float(
        np.sum(
            [np.sum((group - mean) ** 2) for group, mean in zip(clean_groups, group_means)]
        )
    )

    degrees_of_freedom_between = k - 1
    degrees_of_freedom_within = int(n_total - k)
    if degrees_of_freedom_within <= 0:
        return float("nan")

    mean_square_between = sum_squares_between / degrees_of_freedom_between
    mean_square_within = sum_squares_within / degrees_of_freedom_within
    if mean_square_within <= 0:
        return float("inf") if mean_square_between > 0 else float("nan")

    return float(mean_square_between / mean_square_within)


def _zscore_rows(x: np.ndarray) -> np.ndarray:
    """Z-score each row independently, ignoring NaNs."""

    mu = np.nanmean(x, axis=1, keepdims=True)
    sd = np.nanstd(x, axis=1, keepdims=True)
    sd = np.where(sd == 0, 1.0, sd)
    return (x - mu) / sd


def _corr_distance_matrix(x: np.ndarray) -> np.ndarray:
    """Correlation distance matrix (1 - Pearson correlation).

    Returns an NxN matrix where N is the number of rows in x.
    """

    c = np.corrcoef(x)
    c = np.nan_to_num(c, nan=0.0)
    c = np.clip(c, -1.0, 1.0)
    d = 1.0 - c
    np.fill_diagonal(d, 0.0)
    return d


def _cluster_pearson(x: np.ndarray, n_clusters: int) -> np.ndarray:
    """Hierarchical clustering using correlation distance and average linkage."""

    try:
        from scipy.spatial.distance import squareform
        from scipy.cluster.hierarchy import fcluster, linkage
    except ImportError as e:  # pragma: no cover
        raise ImportError("scipy is required for clustering") from e

    d = _corr_distance_matrix(x)
    dcond = squareform(d, checks=False)
    z = linkage(dcond, method="average")
    labels = fcluster(z, t=int(n_clusters), criterion="maxclust")
    return labels.astype(int)


def categorize_gene_trajectories(
    expr: pd.DataFrame,
    n_clusters: int,
    f_threshold: float = 10.0,
    scale_for_clustering: bool = True,
    scale_for_plots: bool = True,
    pdf_path: str = "trajectory_patterns.pdf",
) -> Tuple[pd.DataFrame, str]:
    """Categorize gene expression trajectories across days.

    Overview:
    1) Parse days from column names (expects columns containing "day<INT>").
    2) For each gene (row), compute a 1-way ANOVA F-statistic across days.
       Genes with F <= f_threshold are treated as "constant" (not dynamic).
    3) For non-constant genes, compute a per-day mean trajectory and cluster
       trajectories using hierarchical clustering on correlation distance.
    4) Produce a PDF with one page per group ("constant", "C1", ...).

    Dependencies:
    - Clustering requires `scipy`.
    - PDF output requires `matplotlib`.

    Parameters
    - expr: rows are genes; columns are measurements/replicates. Column names must
      include day labels like "day0", "day3", etc.
    - n_clusters: number of trajectory clusters for non-constant genes.
    - f_threshold: ANOVA F-statistic cutoff for calling a gene "constant".
    - scale_for_clustering: if True, z-score each gene trajectory before clustering.
    - scale_for_plots: if True, z-score trajectories in the PDF for readability.
    - pdf_path: output PDF path.

    Returns
    - grouping_df: DataFrame with columns [gene, group, is_constant, F_stat].
    - pdf_path: the same path that was written.
    """

    # Map each column to a day integer (or None if unparseable).
    column_to_day = {col: _parse_day_from_col(str(col)) for col in expr.columns}
    unique_days_sorted = sorted({day for day in column_to_day.values() if day is not None})
    if not unique_days_sorted:
        raise ValueError("No day columns found; expected columns containing 'day<INT>'")

    # Group columns by day.
    columns_by_day = {
        day: [col for col in expr.columns if column_to_day[col] == day]
        for day in unique_days_sorted
    }

    gene_names = expr.index.astype(str).to_list()
    expression_values = expr.to_numpy(dtype=float)

    # Precompute index positions for each day group (faster than repeated loc).
    column_indices_by_day = {
        day: [expr.columns.get_loc(col) for col in cols]
        for day, cols in columns_by_day.items()
        if len(cols) > 0
    }

    # --- Step 1: ANOVA F-statistic per gene ---
    f_statistics = np.full(len(gene_names), np.nan, dtype=float)
    for gene_idx in range(expression_values.shape[0]):
        groups_by_day = [
            expression_values[gene_idx, column_indices_by_day[day]]
            for day in unique_days_sorted
        ]
        f_statistics[gene_idx] = _anova_f_oneway(groups_by_day)

    is_constant = np.isfinite(f_statistics) & (f_statistics <= float(f_threshold))

    # --- Step 2: build per-day mean trajectory per gene ---
    trajectories_by_day = np.zeros((len(gene_names), len(unique_days_sorted)), dtype=float)
    for day_idx, day in enumerate(unique_days_sorted):
        cols_for_day = columns_by_day[day]
        trajectories_by_day[:, day_idx] = expr[cols_for_day].mean(axis=1).to_numpy(
            dtype=float
        )

    # --- Step 3: cluster the non-constant genes ---
    cluster_labels = np.zeros(len(gene_names), dtype=int)
    non_constant_gene_indices = np.where(~is_constant)[0]
    if non_constant_gene_indices.size > 0:
        trajectories_for_clustering = trajectories_by_day[non_constant_gene_indices, :].copy()
        if scale_for_clustering:
            trajectories_for_clustering = _zscore_rows(trajectories_for_clustering)
        cluster_labels[non_constant_gene_indices] = _cluster_pearson(
            trajectories_for_clustering, n_clusters=int(n_clusters)
        )

    # --- Step 4: assemble output groups ---
    group_labels = np.array(["constant"] * len(gene_names), dtype=object)
    for k in range(1, int(n_clusters) + 1):
        group_labels[cluster_labels == k] = f"C{k}"

    grouping_df = pd.DataFrame(
        {
            "gene": gene_names,
            "group": group_labels,
            "is_constant": is_constant,
            "F_stat": f_statistics,
        }
    )

    # --- Step 5: PDF output ---
    # Matplotlib is only needed for PDF output; keep it as a runtime dependency
    # so importing this module doesn't fail in lightweight environments.
    try:
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "matplotlib is required to write the PDF output; install it or run in an environment with matplotlib"
        ) from e

    trajectories_for_plotting = (
        _zscore_rows(trajectories_by_day.copy())
        if scale_for_plots
        else trajectories_by_day.copy()
    )
    with PdfPages(pdf_path) as pdf:
        for gname in sorted(grouping_df["group"].unique()):
            gene_indices_in_group = np.where(grouping_df["group"].to_numpy() == gname)[0]
            if gene_indices_in_group.size == 0:
                continue

            group_trajectories = trajectories_for_plotting[gene_indices_in_group, :]
            mean_trajectory = np.nanmean(group_trajectories, axis=0)

            fig = plt.figure(figsize=(11, 8.5))
            ax = fig.add_subplot(111)

            # Individual gene trajectories.
            for r in range(group_trajectories.shape[0]):
                ax.plot(unique_days_sorted, group_trajectories[r, :], linewidth=1)

            # Cluster/group mean.
            ax.plot(unique_days_sorted, mean_trajectory, linewidth=3)
            ax.set_title(f"{gname} (n={gene_indices_in_group.size})")
            ax.set_xticks(unique_days_sorted)
            ax.set_ylim(min(ax.get_ylim()[0], 0), ax.get_ylim()[1])

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    return grouping_df, pdf_path
