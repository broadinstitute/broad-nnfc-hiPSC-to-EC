"""Rewrite the guide design inside the day objects, from a design folder.

Nothing else of a day object depends on the design, so only three elements are
rewritten, as `prepare_channels.py` and `assemble_days.py` built them:

- ``mod/guide/var``: `mudata_builder.guide_var` on the design's
  ``igvf/guide_metadata.tsv``, each column keeping its dtype;
- ``mod/guide/uns/annotation`` and ``mod/guide/uns/element_mapping``: the
  design's gene annotation and element mapping tables.

Each day is copied, patched, read back and checked against the design, and
only then moved over the original. Run it with the v2 folder to put v2 back.

Usage::

    pixi run --manifest-path ~/GitHub/crispr-analysis-utils/pyproject.toml -e dev \
        python patch_design.py --design results/2026-11-18_guide_design_v3 [--days 0 2 4]
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import h5py
import mudata as md
import pandas as pd
from anndata.io import read_elem, write_elem
from crispr_analysis_utils.mudata_builder import guide_var

REPO = Path(__file__).resolve().parents[2]
DAYS = REPO / "results" / "2026-11-12_day_objects" / "days"
GZIP = {"compression": "gzip"}  # as `prepare_channels.py` wrote the objects
TABLES = {
    "annotation": Path("annotations") / "guide_gene_annotation.tsv",
    "element_mapping": Path("elements") / "guide_element_mapping.tsv",
}


def as_text(frame: pd.DataFrame) -> pd.DataFrame:
    """Every value as a string and every gap as one empty string, to compare."""
    return frame.astype(object).where(frame.notna(), "").astype(str)


def expected_var(design: Path, old: pd.DataFrame) -> pd.DataFrame:
    """The design's guide var, in the object's guide order and dtypes."""
    new = guide_var(design / "igvf" / "guide_metadata.tsv", old.index)
    new.index = old.index
    if list(new.columns) != list(old.columns):
        raise SystemExit(f"columns differ: {list(old.columns)} vs {list(new.columns)}")
    for column in new.columns:
        if isinstance(old[column].dtype, pd.CategoricalDtype):
            new[column] = new[column].astype("category")
        else:
            new[column] = new[column].astype(old[column].dtype)
    return new


def patch(path: Path, design: Path) -> None:
    tables = {
        key: pd.read_csv(design / table, sep="\t", dtype=str)
        for key, table in TABLES.items()
    }
    partial = path.with_suffix(".patching.h5mu")
    shutil.copyfile(path, partial)
    with h5py.File(partial, "r+") as f:
        old = read_elem(f["mod/guide/var"])
        new = expected_var(design, old)
        del f["mod/guide/var"]
        write_elem(f["mod/guide"], "var", new, dataset_kwargs=GZIP)
        for key, table in tables.items():
            del f[f"mod/guide/uns/{key}"]
            write_elem(f["mod/guide/uns"], key, table, dataset_kwargs=GZIP)

    # Read back from the patched file, through mudata, before trusting it.
    data = md.read_h5mu(partial, backed="r")
    var = data["guide"].var
    # Values as text, and each column's kind as the object had it (an empty
    # column's categories read back as object whatever they were written as).
    kinds = [str(dtype) for dtype in var.dtypes]
    if not as_text(var).equals(as_text(new)) or kinds != [str(d) for d in old.dtypes]:
        raise SystemExit(f"{partial}: the guide var did not read back as written")
    for key, table in tables.items():
        if not as_text(data["guide"].uns[key]).equals(as_text(table)):
            raise SystemExit(f"{partial}: uns[{key!r}] did not read back as written")
    data.file.close()

    moved = (as_text(old) != as_text(new)).any(axis=1)
    columns = (as_text(old) != as_text(new)).sum()
    print(
        f"{path.name}: {int(moved.sum()):,} of {len(new):,} guides changed",
        dict(columns[columns > 0]),
    )
    partial.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--days", type=int, nargs="+", default=[0, 2, 4])
    args = parser.parse_args()
    design = args.design if args.design.is_absolute() else REPO / args.design
    for day in args.days:
        patch(DAYS / f"day{day}.h5mu", design)


if __name__ == "__main__":
    main()
