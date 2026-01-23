#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import os
import re
from pathlib import Path
from typing import Optional, Tuple

# Example filename:
# rep03_day00_chip_a_batch04_channel01_web_summary.html
PATTERN = re.compile(
    r"^rep(?P<rep>\d+)_day(?P<day>\d+)_chip_(?P<chip>[A-Za-z]+)_batch(?P<batch>\d+)_channel(?P<channel>\d+)_web_summary\.html$"
)

def parse_filename(name: str) -> Optional[Tuple[int, int, str, int, int]]:
    m = PATTERN.match(name)
    if not m:
        return None
    rep = int(m.group("rep"))
    day = int(m.group("day"))
    chip = m.group("chip")
    batch = int(m.group("batch"))
    channel = int(m.group("channel"))
    return rep, day, chip, batch, channel

def build_html(rows, title: str) -> str:
    # Minimal styling + sortable columns (click headers)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 24px; }}
    h1 {{ margin: 0 0 12px 0; font-size: 20px; }}
    .note {{ color: #555; margin-bottom: 16px; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #ddd; padding: 8px 10px; }}
    th {{ background: #f6f6f6; cursor: pointer; user-select: none; position: sticky; top: 0; }}
    tr:nth-child(even) {{ background: #fafafa; }}
    a {{ color: #0b5fff; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .small {{ font-size: 12px; color: #666; }}
  </style>
</head>
<body>
  <h1>{html.escape(title)}</h1>
  <div class="note small">Click a header to sort. Generated automatically.</div>

  <table id="tbl">
    <thead>
      <tr>
        <th data-type="num">Rep</th>
        <th data-type="num">Day</th>
        <th data-type="str">Chip</th>
        <th data-type="num">Batch</th>
        <th data-type="num">Channel</th>
        <th data-type="str">Link</th>
      </tr>
    </thead>
    <tbody>
      {rows}
    </tbody>
  </table>

<script>
(function() {{
  const table = document.getElementById("tbl");
  const getCellValue = (tr, idx) => tr.children[idx].innerText;

  const comparer = (idx, type, asc) => (a, b) => {{
    const v1 = getCellValue(asc ? a : b, idx).trim();
    const v2 = getCellValue(asc ? b : a, idx).trim();
    if (type === "num") {{
      const n1 = parseFloat(v1.replace(/[^0-9.\\-]/g, "")) || 0;
      const n2 = parseFloat(v2.replace(/[^0-9.\\-]/g, "")) || 0;
      return n1 - n2;
    }}
    return v1.localeCompare(v2);
  }};

  Array.from(table.querySelectorAll("th")).forEach((th, idx) => {{
    let asc = true;
    th.addEventListener("click", () => {{
      const type = th.getAttribute("data-type") || "str";
      const tbody = table.tBodies[0];
      const rows = Array.from(tbody.querySelectorAll("tr"));
      rows.sort(comparer(idx, type, asc));
      asc = !asc;
      rows.forEach(r => tbody.appendChild(r));
    }});
  }});
}})();
</script>

</body>
</html>
"""

def main():
    ap = argparse.ArgumentParser(
        description="Generate an HTML index table from Cell Ranger web summary filenames."
    )
    ap.add_argument("folder", help="Folder containing *_web_summary.html files")
    ap.add_argument(
        "--base-url",
        default="https://mitra.stanford.edu/engreitz/oak/public/emattei/hiPSC-EC/cellranger_reports/",
        help="URL prefix to prepend to each filename for the Link column",
    )
    ap.add_argument("--out", default="reports.html", help="Output HTML filename")
    ap.add_argument("--title", default="Cell Ranger Reports", help="HTML page title")
    args = ap.parse_args()

    folder = Path(args.folder)
    base_url = args.base_url.rstrip("/") + "/"

    entries = []
    for p in sorted(folder.iterdir()):
        if not p.is_file():
            continue
        parsed = parse_filename(p.name)
        if not parsed:
            continue
        rep, day, chip, batch, channel = parsed
        entries.append((rep, day, chip, batch, channel, p.name))

    # Sort naturally by Rep, Day, Chip, Batch, Channel
    entries.sort(key=lambda x: (x[0], x[1], x[2].lower(), x[3], x[4]))

    row_html = []
    for rep, day, chip, batch, channel, fname in entries:
        link = base_url + fname
        row_html.append(
            "<tr>"
            f"<td>{rep:02d}</td>"
            f"<td>{day:02d}</td>"
            f"<td>{html.escape(chip)}</td>"
            f"<td>{batch:02d}</td>"
            f"<td>{channel:02d}</td>"
            f"<td><a href=\"{html.escape(link)}\" target=\"_blank\" rel=\"noopener\">{html.escape(fname)}</a></td>"
            "</tr>"
        )

    out_html = build_html("\n      ".join(row_html), args.title)
    out_path = folder / args.out
    out_path.write_text(out_html, encoding="utf-8")
    print(f"Wrote: {out_path} ({len(entries)} rows)")

if __name__ == "__main__":
    main()

