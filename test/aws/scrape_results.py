#!/usr/bin/env python3
"""Scrape wall-clock + Max-RSS from the benchmark logs into a comparison table.

Usage:  python3 scrape_results.py <results_dir>

Reads the GNU `/usr/bin/time -v` blocks in:
  build_star.log build_rustar.log   (one-time index builds)
  rustar_<rep>.log star_<rep>.log cr_<rep>.log   (the align+count comparison)
Prints a Markdown table + writes results.json. rustar's Max-RSS includes
reclaimable mmap pages (see note), so it overstates real memory pressure.
"""
import json
import re
import statistics
import sys
from pathlib import Path

WALL = re.compile(r"Elapsed \(wall clock\) time.*?:\s*([0-9:.]+)")
RSS = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")


def to_seconds(s: str) -> float:
    parts = s.split(":")
    parts = [float(p) for p in parts]
    if len(parts) == 3:
        h, m, sec = parts
    elif len(parts) == 2:
        h, m, sec = 0, parts[0], parts[1]
    else:
        h, m, sec = 0, 0, parts[0]
    return h * 3600 + m * 60 + sec


def parse(path: Path):
    if not path.exists():
        return None
    text = path.read_text(errors="replace")
    w, r = WALL.search(text), RSS.search(text)
    if not w or not r:
        return None
    return {"wall_s": round(to_seconds(w.group(1)), 1),
            "max_rss_gb": round(int(r.group(1)) / 1048576, 2)}


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    d = Path(sys.argv[1])
    tools = {"rustar": "rustar", "STARsolo": "star", "CellRanger": "cr"}
    out = {"index_build": {}, "align_count": {}}

    for label, stem in {"STARsolo": "build_star", "rustar": "build_rustar"}.items():
        m = parse(d / f"{stem}.log")
        if m:
            out["index_build"][label] = m

    for label, stem in tools.items():
        reps = [parse(d / f"{stem}_{r}.log") for r in (1, 2, 3)]
        reps = [x for x in reps if x]
        if not reps:
            continue
        walls = [x["wall_s"] for x in reps]
        rsss = [x["max_rss_gb"] for x in reps]
        out["align_count"][label] = {
            "reps": reps,
            "wall_median_s": round(statistics.median(walls), 1),
            "wall_min_s": min(walls), "wall_max_s": max(walls),
            "max_rss_gb": max(rsss),
        }

    (d / "results.json").write_text(json.dumps(out, indent=2))

    print("\n### Index build (one-time)\n")
    print("| tool | wall (s) | Max RSS (GB) |")
    print("|---|--:|--:|")
    for label, m in out["index_build"].items():
        print(f"| {label} | {m['wall_s']:.0f} | {m['max_rss_gb']:.1f} |")

    print("\n### Align + count → matrix (native x86_64, fixed threads, no BAM)\n")
    print("| tool | wall median (s) | wall min–max | Max RSS (GB) |")
    print("|---|--:|--:|--:|")
    for label in tools:
        m = out["align_count"].get(label)
        if not m:
            print(f"| {label} | — (no log) | | |")
            continue
        print(f"| {label} | {m['wall_median_s']:.0f} | "
              f"{m['wall_min_s']:.0f}–{m['wall_max_s']:.0f} | {m['max_rss_gb']:.1f} |")

    print("\n_Note: rustar mmaps its index, so its Max RSS includes reclaimable "
          "file-backed pages and overstates real memory pressure vs STAR's "
          "anonymous read-into-RAM. CellRanger `count` is end-to-end "
          "(--create-bam=false --nosecondary)._")
    print(f"\nWrote {d/'results.json'}")


if __name__ == "__main__":
    main()
