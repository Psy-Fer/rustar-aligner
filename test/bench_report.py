#!/usr/bin/env python3
"""Turn divan's output into a Markdown summary, optionally against a baseline.

Divan prints a tree; this pulls out the leaf rows (name, median) so a job
summary shows numbers rather than box drawing, and flags anything that moved
by more than 10% against the baseline run.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

# A leaf row: name, then four "value unit" columns separated by │. The third
# is the median, which is the one worth reporting.
ROW = re.compile(
    r"^[\s│├╰─]*([A-Za-z0-9_.]+)\s+"
    r"[\d.]+\s*(?:ns|µs|ms|s)\s*│\s*"
    r"[\d.]+\s*(?:ns|µs|ms|s)\s*│\s*"
    r"([\d.]+)\s*(ns|µs|ms|s)"
)
# A group header: a name with the columns empty.
GROUP = re.compile(r"^[\s│├╰─]*([A-Za-z0-9_]+)\s+│\s*│\s*│\s*│")
SCALE = {"ns": 1e-9, "µs": 1e-6, "ms": 1e-3, "s": 1.0}


def parse(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    if not path.exists():
        return out
    group = ""
    for line in path.read_text(errors="replace").splitlines():
        g = GROUP.match(line)
        if g:
            group = g.group(1)
            continue
        m = ROW.match(line)
        if not m:
            continue
        name, median, unit = m.groups()
        key = f"{group}/{name}" if group and group != name else name
        out[key] = float(median) * SCALE[unit]
    return out


def fmt(seconds: float) -> str:
    for unit, scale in (("s", 1.0), ("ms", 1e-3), ("µs", 1e-6), ("ns", 1e-9)):
        if seconds >= scale:
            return f"{seconds / scale:.2f} {unit}"
    return f"{seconds * 1e9:.2f} ns"


def main() -> int:
    head = parse(Path(sys.argv[1])) if len(sys.argv) > 1 else {}
    base = parse(Path(sys.argv[2])) if len(sys.argv) > 2 else {}

    if not head:
        print("### Benchmarks\n\nNo benchmark rows parsed; see the uploaded raw output.")
        return 0

    print("### Benchmarks (median)\n")
    if base:
        print("| benchmark | this ref | baseline | change |")
        print("|---|---|---|---|")
        for name in sorted(head):
            h = head[name]
            b = base.get(name)
            if b is None:
                print(f"| {name} | {fmt(h)} | — | new |")
                continue
            delta = (h - b) / b if b else 0.0
            flag = " ⚠️" if abs(delta) > 0.10 else ""
            print(f"| {name} | {fmt(h)} | {fmt(b)} | {delta:+.1%}{flag} |")
        print("\nA change above 10% is flagged; on a shared runner treat it as a")
        print("prompt to re-run rather than as a verdict.")
    else:
        print("| benchmark | median |")
        print("|---|---|")
        for name in sorted(head):
            print(f"| {name} | {fmt(head[name])} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
