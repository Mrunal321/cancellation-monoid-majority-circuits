#!/usr/bin/env python3
"""Merge run_cancel_tree_eval outputs with ABC optimized/mapped metrics.

For each provided run directory (which must contain master_results.csv), this script:
1) Picks one preferred row per (logic_repr, method, n) using flow priority
   post_route > synth_only > ntk_only.
2) Runs ABC on the row's netlist file to collect:
   - original AIG stats after strash
   - optimized AIG stats after dc2
   - mapped LUT-style stats after if -K 6
3) Writes one merged CSV with all original row fields + ABC fields.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import shlex
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
AND_RE = re.compile(r"\band\s*=\s*(\d+)\b.*\blev\s*=\s*(\d+)\b")
MAP_RE = re.compile(r"\bnd\s*=\s*(\d+)\b.*\bedge\s*=\s*(\d+)\b.*\baig\s*=\s*(\d+)\b.*\blev\s*=\s*(\d+)\b")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Merge evaluation CSVs with ABC optimized/mapped metrics")
    ap.add_argument("--run-dir", action="append", required=True, help="Run directory containing master_results.csv")
    ap.add_argument("--out-csv", required=True, help="Path to output merged CSV")
    ap.add_argument("--abc-bin", default="abc", help="ABC binary path (default: abc)")
    ap.add_argument("--k-lut", type=int, default=6, help="K for ABC if mapping (default: 6)")
    return ap.parse_args()


def strip_ansi(text: str) -> str:
    return ANSI_RE.sub("", text)


def parse_abc_output(text: str) -> Dict[str, str]:
    lines = [strip_ansi(line).strip() for line in text.splitlines() if "i/o" in line and "lev" in line]
    and_lines: List[Tuple[int, int]] = []
    map_line: Optional[Tuple[int, int, int, int]] = None

    for line in lines:
        mm = MAP_RE.search(line)
        if mm:
            map_line = (int(mm.group(1)), int(mm.group(2)), int(mm.group(3)), int(mm.group(4)))
            continue
        ma = AND_RE.search(line)
        if ma:
            and_lines.append((int(ma.group(1)), int(ma.group(2))))

    out: Dict[str, str] = {
        "abc_orig_and": "",
        "abc_orig_lev": "",
        "abc_opt_and": "",
        "abc_opt_lev": "",
        "abc_map_nd": "",
        "abc_map_edge": "",
        "abc_map_aig": "",
        "abc_map_lev": "",
    }

    if len(and_lines) >= 1:
        out["abc_orig_and"] = str(and_lines[0][0])
        out["abc_orig_lev"] = str(and_lines[0][1])
    if len(and_lines) >= 2:
        out["abc_opt_and"] = str(and_lines[1][0])
        out["abc_opt_lev"] = str(and_lines[1][1])
    if map_line is not None:
        out["abc_map_nd"] = str(map_line[0])
        out["abc_map_edge"] = str(map_line[1])
        out["abc_map_aig"] = str(map_line[2])
        out["abc_map_lev"] = str(map_line[3])

    return out


def run_abc_metrics(abc_bin: str, verilog_path: pathlib.Path, k_lut: int) -> Dict[str, str]:
    cmd_script = (
        f"read_verilog {shlex.quote(str(verilog_path))}; "
        "strash; print_stats; "
        "dc2; print_stats; "
        f"if -K {k_lut}; print_stats"
    )
    proc = subprocess.run([abc_bin, "-c", cmd_script], capture_output=True, text=True, check=False)

    out = {
        "abc_status": "ok" if proc.returncode == 0 else "error",
        "abc_error": "",
    }
    if proc.returncode != 0:
        err = (proc.stderr or "").strip()
        if not err:
            err = (proc.stdout or "").strip().splitlines()[-1] if (proc.stdout or "").strip() else f"abc rc={proc.returncode}"
        out["abc_error"] = err
        out.update(
            {
                "abc_orig_and": "",
                "abc_orig_lev": "",
                "abc_opt_and": "",
                "abc_opt_lev": "",
                "abc_map_nd": "",
                "abc_map_edge": "",
                "abc_map_aig": "",
                "abc_map_lev": "",
            }
        )
        return out

    parsed = parse_abc_output((proc.stdout or "") + "\n" + (proc.stderr or ""))
    out.update(parsed)
    return out


def priority(flow_mode: str) -> int:
    return {"post_route": 3, "synth_only": 2, "ntk_only": 1}.get(flow_mode, 0)


def load_preferred_rows(master_csv: pathlib.Path) -> List[Dict[str, str]]:
    by_key: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    with master_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("status", "") != "ok":
                continue
            key = (row.get("logic_repr", ""), row.get("method", ""), row.get("n", ""))
            if key not in by_key or priority(row.get("flow_mode", "")) > priority(by_key[key].get("flow_mode", "")):
                by_key[key] = row
    return list(by_key.values())


def main() -> int:
    args = parse_args()
    run_dirs = [pathlib.Path(d).resolve() for d in args.run_dir]

    merged_rows: List[Dict[str, str]] = []
    for run_dir in run_dirs:
        master_csv = run_dir / "master_results.csv"
        if not master_csv.exists():
            print(f"[WARN] Missing master_results.csv: {master_csv}", file=sys.stderr)
            continue

        rows = load_preferred_rows(master_csv)
        for row in rows:
            row_out = dict(row)
            row_out["source_run_dir"] = str(run_dir)
            verilog_file = pathlib.Path(row.get("file", ""))
            if not verilog_file.exists():
                row_out.update(
                    {
                        "abc_status": "error",
                        "abc_error": f"missing netlist: {verilog_file}",
                        "abc_orig_and": "",
                        "abc_orig_lev": "",
                        "abc_opt_and": "",
                        "abc_opt_lev": "",
                        "abc_map_nd": "",
                        "abc_map_edge": "",
                        "abc_map_aig": "",
                        "abc_map_lev": "",
                    }
                )
            else:
                row_out.update(run_abc_metrics(args.abc_bin, verilog_file, args.k_lut))
            merged_rows.append(row_out)

    merged_rows.sort(key=lambda r: (r.get("logic_repr", ""), r.get("method", ""), int(r.get("n", "0"))))

    out_csv = pathlib.Path(args.out_csv).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames: List[str] = []
    for row in merged_rows:
        for k in row.keys():
            if k not in fieldnames:
                fieldnames.append(k)

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in merged_rows:
            writer.writerow(row)

    print(f"[INFO] Wrote merged CSV: {out_csv}")
    print(f"[INFO] Rows: {len(merged_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
