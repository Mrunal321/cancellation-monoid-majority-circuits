#!/usr/bin/env python3
"""Run Vivado OOC evaluation for arbitrary Verilog netlists.

This is a lightweight wrapper around scripts/vivado_eval.tcl for cases where
we already have standalone Verilog netlists and want synth/post-route reports
without going through the full majority sweep flow.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import glob
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence


MODULE_RE = re.compile(r"^\s*module\s+([A-Za-z_][A-Za-z0-9_$]*)\b")

CSV_COLUMNS = [
    "label",
    "netlist",
    "top_module",
    "part",
    "flow_mode",
    "status",
    "lut_count",
    "ff_count",
    "carry_count",
    "wns_ns",
    "crit_delay_ns",
    "power_w",
    "util_rpt",
    "timing_rpt",
    "power_rpt",
    "vivado_log",
    "error",
]


def run_cmd(cmd: Sequence[str], capture: bool = False) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(
            list(cmd),
            check=False,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.PIPE if capture else None,
            text=True,
        )
    except FileNotFoundError as err:
        return subprocess.CompletedProcess(
            args=list(cmd),
            returncode=127,
            stdout="",
            stderr=str(err),
        )


def parse_int(text: str) -> Optional[int]:
    text = (text or "").strip().replace(",", "")
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def parse_float(text: str) -> Optional[float]:
    text = (text or "").strip().replace(",", "")
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def to_csv_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def parse_utilization_report(path: Path) -> Dict[str, Optional[int]]:
    text = path.read_text(encoding="utf-8", errors="ignore")

    def pick(patterns: Sequence[str]) -> Optional[int]:
        for pattern in patterns:
            match = re.search(pattern, text, flags=re.MULTILINE)
            if match:
                return parse_int(match.group(1))
        return None

    lut = pick(
        [
            r"^\|\s*CLB LUTs\*?\s*\|\s*([0-9,]+)\s*\|",
            r"^\|\s*Slice LUTs\*?\s*\|\s*([0-9,]+)\s*\|",
        ]
    )
    ff = pick(
        [
            r"^\|\s*CLB Registers\*?\s*\|\s*([0-9,]+)\s*\|",
            r"^\|\s*Slice Registers\*?\s*\|\s*([0-9,]+)\s*\|",
        ]
    )
    carry = pick(
        [
            r"^\|\s*CARRY4\*?\s*\|\s*([0-9,]+)\s*\|",
            r"^\|\s*CARRY8\*?\s*\|\s*([0-9,]+)\s*\|",
        ]
    )

    return {
        "lut_count": 0 if lut is None else lut,
        "ff_count": 0 if ff is None else ff,
        "carry_count": 0 if carry is None else carry,
    }


def parse_timing_report(path: Path, period_ns: float) -> Dict[str, Optional[float]]:
    text = path.read_text(encoding="utf-8", errors="ignore")

    wns: Optional[float] = None
    crit: Optional[float] = None

    m_slack = re.search(r"Slack\s*\([^\)]*\)\s*:\s*([-+]?\d+(?:\.\d+)?)ns", text)
    if m_slack:
        wns = parse_float(m_slack.group(1))

    m_data = re.search(r"Data Path Delay\s*:\s*([-+]?\d+(?:\.\d+)?)ns", text)
    if m_data:
        crit = parse_float(m_data.group(1))

    if crit is None and wns is not None:
        crit = period_ns - wns

    return {
        "wns_ns": wns,
        "crit_delay_ns": crit,
    }


def parse_power_report(path: Path) -> Dict[str, Optional[float]]:
    text = path.read_text(encoding="utf-8", errors="ignore")

    m = re.search(r"^\|\s*Total On-Chip Power \(W\)\s*\|\s*([-+]?\d+(?:\.\d+)?)\s*\|", text, flags=re.MULTILINE)
    if m:
        return {"power_w": parse_float(m.group(1))}

    m2 = re.search(r"Total On-Chip Power\s*\(W\)\s*[:|]\s*([-+]?\d+(?:\.\d+)?)", text)
    if m2:
        return {"power_w": parse_float(m2.group(1))}

    return {"power_w": None}


def parse_csv_list(text: str) -> List[str]:
    return [tok.strip() for tok in text.split(",") if tok.strip()]


def infer_top_module(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            match = MODULE_RE.match(line)
            if match:
                return match.group(1)
    raise RuntimeError(f"Could not infer top module from {path}")


def collect_netlists(explicit: List[str], patterns: List[str]) -> List[Path]:
    seen = set()
    out: List[Path] = []

    for item in explicit:
        path = Path(item).resolve()
        if not path.exists():
            raise RuntimeError(f"Netlist does not exist: {path}")
        if path not in seen:
            seen.add(path)
            out.append(path)

    for pattern in patterns:
        matches = sorted(Path(m).resolve() for m in glob.glob(pattern))
        if not matches:
            raise RuntimeError(f"No files matched --netlist-glob '{pattern}'")
        for path in matches:
            if path not in seen:
                seen.add(path)
                out.append(path)

    if not out:
        raise RuntimeError("Provide at least one --netlist or --netlist-glob")

    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Vivado synth/post-route evaluation for arbitrary Verilog netlists")
    parser.add_argument("--netlist", action="append", default=[], help="Path to a Verilog netlist; can be repeated")
    parser.add_argument("--netlist-glob", action="append", default=[], help="Glob pattern for Verilog netlists; can be repeated")
    parser.add_argument("--out-dir", default="results/vivado_netlist_eval", help="Parent output directory for this run")
    parser.add_argument("--flow-modes", default="post_route", help="Comma-separated flow modes: synth_only,post_route")
    parser.add_argument("--part", default="xc7a200tfbg484-1")
    parser.add_argument("--period-ns", type=float, default=5.0)
    parser.add_argument("--vivado-bin", default=os.environ.get("VIVADO_BIN", "vivado"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    flow_modes = parse_csv_list(args.flow_modes)
    if not flow_modes:
        raise RuntimeError("No flow modes selected")
    for flow_mode in flow_modes:
        if flow_mode not in {"synth_only", "post_route"}:
            raise RuntimeError(f"Unsupported flow mode '{flow_mode}'")

    vivado_probe = run_cmd([args.vivado_bin, "-version"], capture=True)
    if vivado_probe.returncode != 0:
        raise RuntimeError(
            f"Vivado not available via --vivado-bin '{args.vivado_bin}'. "
            f"stderr: {(vivado_probe.stderr or '').strip()}"
        )

    netlists = collect_netlists(args.netlist, args.netlist_glob)

    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = (repo_root / args.out_dir / f"run_{stamp}").resolve()
    reports_root = out_dir / "reports"
    logs_root = out_dir / "logs"
    out_dir.mkdir(parents=True, exist_ok=True)
    reports_root.mkdir(parents=True, exist_ok=True)
    logs_root.mkdir(parents=True, exist_ok=True)

    alias_root = Path(f"/tmp/vivado_netlist_eval_repo_alias_{os.getpid()}")
    if alias_root.exists() or alias_root.is_symlink():
        alias_root.unlink()
    alias_root.symlink_to(repo_root, target_is_directory=True)

    def alias_path(path: Path) -> Path:
        try:
            rel = path.resolve().relative_to(repo_root.resolve())
            return alias_root / rel
        except Exception:
            return path

    rows: List[Dict[str, str]] = []

    try:
        for netlist in netlists:
            top_module = infer_top_module(netlist)
            label = netlist.stem

            for flow_mode in flow_modes:
                report_dir = reports_root / label / flow_mode
                report_dir.mkdir(parents=True, exist_ok=True)

                util_rpt = report_dir / "utilization.rpt"
                timing_rpt = report_dir / "timing.rpt"
                power_rpt = report_dir / "power.rpt"
                vivado_log = logs_root / f"{label}_{flow_mode}.log"

                cmd = [
                    args.vivado_bin,
                    "-mode",
                    "batch",
                    "-source",
                    str(alias_path(repo_root / "scripts" / "vivado_eval.tcl")),
                    "-notrace",
                    "-tclargs",
                    str(alias_path(netlist)),
                    top_module,
                    args.part,
                    flow_mode,
                    f"{args.period_ns}",
                    str(alias_path(report_dir)),
                ]

                proc = run_cmd(cmd, capture=True)
                vivado_log.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")

                row = {
                    "label": label,
                    "netlist": str(netlist),
                    "top_module": top_module,
                    "part": args.part,
                    "flow_mode": flow_mode,
                    "status": "ok",
                    "lut_count": "",
                    "ff_count": "",
                    "carry_count": "",
                    "wns_ns": "",
                    "crit_delay_ns": "",
                    "power_w": "",
                    "util_rpt": str(util_rpt),
                    "timing_rpt": str(timing_rpt),
                    "power_rpt": str(power_rpt),
                    "vivado_log": str(vivado_log),
                    "error": "",
                }

                if proc.returncode != 0:
                    row["status"] = "error"
                    row["error"] = f"vivado failed rc={proc.returncode}"
                    rows.append(row)
                    continue

                missing = [p for p in (util_rpt, timing_rpt, power_rpt) if not p.exists()]
                if missing:
                    row["status"] = "error"
                    row["error"] = "missing reports: " + ", ".join(str(p) for p in missing)
                    rows.append(row)
                    continue

                util = parse_utilization_report(util_rpt)
                timing = parse_timing_report(timing_rpt, args.period_ns)
                power = parse_power_report(power_rpt)

                row["lut_count"] = to_csv_value(util["lut_count"])
                row["ff_count"] = to_csv_value(util["ff_count"])
                row["carry_count"] = to_csv_value(util["carry_count"])
                row["wns_ns"] = to_csv_value(timing["wns_ns"])
                row["crit_delay_ns"] = to_csv_value(timing["crit_delay_ns"])
                row["power_w"] = to_csv_value(power["power_w"])
                rows.append(row)

    finally:
        if alias_root.exists() or alias_root.is_symlink():
            alias_root.unlink()

    summary_csv = out_dir / "summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in CSV_COLUMNS})

    print(f"[INFO] Output directory: {out_dir}")
    print(f"[INFO] Summary CSV: {summary_csv}")
    print(f"[INFO] Rows: {len(rows)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as err:
        print(f"[ERROR] {err}", file=sys.stderr)
        raise SystemExit(1)
