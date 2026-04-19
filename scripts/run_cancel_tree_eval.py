#!/usr/bin/env python3
"""End-to-end cancel_tree evaluation flow.

Steps:
1) Build tools.
2) Generate majority netlists for odd n range and selected methods.
3) Collect Mockturtle-level stats into master CSV.
4) Optionally run Vivado (synth_only/post_route) and enrich CSV.
5) Generate plots, LaTeX/CSV tables, and a short markdown report (+PDF if pandoc is available).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import numpy as np
except Exception:  # pragma: no cover - runtime dependency check
    np = None

CSV_COLUMNS = [
    "method",
    "logic_repr",
    "n",
    "module",
    "file",
    "part",
    "flow_mode",
    "fa_count",
    "mig_maj_count",
    "aig_and_count",
    "inv_count",
    "popcount_levels",
    "strict_schedule_mode",
    "strict_scaffold_p",
    "strict_scaffold_inputs",
    "strict_scaffold_threshold",
    "strict_comparator_width",
    "strict_num_fixed_pairs",
    "strict_csa_fa_count",
    "strict_comparator_fa_count",
    "strict_total_fa_count",
    "strict_csa_levels",
    "strict_total_levels",
    "status",
    "cell_count",
    "lut_count",
    "ff_count",
    "carry_count",
    "wns_ns",
    "crit_delay_ns",
    "power_w",
    "util_rpt",
    "timing_rpt",
    "power_rpt",
    "error",
    "cancel_merge_depth",
    "cancel_max_k_width",
    "cancel_num_merges",
    "xag_xor_count",
    "xag_and_count",
    "ntk_depth_logic",
    "approx_group_size",
    "approx_compressed_n",
    "approx_groups3",
    "approx_passthrough_bits",
]

ACCURACY_COLUMNS = [
    "method",
    "n",
    "p_input",
    "trials",
    "accuracy",
    "error_rate",
]


@dataclass
class GeneratedDesign:
    method: str
    logic_repr: str
    n: int
    module: str
    verilog_file: Path
    stats_json: Path


def run_cmd(cmd: Sequence[str], cwd: Optional[Path] = None, capture: bool = False) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(
            list(cmd),
            cwd=str(cwd) if cwd else None,
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


def odd_values(n_min: int, n_max: int) -> List[int]:
    n0 = n_min if n_min % 2 == 1 else n_min + 1
    return list(range(n0, n_max + 1, 2))


def parse_probability_list(text: str) -> List[float]:
    ps: List[float] = []
    for tok in text.split(","):
        t = tok.strip()
        if not t:
            continue
        v = float(t)
        if v <= 0.0 or v >= 1.0:
            raise ValueError(f"Probability must be in (0,1): {v}")
        ps.append(v)
    if not ps:
        raise ValueError("No probabilities provided")
    return ps


def parse_csv_list(text: str) -> List[str]:
    return [tok.strip() for tok in text.split(",") if tok.strip()]


def row_logic_repr(row: Dict[str, str]) -> str:
    return row.get("logic_repr", "").strip()


def variant_label(method: str, logic_repr: str) -> str:
    lr = logic_repr.strip()
    if not lr:
        return method
    return f"{method}[{lr}]"


def row_variant_label(row: Dict[str, str]) -> str:
    return variant_label(row.get("method", ""), row_logic_repr(row))


def is_exact_method(method: str) -> bool:
    return method in {
        "popcount",
        "popcount_strict",
        "baseline_strict",
        "cancel_tree",
        "boyermoore_tree",
        "cancel_tree_v2",
        "boyermoore_tree_v2",
    }


def is_approx_block3_method(method: str) -> bool:
    return method in {"approx_block3", "approx_block3_popcount"}


def stable_method_seed(method: str) -> int:
    h = hashlib.sha256(method.encode("utf-8")).digest()
    return int.from_bytes(h[:8], byteorder="little", signed=False)


def estimate_accuracy_for_method(method: str, n: int, p: float, trials: int, seed: int, batch_size: int) -> float:
    if is_exact_method(method):
        return 1.0

    if not is_approx_block3_method(method):
        raise RuntimeError(
            f"Accuracy model not implemented for method '{method}'. "
            "Use only exact methods and/or approx_block3 methods for --run-accuracy."
        )

    if np is None:
        raise RuntimeError("numpy is required for --run-accuracy but is not installed")

    # Deterministic per (method, n, p) so accuracy sweeps are reproducible.
    local_seed = (seed + stable_method_seed(method) + 1000003 * n + int(round(p * 100000.0))) & ((1 << 63) - 1)
    rng = np.random.default_rng(local_seed)

    hits = 0
    done = 0

    groups = n // 3
    rem = n % 3
    maj_threshold = n // 2

    while done < trials:
        cur = min(batch_size, trials - done)
        x = (rng.random((cur, n)) < p).astype(np.uint8)

        y_exact = x.sum(axis=1) > maj_threshold

        if groups > 0:
            main = x[:, : 3 * groups].reshape(cur, groups, 3)
            compressed = (main.sum(axis=2) >= 2).astype(np.uint8)
        else:
            compressed = np.zeros((cur, 0), dtype=np.uint8)

        if rem > 0:
            compressed = np.concatenate([compressed, x[:, 3 * groups :]], axis=1)

        y_approx = compressed.sum(axis=1) > (compressed.shape[1] // 2)
        hits += int(np.sum(y_exact == y_approx))
        done += cur

    return float(hits) / float(trials)


def run_accuracy_sweep(methods: List[str], ns: List[int], ps: List[float], trials: int, seed: int, batch_size: int) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for method in methods:
        for n in ns:
            for p in ps:
                acc = estimate_accuracy_for_method(method=method, n=n, p=p, trials=trials, seed=seed, batch_size=batch_size)
                rows.append(
                    {
                        "method": method,
                        "n": str(n),
                        "p_input": f"{p:.6f}",
                        "trials": str(trials),
                        "accuracy": f"{acc:.6f}",
                        "error_rate": f"{(1.0 - acc):.6f}",
                    }
                )
    return rows


def write_accuracy_csv(rows: List[Dict[str, str]], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=ACCURACY_COLUMNS)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in ACCURACY_COLUMNS})


def save_accuracy_vs_n_plots(accuracy_rows: List[Dict[str, str]], out_dir: Path) -> None:
    if not accuracy_rows:
        return

    grouped: Dict[Tuple[float, str], List[Tuple[int, float]]] = {}
    p_values = sorted({float(r["p_input"]) for r in accuracy_rows})
    methods = sorted({r["method"] for r in accuracy_rows})

    for r in accuracy_rows:
        p = float(r["p_input"])
        method = r["method"]
        grouped.setdefault((p, method), []).append((int(r["n"]), float(r["accuracy"])))

    for p in p_values:
        plt.figure(figsize=(9, 5))
        drawn = False
        for method in methods:
            pts = sorted(grouped.get((p, method), []), key=lambda t: t[0])
            if not pts:
                continue
            xs = [x for x, _ in pts]
            ys = [y for _, y in pts]
            plt.plot(xs, ys, marker="o", linewidth=1.5, markersize=3, label=method)
            drawn = True
        if not drawn:
            plt.close()
            continue
        plt.xlabel("n")
        plt.ylabel("Accuracy")
        plt.title(f"Accuracy vs n (p={p:.2f})")
        plt.ylim(0.0, 1.01)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        out_dir.mkdir(parents=True, exist_ok=True)
        tag = f"{int(round(p * 100)):02d}"
        plt.savefig(str((out_dir / f"accuracy_vs_n_p{tag}").with_suffix(".png")), dpi=180)
        plt.savefig(str((out_dir / f"accuracy_vs_n_p{tag}").with_suffix(".pdf")))
        plt.close()


def save_accuracy_vs_p_plot(accuracy_rows: List[Dict[str, str]], out_prefix: Path) -> None:
    if not accuracy_rows:
        return

    grouped: Dict[Tuple[str, float], List[float]] = {}
    methods = sorted({r["method"] for r in accuracy_rows})
    p_values = sorted({float(r["p_input"]) for r in accuracy_rows})

    for r in accuracy_rows:
        method = r["method"]
        p = float(r["p_input"])
        grouped.setdefault((method, p), []).append(float(r["accuracy"]))

    plt.figure(figsize=(9, 5))
    drawn = False
    for method in methods:
        xs: List[float] = []
        ys: List[float] = []
        for p in p_values:
            vals = grouped.get((method, p), [])
            if not vals:
                continue
            xs.append(p)
            ys.append(statistics.mean(vals))
        if xs:
            plt.plot(xs, ys, marker="o", linewidth=1.5, markersize=4, label=method)
            drawn = True
    if not drawn:
        plt.close()
        return
    plt.xlabel("Input Bernoulli p")
    plt.ylabel("Mean Accuracy Across n")
    plt.title("Accuracy vs Input Bias")
    plt.ylim(0.0, 1.01)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_prefix.with_suffix(".png")), dpi=180)
    plt.savefig(str(out_prefix.with_suffix(".pdf")))
    plt.close()


def build_accuracy_tables(accuracy_rows: List[Dict[str, str]], out_csv: Path, out_tex: Path) -> List[Dict[str, str]]:
    if not accuracy_rows:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        with out_csv.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["method", "p_input", "mean_accuracy", "min_accuracy", "max_accuracy"])
            writer.writeheader()
        out_tex.parent.mkdir(parents=True, exist_ok=True)
        with out_tex.open("w", encoding="utf-8") as f:
            f.write("\\begin{tabular}{llrrr}\n\\hline\nMethod & p & Mean Acc & Min Acc & Max Acc \\\\ \n\\hline\n\\hline\n\\end{tabular}\n")
        return []

    methods = sorted({r["method"] for r in accuracy_rows})
    ps = sorted({float(r["p_input"]) for r in accuracy_rows})

    rows = []
    for method in methods:
        for p in ps:
            vals = [float(r["accuracy"]) for r in accuracy_rows if r["method"] == method and float(r["p_input"]) == p]
            if not vals:
                continue
            rows.append(
                {
                    "method": method,
                    "p_input": f"{p:.6f}",
                    "mean_accuracy": statistics.mean(vals),
                    "min_accuracy": min(vals),
                    "max_accuracy": max(vals),
                }
            )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["method", "p_input", "mean_accuracy", "min_accuracy", "max_accuracy"])
        writer.writeheader()
        for r in rows:
            writer.writerow({k: to_csv_value(v) for k, v in r.items()})

    out_tex.parent.mkdir(parents=True, exist_ok=True)
    with out_tex.open("w", encoding="utf-8") as f:
        f.write("\\begin{tabular}{llrrr}\n")
        f.write("\\hline\n")
        f.write("Method & p & Mean Acc & Min Acc & Max Acc \\\\ \n")
        f.write("\\hline\n")
        for r in rows:
            f.write(
                f"{r['method']} & {float(r['p_input']):.2f} & {to_csv_value(r['mean_accuracy'])} "
                f"& {to_csv_value(r['min_accuracy'])} & {to_csv_value(r['max_accuracy'])} \\\\ \n"
            )
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")

    return rows

def parse_int(text: str) -> Optional[int]:
    if text is None:
        return None
    text = text.strip().replace(",", "")
    if text == "":
        return None
    try:
        return int(text)
    except ValueError:
        return None


def parse_float(text: str) -> Optional[float]:
    if text is None:
        return None
    text = text.strip().replace(",", "")
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def to_csv_value(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        if math.isnan(x) or math.isinf(x):
            return ""
        return f"{x:.6f}"
    return str(x)


def infer_vivado_targets(all_ns: List[int], subset_raw: str, full_sweep: bool) -> List[int]:
    if full_sweep:
        return all_ns
    subset = []
    if subset_raw.strip():
        for tok in subset_raw.split(","):
            tok = tok.strip()
            if not tok:
                continue
            subset.append(int(tok))
    subset = sorted(set(subset))
    return [n for n in subset if n in all_ns]


def ensure_tools_built(repo_root: Path, build_dir: Path) -> None:
    build_dir.mkdir(parents=True, exist_ok=True)
    c1 = run_cmd(["cmake", "-S", str(repo_root), "-B", str(build_dir)])
    if c1.returncode != 0:
        raise RuntimeError("cmake configure failed")
    c2 = run_cmd(["cmake", "--build", str(build_dir), "-j"])
    if c2.returncode != 0:
        raise RuntimeError("cmake build failed")


def read_json(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def default_row() -> Dict[str, str]:
    return {k: "" for k in CSV_COLUMNS}


def parse_utilization_report(path: Path) -> Dict[str, Optional[int]]:
    text = path.read_text(encoding="utf-8", errors="ignore")

    def pick(patterns: Sequence[str]) -> Optional[int]:
        for pattern in patterns:
            m = re.search(pattern, text, flags=re.MULTILINE)
            if m:
                return parse_int(m.group(1))
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

    if lut is None:
        lut = 0
    if ff is None:
        ff = 0
    if carry is None:
        carry = 0

    return {
        "lut_count": lut,
        "ff_count": ff,
        "carry_count": carry,
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

    if wns is None:
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if "WNS(ns)" in line and i + 2 < len(lines):
                vals = re.findall(r"[-+]?\d+(?:\.\d+)?", lines[i + 2])
                if vals:
                    wns = parse_float(vals[0])
                    break

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


def run_vivado_one(
    repo_root: Path,
    repo_alias_root: Optional[Path],
    design: GeneratedDesign,
    part: str,
    flow_mode: str,
    period_ns: float,
    reports_dir: Path,
    vivado_bin: str,
) -> Tuple[bool, Dict[str, str]]:
    out = {
        "util_rpt": "",
        "timing_rpt": "",
        "power_rpt": "",
        "error": "",
        "status": "ok",
    }

    rpt_dir = reports_dir / design.logic_repr / design.method / f"n{design.n}" / flow_mode
    rpt_dir.mkdir(parents=True, exist_ok=True)

    util_rpt = rpt_dir / "utilization.rpt"
    timing_rpt = rpt_dir / "timing.rpt"
    power_rpt = rpt_dir / "power.rpt"
    vivado_log = rpt_dir / "vivado.log"

    def alias_path(p: Path) -> Path:
        if repo_alias_root is None:
            return p
        try:
            rel = p.absolute().relative_to(repo_root.absolute())
            return repo_alias_root / rel
        except Exception:
            return p

    cmd = [
        vivado_bin,
        "-mode",
        "batch",
        "-source",
        str(alias_path(repo_root / "scripts" / "vivado_eval.tcl")),
        "-notrace",
        "-tclargs",
        str(alias_path(design.verilog_file)),
        design.module,
        part,
        flow_mode,
        f"{period_ns}",
        str(alias_path(rpt_dir)),
    ]

    proc = run_cmd(cmd, capture=True)
    vivado_log.write_text((proc.stdout or "") + "\n" + (proc.stderr or ""), encoding="utf-8")

    out["util_rpt"] = str(util_rpt)
    out["timing_rpt"] = str(timing_rpt)
    out["power_rpt"] = str(power_rpt)

    if proc.returncode != 0:
        out["status"] = "error"
        out["error"] = f"vivado failed rc={proc.returncode}"
        return False, out

    missing = [p for p in (util_rpt, timing_rpt, power_rpt) if not p.exists()]
    if missing:
        out["status"] = "error"
        out["error"] = "missing reports: " + ", ".join(str(p) for p in missing)
        return False, out

    return True, out


def pick_best_rows_for_plot(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    by_key: Dict[Tuple[str, str, int], Dict[str, str]] = {}
    priority = {"post_route": 3, "synth_only": 2, "ntk_only": 1}

    for r in rows:
        try:
            key = (r["method"], row_logic_repr(r), int(r["n"]))
        except Exception:
            continue
        p = priority.get(r.get("flow_mode", ""), 0)
        if key not in by_key or p > priority.get(by_key[key].get("flow_mode", ""), 0):
            by_key[key] = r

    return list(by_key.values())


def rows_by_method(rows: List[Dict[str, str]]) -> Dict[str, List[Dict[str, str]]]:
    out: Dict[str, List[Dict[str, str]]] = {}
    for r in rows:
        out.setdefault(r["method"], []).append(r)
    for k in out:
        out[k].sort(key=lambda x: int(x["n"]))
    return out


def rows_by_variant(rows: List[Dict[str, str]]) -> Dict[str, List[Dict[str, str]]]:
    out: Dict[str, List[Dict[str, str]]] = {}
    for r in rows:
        out.setdefault(row_variant_label(r), []).append(r)
    for k in out:
        out[k].sort(key=lambda x: int(x["n"]))
    return out


def get_numeric(r: Dict[str, str], key: str) -> Optional[float]:
    v = r.get(key, "")
    if v is None or v == "":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def save_plot_metric(rows: List[Dict[str, str]], metric: str, ylabel: str, title: str, out_prefix: Path) -> None:
    plt.figure(figsize=(9, 5))
    grouped = rows_by_variant(rows)

    drawn = False
    for label, rs in grouped.items():
        xs: List[int] = []
        ys: List[float] = []
        for r in rs:
            y = get_numeric(r, metric)
            if y is None:
                continue
            xs.append(int(r["n"]))
            ys.append(y)
        if xs:
            plt.plot(xs, ys, marker="o", linewidth=1.5, markersize=3, label=label)
            drawn = True

    if not drawn:
        plt.close()
        return

    plt.xlabel("n")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_prefix.with_suffix(".png")), dpi=180)
    plt.savefig(str(out_prefix.with_suffix(".pdf")))
    plt.close()


def save_depth_plot(rows: List[Dict[str, str]], out_prefix: Path) -> None:
    plt.figure(figsize=(9, 5))
    grouped = rows_by_variant(rows)

    for label, rs in grouped.items():
        method = rs[0].get("method", "") if rs else ""
        xs: List[int] = []
        ys_logic: List[float] = []
        for r in rs:
            n = int(r["n"])
            d_logic = get_numeric(r, "ntk_depth_logic")
            if d_logic is not None:
                xs.append(n)
                ys_logic.append(d_logic)

        if xs and ys_logic:
            plt.plot(xs, ys_logic, marker="o", linewidth=1.5, markersize=3, label=f"{label} logic depth")
        if method in {"cancel_tree", "cancel_tree_v2", "boyermoore_tree", "boyermoore_tree_v2"}:
            xs2 = [int(r["n"]) for r in rs if get_numeric(r, "cancel_merge_depth") is not None]
            ys2 = [get_numeric(r, "cancel_merge_depth") for r in rs if get_numeric(r, "cancel_merge_depth") is not None]
            if xs2 and ys2:
                plt.plot(xs2, ys2, linestyle="--", linewidth=1.5, label=f"{label} merge depth")

    plt.xlabel("n")
    plt.ylabel("depth")
    plt.title("Depth vs n")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_prefix.with_suffix(".png")), dpi=180)
    plt.savefig(str(out_prefix.with_suffix(".pdf")))
    plt.close()


def save_area_delay_scatter(rows: List[Dict[str, str]], out_prefix: Path) -> None:
    plt.figure(figsize=(7, 6))
    grouped = rows_by_variant(rows)

    drawn = False
    for label, rs in grouped.items():
        xs = []
        ys = []
        for r in rs:
            lut = get_numeric(r, "lut_count")
            delay = get_numeric(r, "crit_delay_ns")
            if lut is None or delay is None:
                continue
            xs.append(lut)
            ys.append(delay)
        if xs:
            plt.scatter(xs, ys, label=label, s=25, alpha=0.8)
            drawn = True

    if not drawn:
        plt.close()
        return

    plt.xlabel("LUT count")
    plt.ylabel("Critical delay (ns)")
    plt.title("Area-Delay Tradeoff")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_prefix.with_suffix(".png")), dpi=180)
    plt.savefig(str(out_prefix.with_suffix(".pdf")))
    plt.close()


def get_area_metric_value(row: Dict[str, str]) -> Tuple[Optional[str], Optional[float]]:
    lut = get_numeric(row, "lut_count")
    if lut is not None:
        return "lut_count", lut

    logic_repr = row_logic_repr(row)
    if logic_repr == "aig":
        aig_and = get_numeric(row, "aig_and_count")
        if aig_and is not None:
            return "aig_and_count", aig_and
    elif logic_repr == "mig":
        mig_maj = get_numeric(row, "mig_maj_count")
        if mig_maj is not None:
            return "mig_maj_count", mig_maj

    and_count = get_numeric(row, "xag_and_count")
    xor_count = get_numeric(row, "xag_xor_count")
    if and_count is not None and xor_count is not None:
        return "xag_gate_total", and_count + xor_count

    aig_and = get_numeric(row, "aig_and_count")
    if aig_and is not None:
        return "aig_and_count", aig_and

    mig_maj = get_numeric(row, "mig_maj_count")
    if mig_maj is not None:
        return "mig_maj_count", mig_maj

    return None, None


def build_half_resource_summary_rows(rows: List[Dict[str, str]], baseline_method: str) -> List[Dict[str, str]]:
    grouped = rows_by_variant(rows)
    base_map: Dict[Tuple[str, int], Dict[str, str]] = {}
    for r in rows:
        if r.get("method") != baseline_method:
            continue
        base_map[(row_logic_repr(r), int(r["n"]))] = r

    out_rows: List[Dict[str, str]] = []

    for _, rs in grouped.items():
        if not rs:
            continue
        method = rs[0].get("method", "")
        logic_repr = row_logic_repr(rs[0])
        if method == baseline_method:
            continue

        ratios: List[float] = []
        lut_points = 0
        proxy_points = 0

        for r in rs:
            n = int(r["n"])
            b = base_map.get((logic_repr, n))
            if b is None:
                continue

            if get_numeric(r, "lut_count") is not None and get_numeric(b, "lut_count") is not None:
                a = float(get_numeric(r, "lut_count"))
                bb = float(get_numeric(b, "lut_count"))
                metric = "lut_count"
            else:
                a_metric, a = get_area_metric_value(r)
                b_metric, bb = get_area_metric_value(b)
                if a is None or bb is None or a_metric != b_metric:
                    continue
                metric = a_metric

            if bb <= 0:
                continue

            ratio = a / bb
            ratios.append(ratio)
            if metric == "lut_count":
                lut_points += 1
            else:
                proxy_points += 1

        if not ratios:
            continue

        out_rows.append(
            {
                "method": method,
                "logic_repr": logic_repr,
                "points_total": len(ratios),
                "points_lut": lut_points,
                "points_proxy": proxy_points,
                "mean_area_ratio": statistics.mean(ratios),
                "min_area_ratio": min(ratios),
                "max_area_ratio": max(ratios),
                "fraction_le_half": sum(1 for v in ratios if v <= 0.5) / float(len(ratios)),
            }
        )

    out_rows.sort(key=lambda r: r["method"])
    return out_rows


def write_half_resource_tables(rows: List[Dict[str, str]], out_csv: Path, out_tex: Path) -> None:
    keys = [
        "method",
        "logic_repr",
        "points_total",
        "points_lut",
        "points_proxy",
        "mean_area_ratio",
        "min_area_ratio",
        "max_area_ratio",
        "fraction_le_half",
    ]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: to_csv_value(r.get(k, "")) for k in keys})

    out_tex.parent.mkdir(parents=True, exist_ok=True)
    with out_tex.open("w", encoding="utf-8") as f:
        f.write("\\begin{tabular}{llrrrrrrr}\n")
        f.write("\\hline\n")
        f.write("Method & Repr & Points & LUT Pts & Proxy Pts & Mean Ratio & Min Ratio & Max Ratio & Frac $\\le 0.5$ \\\\ \n")
        f.write("\\hline\n")
        for r in rows:
            f.write(
                f"{r['method']} & {r['logic_repr']} & {r['points_total']} & {r['points_lut']} & {r['points_proxy']} "
                f"& {to_csv_value(r['mean_area_ratio'])} & {to_csv_value(r['min_area_ratio'])} "
                f"& {to_csv_value(r['max_area_ratio'])} & {to_csv_value(r['fraction_le_half'])} \\\\ \n"
            )
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")


def save_area_ratio_plot(rows: List[Dict[str, str]], baseline_method: str, out_prefix: Path) -> None:
    grouped = rows_by_variant(rows)
    baseline: Dict[Tuple[str, int], Dict[str, str]] = {}
    for r in rows:
        if r.get("method") != baseline_method:
            continue
        baseline[(row_logic_repr(r), int(r["n"]))] = r

    plt.figure(figsize=(9, 5))
    drawn = False
    for label, rs in grouped.items():
        if not rs:
            continue
        method = rs[0].get("method", "")
        logic_repr = row_logic_repr(rs[0])
        if method == baseline_method:
            continue
        xs: List[int] = []
        ys: List[float] = []
        for r in rs:
            n = int(r["n"])
            b = baseline.get((logic_repr, n))
            if b is None:
                continue
            if get_numeric(r, "lut_count") is not None and get_numeric(b, "lut_count") is not None:
                a = float(get_numeric(r, "lut_count"))
                bb = float(get_numeric(b, "lut_count"))
            else:
                ma, a = get_area_metric_value(r)
                mb, bb = get_area_metric_value(b)
                if ma != mb or a is None or bb is None:
                    continue
            if bb <= 0:
                continue
            xs.append(n)
            ys.append(a / bb)
        if xs:
            plt.plot(xs, ys, marker="o", linewidth=1.5, markersize=3, label=label)
            drawn = True

    if not drawn:
        plt.close()
        return

    plt.axhline(0.5, color="red", linestyle="--", linewidth=1.2, label="half-resource target")
    plt.xlabel("n")
    plt.ylabel("Area Ratio vs Baseline")
    plt.title(f"Area Ratio vs {baseline_method}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_prefix.with_suffix(".png")), dpi=180)
    plt.savefig(str(out_prefix.with_suffix(".pdf")))
    plt.close()


def geometric_mean(values: List[float]) -> Optional[float]:
    vals = [v for v in values if v > 0]
    if not vals:
        return None
    return math.exp(sum(math.log(v) for v in vals) / len(vals))


def build_tables(rows: List[Dict[str, str]], out_csv: Path, out_tex: Path, baseline_method: str) -> None:
    buckets = [(3, 63), (65, 127), (129, 255), (257, 511)]
    grouped = rows_by_variant(rows)

    table_rows = []

    # Pre-index baseline by (logic_repr, n) for fair speedup comparisons.
    baseline_by_key: Dict[Tuple[str, int], Dict[str, str]] = {}
    baseline_rows = [r for r in rows if r.get("method") == baseline_method]
    for r in baseline_rows:
        baseline_by_key[(row_logic_repr(r), int(r["n"]))] = r

    for lo, hi in buckets:
        for _, rs in grouped.items():
            if not rs:
                continue
            logic_repr = row_logic_repr(rs[0])
            method = rs[0].get("method", "")
            in_range = [r for r in rs if lo <= int(r["n"]) <= hi]
            luts = [get_numeric(r, "lut_count") for r in in_range]
            delays = [get_numeric(r, "crit_delay_ns") for r in in_range]
            luts = [x for x in luts if x is not None]
            delays = [x for x in delays if x is not None]

            mean_lut = statistics.mean(luts) if luts else None
            max_lut = max(luts) if luts else None
            mean_delay = statistics.mean(delays) if delays else None

            speedup_terms = []
            for r in in_range:
                n = int(r["n"])
                d_method = get_numeric(r, "crit_delay_ns")
                d_base = get_numeric(baseline_by_key.get((logic_repr, n), {}), "crit_delay_ns")
                if d_method is None or d_base is None or d_method <= 0:
                    continue
                speedup_terms.append(d_base / d_method)
            gm_speedup = geometric_mean(speedup_terms)

            baseline_range = [
                r
                for r in baseline_rows
                if row_logic_repr(r) == logic_repr and lo <= int(r["n"]) <= hi
            ]
            baseline_luts = [get_numeric(r, "lut_count") for r in baseline_range]
            baseline_luts = [x for x in baseline_luts if x is not None]
            base_mean_lut = statistics.mean(baseline_luts) if baseline_luts else None

            lut_reduction = None
            if mean_lut is not None and base_mean_lut is not None and base_mean_lut > 0:
                lut_reduction = 1.0 - (mean_lut / base_mean_lut)

            table_rows.append(
                {
                    "range": f"{lo}-{hi}",
                    "method": method,
                    "logic_repr": logic_repr,
                    "mean_lut": mean_lut,
                    "max_lut": max_lut,
                    "mean_delay": mean_delay,
                    "geom_speedup_vs_baseline": gm_speedup,
                    "lut_reduction_ratio": lut_reduction,
                }
            )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(table_rows[0].keys()) if table_rows else ["range", "method", "logic_repr"])
        writer.writeheader()
        for r in table_rows:
            writer.writerow({k: to_csv_value(v) for k, v in r.items()})

    out_tex.parent.mkdir(parents=True, exist_ok=True)
    with out_tex.open("w", encoding="utf-8") as f:
        f.write("\\begin{tabular}{lllrrrrr}\n")
        f.write("\\hline\n")
        f.write("Range & Method & Repr & Mean LUT & Max LUT & Mean Delay (ns) & Geo Speedup & LUT Reduction \\\\ \n")
        f.write("\\hline\n")
        for r in table_rows:
            f.write(
                f"{r['range']} & {r['method']} & {r['logic_repr']} & {to_csv_value(r['mean_lut'])} & {to_csv_value(r['max_lut'])} "
                f"& {to_csv_value(r['mean_delay'])} & {to_csv_value(r['geom_speedup_vs_baseline'])} "
                f"& {to_csv_value(r['lut_reduction_ratio'])} \\\\ \n"
            )
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")


def get_tool_version(cmd: Sequence[str]) -> str:
    proc = run_cmd(cmd, capture=True)
    txt = (proc.stdout or "") + "\n" + (proc.stderr or "")
    for line in txt.splitlines():
        s = line.strip()
        if s:
            return s
    if proc.returncode != 0:
        return "unavailable"
    return "unavailable"


def find_global_stats(rows: List[Dict[str, str]], method: str, key: str) -> List[float]:
    vals = []
    for r in rows:
        if r.get("method") != method:
            continue
        v = get_numeric(r, key)
        if v is not None:
            vals.append(v)
    return vals


def write_report(
    repo_root: Path,
    out_md: Path,
    selected_rows: List[Dict[str, str]],
    all_rows: List[Dict[str, str]],
    baseline_method: str,
    part: str,
    period_ns: float,
    vivado_bin: str,
    accuracy_summary_rows: Optional[List[Dict[str, str]]] = None,
    half_resource_rows: Optional[List[Dict[str, str]]] = None,
) -> None:
    accuracy_summary_rows = accuracy_summary_rows or []
    half_resource_rows = half_resource_rows or []
    methods = sorted({r.get("method", "") for r in selected_rows if r.get("method", "")})
    logic_reprs = sorted({row_logic_repr(r) for r in selected_rows if row_logic_repr(r)})
    obs_lines = []

    grouped = rows_by_variant(selected_rows)
    baseline_stats_by_repr: Dict[str, Dict[str, Optional[float]]] = {}
    for logic_repr in logic_reprs:
        base_rows = [r for r in selected_rows if r.get("method") == baseline_method and row_logic_repr(r) == logic_repr]
        base_luts = [v for v in [get_numeric(r, "lut_count") for r in base_rows] if v is not None]
        base_delays = [v for v in [get_numeric(r, "crit_delay_ns") for r in base_rows] if v is not None]
        base_powers = [v for v in [get_numeric(r, "power_w") for r in base_rows] if v is not None]
        base_lut = statistics.mean(base_luts) if base_luts else None
        base_delay = statistics.mean(base_delays) if base_delays else None
        base_power = statistics.mean(base_powers) if base_powers else None
        baseline_stats_by_repr[logic_repr] = {
            "lut": base_lut,
            "delay": base_delay,
            "power": base_power,
        }

    for label, rs in grouped.items():
        if not rs:
            continue
        method = rs[0].get("method", "")
        logic_repr = row_logic_repr(rs[0])
        if method == baseline_method:
            continue
        luts = [v for v in [get_numeric(r, "lut_count") for r in rs] if v is not None]
        delays = [v for v in [get_numeric(r, "crit_delay_ns") for r in rs] if v is not None]
        powers = [v for v in [get_numeric(r, "power_w") for r in rs] if v is not None]

        base_stats = baseline_stats_by_repr.get(logic_repr, {})
        base_label = variant_label(baseline_method, logic_repr)
        base_lut = base_stats.get("lut")
        base_delay = base_stats.get("delay")
        base_power = base_stats.get("power")

        if luts and base_lut is not None and base_lut > 0:
            rel = 100.0 * (1.0 - statistics.mean(luts) / base_lut)
            obs_lines.append(f"- `{label}` vs `{base_label}` mean LUT delta: {rel:.2f}% (positive means smaller).")
        if delays and base_delay is not None and base_delay > 0:
            rel = 100.0 * (1.0 - statistics.mean(delays) / base_delay)
            obs_lines.append(f"- `{label}` vs `{base_label}` mean delay delta: {rel:.2f}% (positive means faster).")
        if powers and base_power is not None and base_power > 0:
            rel = 100.0 * (1.0 - statistics.mean(powers) / base_power)
            obs_lines.append(f"- `{label}` vs `{base_label}` mean power delta: {rel:.2f}% (positive means lower power).")

    for row in half_resource_rows:
        row_label = variant_label(str(row.get("method", "")), str(row.get("logic_repr", "")))
        obs_lines.append(
            f"- `{row_label}` half-resource hit rate (area ratio <= 0.5): {100.0 * float(row['fraction_le_half']):.2f}% "
            f"over {row['points_total']} n-points."
        )

    # Accuracy at p=0.5 is the key stress case for majority approximation.
    p050_rows = [r for r in accuracy_summary_rows if abs(float(r["p_input"]) - 0.5) < 1e-9]
    for r in p050_rows:
        obs_lines.append(
            f"- Accuracy at p=0.50 for `{r['method']}`: mean={float(r['mean_accuracy']):.4f}, "
            f"min={float(r['min_accuracy']):.4f}, max={float(r['max_accuracy']):.4f}."
        )

    all_wns = [v for r in selected_rows for v in [get_numeric(r, "wns_ns")] if v is not None]
    if all_wns:
        obs_lines.append(f"- WNS range across selected rows: [{min(all_wns):.3f}, {max(all_wns):.3f}] ns.")

    if not obs_lines:
        obs_lines.append("- No timing/area/accuracy observations available yet; run full sweep options to populate results.")

    # Reproducibility metadata.
    git_commit = "N/A (no git repo)"
    if (repo_root / ".git").exists():
        proc = run_cmd(["git", "rev-parse", "HEAD"], cwd=repo_root, capture=True)
        if proc.returncode == 0:
            git_commit = (proc.stdout or "").strip() or git_commit

    tool_versions = {
        "majgen": get_tool_version([str(repo_root / "build" / "majgen"), "--version"]),
        "vivado": get_tool_version([vivado_bin, "-version"]),
        "yosys": get_tool_version(["yosys", "-V"]),
        "abc": get_tool_version(["abc", "-h"]),
    }

    out_md.parent.mkdir(parents=True, exist_ok=True)
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Majority Generator Evaluation\n\n")
        f.write("## Method\n")
        f.write(
            "This run compares exact and/or approximate odd-`n` majority generators using a shared flow "
            "(Mockturtle network generation, optional Vivado OOC synth/post-route, and uniform CSV/plot/reporting outputs).\n\n"
        )
        if any(is_approx_block3_method(m) for m in methods):
            f.write(
                "Approximation used here: `approx_block3[_popcount]`, which applies one compression stage of 3-input majority "
                "(`maj3`) on disjoint triples, forwards tail bits unchanged, and then applies exact majority on the compressed vector.\n\n"
            )
        if any(m in {"cancel_tree", "boyermoore_tree", "cancel_tree_v2", "boyermoore_tree_v2"} for m in methods):
            f.write(
                "Exact cancellation variants (`cancel_tree`, `cancel_tree_v2`) implement Boyer-Moore summary cancellation in a balanced tree.\n\n"
            )

        f.write("## Process\n")
        f.write("- Sweep range: odd n values in requested interval.\n")
        f.write("- Methods: " + ", ".join(f"`{m}`" for m in methods) + "\n")
        if logic_reprs:
            f.write("- Logic representations: " + ", ".join(f"`{lr}`" for lr in logic_reprs) + "\n")
        f.write("- Baseline method for ratios/tables: " + f"`{baseline_method}`\n")
        f.write("- Area ratio target tracked: method_area / baseline_area <= 0.5\n\n")

        f.write("## Reproducibility\n")
        f.write(f"- Date: {dt.datetime.now().isoformat()}\n")
        f.write(f"- Repo: {repo_root}\n")
        f.write(f"- Git commit: {git_commit}\n")
        f.write(f"- FPGA part: `{part}`\n")
        f.write(f"- Timing constraint: `set_max_delay {period_ns} -from [all_inputs] -to [all_outputs]`\n")
        for k, v in tool_versions.items():
            f.write(f"- {k}: {v}\n")
        f.write("\n")

        f.write("## Outputs\n")
        f.write("- `master_results.csv` (generation + timing/area/power rows)\n")
        f.write("- `figures/lut_vs_n.(png|pdf)`\n")
        f.write("- `figures/delay_vs_n.(png|pdf)`\n")
        f.write("- `figures/power_vs_n.(png|pdf)`\n")
        f.write("- `figures/depth_vs_n.(png|pdf)`\n")
        f.write("- `figures/area_delay_scatter.(png|pdf)`\n")
        f.write("- `figures/area_ratio_vs_baseline.(png|pdf)`\n")
        f.write("- `tables_cancel_tree.csv`, `tables_cancel_tree.tex`\n")
        f.write("- `tables_half_resource.csv`, `tables_half_resource.tex`\n")
        f.write("- `accuracy_results.csv` (if `--run-accuracy`)\n")
        f.write("- `tables_accuracy.csv`, `tables_accuracy.tex` (if `--run-accuracy`)\n")
        f.write("- `figures/accuracy_vs_n_pXX.(png|pdf)` and `figures/accuracy_vs_p.(png|pdf)` (if `--run-accuracy`)\n\n")

        f.write("## Observations\n")
        for line in obs_lines:
            f.write(f"{line}\n")
        f.write("\n")

        f.write("## TODO\n")
        f.write("- Carry-chain friendly compare/sub implementations.\n")
        f.write("- Optional multi-stage approximation space exploration.\n")
        f.write("- Optional pipelining for Fmax-focused experiments.\n")

        f.write("\n## Cancel-Tree Background\n")
        f.write(
            "The `cancel_tree` method builds a balanced merge tree of Boyer-Moore style summaries `(cand, k)`, "
            "where `cand` is a 1-bit candidate and `k` is the surplus count after pair cancellation. "
            "Each merge follows the exact monoid combine rules with empty-summary priority and deterministic zero-candidate handling for `k=0`. "
            "Widths are assigned as `w(sz)=ceil(log2(sz+1))` and both children are zero-extended to the merge width before add/sub/compare.\n\n"
        )


def maybe_generate_pdf(markdown_path: Path, pdf_path: Path) -> None:
    if shutil.which("pandoc") is None:
        return
    proc = run_cmd(["pandoc", str(markdown_path), "-o", str(pdf_path)], capture=True)
    if proc.returncode != 0:
        # Keep markdown as authoritative output if PDF conversion fails.
        return


def enrich_and_write_csv(rows: List[Dict[str, str]], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: to_csv_value(r.get(k, "")) for k in CSV_COLUMNS})


def main() -> int:
    parser = argparse.ArgumentParser(description="Run cancel_tree generation, evaluation, and reporting flow")
    parser.add_argument("--build-dir", default="build")
    parser.add_argument("--out-dir", default="results/cancel_tree_eval")
    parser.add_argument("--n-min", type=int, default=5)
    parser.add_argument("--n-max", type=int, default=511)
    parser.add_argument("--methods", default="popcount,approx_block3")
    parser.add_argument("--logic-reprs", default="xag")
    parser.add_argument("--part", default="xc7a200tfbg484-1")
    parser.add_argument("--period-ns", type=float, default=5.0)
    parser.add_argument("--run-vivado", action="store_true")
    parser.add_argument("--vivado-bin", default=os.environ.get("VIVADO_BIN", "vivado"))
    parser.add_argument("--flow-modes", default="synth_only,post_route")
    parser.add_argument("--vivado-subset", default="63,127,255,511")
    parser.add_argument("--full-sweep-vivado", action="store_true")
    parser.add_argument("--run-accuracy", action="store_true")
    parser.add_argument("--accuracy-ps", default="0.10,0.25,0.50,0.55,0.75,0.90,0.045")
    parser.add_argument("--accuracy-trials", type=int, default=4000)
    parser.add_argument("--accuracy-seed", type=int, default=20260221)
    parser.add_argument("--accuracy-batch-size", type=int, default=2000)
    parser.add_argument("--skip-build", action="store_true")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    build_dir = (repo_root / args.build_dir).resolve()

    methods = parse_csv_list(args.methods)
    if not methods:
        raise RuntimeError("No methods selected")
    logic_reprs = parse_csv_list(args.logic_reprs)
    if not logic_reprs:
        raise RuntimeError("No logic representations selected")

    valid_methods = {
        "popcount",
        "popcount_strict",
        "baseline_strict",
        "cancel_tree",
        "boyermoore_tree",
        "cancel_tree_v2",
        "boyermoore_tree_v2",
        "approx_block3",
        "approx_block3_popcount",
    }
    for m in methods:
        if m not in valid_methods:
            raise RuntimeError(f"Unsupported method '{m}'")
    valid_logic_reprs = {"xag", "aig", "mig"}
    for lr in logic_reprs:
        if lr not in valid_logic_reprs:
            raise RuntimeError(f"Unsupported logic representation '{lr}'")

    if args.run_accuracy and np is None:
        raise RuntimeError("numpy is required for --run-accuracy")
    if args.accuracy_trials <= 0:
        raise RuntimeError("--accuracy-trials must be > 0")
    if args.accuracy_batch_size <= 0:
        raise RuntimeError("--accuracy-batch-size must be > 0")
    accuracy_ps = parse_probability_list(args.accuracy_ps) if args.run_accuracy else []

    flow_modes = [m.strip() for m in args.flow_modes.split(",") if m.strip()]
    if args.run_vivado:
        for fm in flow_modes:
            if fm not in {"synth_only", "post_route"}:
                raise RuntimeError(f"Unsupported flow mode '{fm}'")
        vivado_probe = run_cmd([args.vivado_bin, "-version"], capture=True)
        if vivado_probe.returncode != 0:
            raise RuntimeError(
                f"Vivado not available via --vivado-bin '{args.vivado_bin}'. "
                f"stderr: {(vivado_probe.stderr or '').strip()}"
            )

    if not args.skip_build:
        ensure_tools_built(repo_root, build_dir)

    majgen = build_dir / "majgen"
    if not majgen.exists():
        raise RuntimeError(f"majgen not found at {majgen}")

    now = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = (repo_root / args.out_dir / f"run_{now}").resolve()
    netlist_dir = out_dir / "netlists"
    stats_dir = out_dir / "stats"
    reports_dir = out_dir / "reports"
    figures_dir = out_dir / "figures"

    for d in [out_dir, netlist_dir, stats_dir, reports_dir, figures_dir]:
        d.mkdir(parents=True, exist_ok=True)

    ns = odd_values(args.n_min, args.n_max)
    vivado_targets = infer_vivado_targets(ns, args.vivado_subset, args.full_sweep_vivado)

    repo_alias_root: Optional[Path] = None
    if args.run_vivado:
        # Vivado's -tclargs handling can break on paths containing spaces.
        # Create a temporary symlink alias with no spaces and pass alias paths.
        repo_alias_root = Path(f"/tmp/cancel_tree_repo_alias_{os.getpid()}")
        if repo_alias_root.exists() or repo_alias_root.is_symlink():
            repo_alias_root.unlink()
        repo_alias_root.symlink_to(repo_root, target_is_directory=True)

    print(f"[INFO] Output directory: {out_dir}")
    print(f"[INFO] Methods: {methods}")
    print(f"[INFO] Logic representations: {logic_reprs}")
    print(f"[INFO] n sweep: {ns[0]}..{ns[-1]} ({len(ns)} odd values)")
    if args.run_vivado:
        print(f"[INFO] Vivado flow modes: {flow_modes}")
        print(f"[INFO] Vivado targets: {vivado_targets}")
    if args.run_accuracy:
        print(f"[INFO] Accuracy sweep p values: {accuracy_ps}")
        print(f"[INFO] Accuracy trials per (method,n,p): {args.accuracy_trials}")

    generated: List[GeneratedDesign] = []
    rows: List[Dict[str, str]] = []

    # Generation phase.
    for logic_repr in logic_reprs:
        for method in methods:
            for n in ns:
                module = f"majority_{method}_{logic_repr}_n{n}"
                verilog_file = netlist_dir / f"maj_{n}_{method}_{logic_repr}.v"
                stats_json = stats_dir / f"maj_{n}_{method}_{logic_repr}.json"

                cmd = [
                    str(majgen),
                    "--maj_impl",
                    method,
                    "--logic_repr",
                    logic_repr,
                    "--n",
                    str(n),
                    "--out",
                    str(verilog_file),
                    "--format",
                    "verilog",
                    "--module",
                    module,
                    "--stats_json",
                    str(stats_json),
                ]

                # Required debug DOT sanity artifact.
                if method == "cancel_tree" and n == 7:
                    if logic_repr == "xag":
                        cmd.extend(["--dot", str(out_dir / "artifacts" / "maj_7_cancel_tree.dot")])
                    else:
                        cmd.extend(["--dot", str(out_dir / "artifacts" / f"maj_7_cancel_tree_{logic_repr}.dot")])

                proc = run_cmd(cmd, capture=True)
                if proc.returncode != 0:
                    r = default_row()
                    r.update(
                        {
                            "method": method,
                            "logic_repr": logic_repr,
                            "n": str(n),
                            "module": module,
                            "file": str(verilog_file),
                            "part": args.part,
                            "flow_mode": "ntk_only",
                            "status": "error",
                            "error": f"majgen rc={proc.returncode}: {(proc.stderr or '').strip()}",
                        }
                    )
                    rows.append(r)
                    continue

                generated.append(
                    GeneratedDesign(
                        method=method,
                        logic_repr=logic_repr,
                        n=n,
                        module=module,
                        verilog_file=verilog_file,
                        stats_json=stats_json,
                    )
                )

    # Enrichment phase.
    for design in generated:
        stats = read_json(design.stats_json)

        base_row = default_row()
        base_row.update(
            {
                "method": design.method,
                "logic_repr": str(stats.get("logic_repr", design.logic_repr)),
                "n": str(design.n),
                "module": design.module,
                "file": str(design.verilog_file),
                "part": args.part,
                "fa_count": "",
                "mig_maj_count": str(stats.get("mig_maj_count", "")),
                "aig_and_count": str(stats.get("aig_and_count", "")),
                "inv_count": str(stats.get("inv_count", "")),
                "popcount_levels": str(stats.get("popcount_levels", "")),
                "strict_schedule_mode": str(stats.get("strict_schedule_mode", "")),
                "strict_scaffold_p": str(stats.get("strict_scaffold_p", "")),
                "strict_scaffold_inputs": str(stats.get("strict_scaffold_inputs", "")),
                "strict_scaffold_threshold": str(stats.get("strict_scaffold_threshold", "")),
                "strict_comparator_width": str(stats.get("strict_comparator_width", "")),
                "strict_num_fixed_pairs": str(stats.get("strict_num_fixed_pairs", "")),
                "strict_csa_fa_count": str(stats.get("strict_csa_fa_count", "")),
                "strict_comparator_fa_count": str(stats.get("strict_comparator_fa_count", "")),
                "strict_total_fa_count": str(stats.get("strict_total_fa_count", "")),
                "strict_csa_levels": str(stats.get("strict_csa_levels", "")),
                "strict_total_levels": str(stats.get("strict_total_levels", "")),
                "cell_count": str(stats.get("node_count", "")),
                "cancel_merge_depth": str(stats.get("cancel_merge_depth", "")),
                "cancel_max_k_width": str(stats.get("cancel_max_k_width", "")),
                "cancel_num_merges": str(stats.get("cancel_num_merges", "")),
                "xag_xor_count": str(stats.get("xag_xor_count", "")),
                "xag_and_count": str(stats.get("xag_and_count", "")),
                "ntk_depth_logic": str(stats.get("ntk_depth_logic", "")),
                "approx_group_size": str(stats.get("approx_group_size", "")),
                "approx_compressed_n": str(stats.get("approx_compressed_n", "")),
                "approx_groups3": str(stats.get("approx_groups3", "")),
                "approx_passthrough_bits": str(stats.get("approx_passthrough_bits", "")),
            }
        )

        do_vivado = args.run_vivado and design.n in vivado_targets
        if not do_vivado:
            r = dict(base_row)
            r.update(
                {
                    "flow_mode": "ntk_only",
                    "status": "ok",
                }
            )
            rows.append(r)
            continue

        for flow_mode in flow_modes:
            ok, vivado_info = run_vivado_one(
                repo_root=repo_root,
                repo_alias_root=repo_alias_root,
                design=design,
                part=args.part,
                flow_mode=flow_mode,
                period_ns=args.period_ns,
                reports_dir=reports_dir,
                vivado_bin=args.vivado_bin,
            )

            r = dict(base_row)
            r["flow_mode"] = flow_mode
            r["status"] = vivado_info["status"]
            r["error"] = vivado_info["error"]
            r["util_rpt"] = vivado_info["util_rpt"]
            r["timing_rpt"] = vivado_info["timing_rpt"]
            r["power_rpt"] = vivado_info["power_rpt"]

            if ok:
                util = parse_utilization_report(Path(vivado_info["util_rpt"]))
                timing = parse_timing_report(Path(vivado_info["timing_rpt"]), args.period_ns)
                power = parse_power_report(Path(vivado_info["power_rpt"]))

                r["lut_count"] = to_csv_value(util["lut_count"])
                r["ff_count"] = to_csv_value(util["ff_count"])
                r["carry_count"] = to_csv_value(util["carry_count"])
                r["wns_ns"] = to_csv_value(timing["wns_ns"])
                r["crit_delay_ns"] = to_csv_value(timing["crit_delay_ns"])
                r["power_w"] = to_csv_value(power["power_w"])

            rows.append(r)

    # Persist master CSV first.
    master_csv = out_dir / "master_results.csv"
    enrich_and_write_csv(rows, master_csv)
    print(f"[INFO] Wrote CSV: {master_csv}")

    # Select preferred rows for plotting/reporting (post_route > synth_only > ntk_only).
    selected = pick_best_rows_for_plot([r for r in rows if r.get("status") == "ok"])

    save_plot_metric(selected, "lut_count", "LUT count", "LUT vs n", figures_dir / "lut_vs_n")
    save_plot_metric(selected, "crit_delay_ns", "Critical delay (ns)", "Critical Delay vs n", figures_dir / "delay_vs_n")
    save_plot_metric(selected, "power_w", "Power (W)", "Power vs n", figures_dir / "power_vs_n")
    save_depth_plot(selected, figures_dir / "depth_vs_n")
    save_area_delay_scatter(selected, figures_dir / "area_delay_scatter")

    baseline = "popcount" if "popcount" in methods else methods[0]
    save_area_ratio_plot(selected, baseline, figures_dir / "area_ratio_vs_baseline")

    tables_csv = out_dir / "tables_cancel_tree.csv"
    tables_tex = out_dir / "tables_cancel_tree.tex"
    build_tables(selected, tables_csv, tables_tex, baseline)

    half_rows = build_half_resource_summary_rows(selected, baseline)
    half_csv = out_dir / "tables_half_resource.csv"
    half_tex = out_dir / "tables_half_resource.tex"
    write_half_resource_tables(half_rows, half_csv, half_tex)

    accuracy_rows: List[Dict[str, str]] = []
    accuracy_summary_rows: List[Dict[str, str]] = []
    if args.run_accuracy:
        accuracy_rows = run_accuracy_sweep(
            methods=methods,
            ns=ns,
            ps=accuracy_ps,
            trials=args.accuracy_trials,
            seed=args.accuracy_seed,
            batch_size=args.accuracy_batch_size,
        )
        accuracy_csv = out_dir / "accuracy_results.csv"
        write_accuracy_csv(accuracy_rows, accuracy_csv)
        save_accuracy_vs_n_plots(accuracy_rows, figures_dir)
        save_accuracy_vs_p_plot(accuracy_rows, figures_dir / "accuracy_vs_p")
        accuracy_summary_rows = build_accuracy_tables(
            accuracy_rows=accuracy_rows,
            out_csv=out_dir / "tables_accuracy.csv",
            out_tex=out_dir / "tables_accuracy.tex",
        )
        print(f"[INFO] Accuracy CSV: {accuracy_csv}")

    report_md = out_dir / "report_cancel_tree.md"
    report_pdf = out_dir / "report_cancel_tree.pdf"
    write_report(
        repo_root,
        report_md,
        selected,
        rows,
        baseline,
        args.part,
        args.period_ns,
        args.vivado_bin,
        accuracy_summary_rows=accuracy_summary_rows,
        half_resource_rows=half_rows,
    )
    maybe_generate_pdf(report_md, report_pdf)

    print(f"[INFO] Figures: {figures_dir}")
    print(f"[INFO] Tables: {tables_csv}, {tables_tex}")
    print(f"[INFO] Half-resource tables: {half_csv}, {half_tex}")
    print(f"[INFO] Report: {report_md}" + (f", {report_pdf}" if report_pdf.exists() else ""))

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
