#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

N_MIN="${1:-5}"
N_MAX="${2:-65}"

METHODS="${METHODS:-popcount,popcount_strict,approx_block3}"
LOGIC_REPRS="${LOGIC_REPRS:-xag,aig,mig}"
FLOW_MODES="${FLOW_MODES:-synth_only,post_route}"

RUN_VIVADO="${RUN_VIVADO:-1}"
RUN_ACCURACY="${RUN_ACCURACY:-1}"

VIVADO_BIN="${VIVADO_BIN:-vivado}"
ABC_BIN="${ABC_BIN:-abc}"
PART="${PART:-xc7a200tfbg484-1}"
PERIOD_NS="${PERIOD_NS:-5.0}"
K_LUT="${K_LUT:-6}"

OUT_PARENT="${OUT_PARENT:-results/paper_full_eval}"

ACCURACY_PS="${ACCURACY_PS:-0.045,0.10,0.25,0.45,0.50,0.55,0.75,0.90}"
ACCURACY_TRIALS="${ACCURACY_TRIALS:-4000}"
ACCURACY_BATCH_SIZE="${ACCURACY_BATCH_SIZE:-2000}"
ACCURACY_SEED="${ACCURACY_SEED:-20260221}"

echo "[INFO] Repo root:  $REPO_ROOT"
echo "[INFO] n range:    $N_MIN..$N_MAX (odd only)"
echo "[INFO] methods:    $METHODS"
echo "[INFO] reprs:      $LOGIC_REPRS"
echo "[INFO] flow modes: $FLOW_MODES"
echo "[INFO] vivado:     $RUN_VIVADO"
echo "[INFO] accuracy:   $RUN_ACCURACY"
echo "[INFO] out parent: $OUT_PARENT"

cd "$REPO_ROOT"

cmake -S "$REPO_ROOT" -B "$REPO_ROOT/build"
cmake --build "$REPO_ROOT/build" -j

CMD=(
  python3 "$REPO_ROOT/scripts/run_cancel_tree_eval.py"
  --n-min "$N_MIN"
  --n-max "$N_MAX"
  --methods "$METHODS"
  --logic-reprs "$LOGIC_REPRS"
  --flow-modes "$FLOW_MODES"
  --out-dir "$OUT_PARENT"
)

if [[ "$RUN_VIVADO" == "1" ]]; then
  CMD+=(
    --run-vivado
    --full-sweep-vivado
    --vivado-bin "$VIVADO_BIN"
    --part "$PART"
    --period-ns "$PERIOD_NS"
  )
fi

if [[ "$RUN_ACCURACY" == "1" ]]; then
  CMD+=(
    --run-accuracy
    --accuracy-ps "$ACCURACY_PS"
    --accuracy-trials "$ACCURACY_TRIALS"
    --accuracy-batch-size "$ACCURACY_BATCH_SIZE"
    --accuracy-seed "$ACCURACY_SEED"
  )
fi

echo "[INFO] Running end-to-end evaluation ..."
"${CMD[@]}"

RUN_DIR="$(ls -dt "$REPO_ROOT/$OUT_PARENT"/run_* | head -n1)"
if [[ -z "$RUN_DIR" ]]; then
  echo "[ERROR] Could not locate generated run directory under $REPO_ROOT/$OUT_PARENT" >&2
  exit 1
fi

echo "[INFO] Run directory: $RUN_DIR"

MERGED_CSV="$RUN_DIR/final_with_abc.csv"
python3 "$REPO_ROOT/scripts/merge_eval_with_abc.py" \
  --run-dir "$RUN_DIR" \
  --out-csv "$MERGED_CSV" \
  --abc-bin "$ABC_BIN" \
  --k-lut "$K_LUT"

export PAPER_RUN_DIR="$RUN_DIR"
export PAPER_MERGED_CSV="$MERGED_CSV"
python3 - <<'PY'
import csv
import os
from pathlib import Path

run_dir = Path(os.environ["PAPER_RUN_DIR"]).resolve()
merged_csv = Path(os.environ["PAPER_MERGED_CSV"]).resolve()
post_route_csv = run_dir / "paper_post_route_rows.csv"
preferred_csv = run_dir / "paper_preferred_rows.csv"
all_flows_csv = run_dir / "paper_all_flows.csv"

with merged_csv.open("r", newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

if not rows:
    raise SystemExit("final_with_abc.csv is empty")

fieldnames = list(rows[0].keys())

post_route_rows = [r for r in rows if r.get("status") == "ok" and r.get("flow_mode") == "post_route"]
with post_route_csv.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(post_route_rows)

with preferred_csv.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

with (run_dir / "master_results.csv").open("r", newline="", encoding="utf-8") as f:
    all_flow_rows = list(csv.DictReader(f))

if all_flow_rows:
    all_fieldnames = list(all_flow_rows[0].keys())
    with all_flows_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=all_fieldnames)
        writer.writeheader()
        writer.writerows(all_flow_rows)
PY

echo "[INFO] Done."
echo "[INFO] Raw physical CSV:      $RUN_DIR/master_results.csv"
echo "[INFO] ABC merged CSV:        $MERGED_CSV"
echo "[INFO] All-flow CSV:          $RUN_DIR/paper_all_flows.csv"
echo "[INFO] Post-route-only CSV:   $RUN_DIR/paper_post_route_rows.csv"
echo "[INFO] Preferred rows CSV:    $RUN_DIR/paper_preferred_rows.csv"
echo "[INFO] Accuracy CSV:          $RUN_DIR/accuracy_results.csv"
echo "[INFO] Report:                $RUN_DIR/report_cancel_tree.md"
echo "[INFO] Figures:               $RUN_DIR/figures"
