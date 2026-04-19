#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

N_MIN="${1:-5}"
N_MAX="${2:-61}"

METHODS="${METHODS:-popcount_strict,approx_block3}"
LOGIC_REPRS="${LOGIC_REPRS:-xag,aig,mig}"
FLOW_MODES="${FLOW_MODES:-synth_only}"

RUN_VIVADO="${RUN_VIVADO:-1}"
VIVADO_BIN="${VIVADO_BIN:-vivado}"
PART="${PART:-xc7a200tfbg484-1}"
PERIOD_NS="${PERIOD_NS:-5.0}"

ABC_BIN="${ABC_BIN:-abc}"
K_LUT="${K_LUT:-6}"

OUT_PARENT="${OUT_PARENT:-results/strict_vs_approx_eval}"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_ROOT="$REPO_ROOT/$OUT_PARENT/run_$STAMP"
mkdir -p "$RUN_ROOT"

echo "[INFO] Repo root: $REPO_ROOT"
echo "[INFO] Run root:  $RUN_ROOT"
echo "[INFO] n range:   $N_MIN..$N_MAX (odd only)"
echo "[INFO] methods:   $METHODS"
echo "[INFO] reprs:     $LOGIC_REPRS"
echo "[INFO] flow:      $FLOW_MODES"
echo "[INFO] vivado:    $RUN_VIVADO"

cmake -S "$REPO_ROOT" -B "$REPO_ROOT/build"
cmake --build "$REPO_ROOT/build" -j

IFS=',' read -r -a REPRS <<< "$LOGIC_REPRS"
RUN_DIRS=()

for REPR_RAW in "${REPRS[@]}"; do
  REPR="$(echo "$REPR_RAW" | tr '[:upper:]' '[:lower:]' | xargs)"
  if [[ -z "$REPR" ]]; then
    continue
  fi
  REPR_DIR="$(echo "$REPR" | tr '[:lower:]' '[:upper:]')"
  OUT_DIR_REL="$OUT_PARENT/run_$STAMP/$REPR_DIR"

  CMD=(
    "$REPO_ROOT/scripts/run_cancel_tree_eval.py"
    --n-min "$N_MIN"
    --n-max "$N_MAX"
    --methods "$METHODS"
    --logic-reprs "$REPR"
    --flow-modes "$FLOW_MODES"
    --out-dir "$OUT_DIR_REL"
    --skip-build
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

  echo "[INFO] Running repr=$REPR ..."
  "${CMD[@]}"

  LATEST_RUN="$(ls -dt "$REPO_ROOT/$OUT_DIR_REL"/run_* | head -n1)"
  RUN_DIRS+=("$LATEST_RUN")
  echo "[INFO] repr=$REPR output: $LATEST_RUN"
done

printf "%s\n" "${RUN_DIRS[@]}" > "$RUN_ROOT/run_dirs.txt"

MERGED_CSV="$RUN_ROOT/final_stats_n${N_MIN}_to_${N_MAX}.csv"
MERGE_CMD=(
  python3 "$REPO_ROOT/scripts/merge_eval_with_abc.py"
  --out-csv "$MERGED_CSV"
  --abc-bin "$ABC_BIN"
  --k-lut "$K_LUT"
)
for D in "${RUN_DIRS[@]}"; do
  MERGE_CMD+=(--run-dir "$D")
done
"${MERGE_CMD[@]}"

echo "[INFO] Done."
echo "[INFO] Final merged CSV: $MERGED_CSV"
echo "[INFO] Run dirs list:     $RUN_ROOT/run_dirs.txt"
