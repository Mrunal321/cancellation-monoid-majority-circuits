# Cancellation-Monoid Tree Project

This project implements odd-`n` majority generators in Mockturtle:
- exact methods (`popcount`, `popcount_strict`, `cancel_tree`, `cancel_tree_v2`)
- an approximate method (`approx_block3`) that compresses input triples with `maj3` and then runs exact majority on the compressed vector

The repository keeps the source tree plus a curated subset of paper-relevant CSVs/netlists.
Bulk generated outputs such as full Vivado run directories, build products, and local logs are intentionally excluded from git.

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

## Unit/Correctness Tests

```bash
ctest --test-dir build --output-on-failure
```

Or directly:

```bash
./build/test_cancel_tree
```

The test executable performs:
- Exhaustive checks for `n = 3,5,7,9,11`
- Random checks (10k vectors each) for `n = 31,63,127`
- Equivalence checks vs baseline exact `popcount`
- Debug summary dumps for `n=7` and `n=11`
- DOT dump `artifacts/maj_7_cancel_tree.dot`

## Single Network Generation

```bash
./build/majgen --maj_impl cancel_tree --logic_repr xag --n 63 --out out/maj_63_cancel_tree_xag.v --format verilog --stats_json out/maj_63_cancel_tree_xag.json
```

Supported methods:
- `cancel_tree` (alias: `boyermoore_tree`)
- `cancel_tree_v2` (alias: `boyermoore_tree_v2`)
- `popcount`
- `popcount_strict` (alias: `baseline_strict`; scaffolded baseline with CSA + threshold ripple carry)
- `approx_block3` (alias: `approx_block3_popcount`)

Supported output formats:
- `verilog`
- `blif`
- `aig`
- `dot`

Supported internal logic representations:
- `xag`
- `aig`
- `mig`

`majgen` now accepts:
- `--logic_repr <xag|aig|mig>` to pick the representation used for emitted netlist and `stats_json`.

## Full Sweep / Evaluation

### Generation + stats only (full odd sweep)

```bash
./scripts/run_cancel_tree_eval.py --n-min 5 --n-max 511 --methods popcount,approx_block3 --logic-reprs xag
```

### Representation sweep (XAG + AIG + MIG)

```bash
./scripts/run_cancel_tree_eval.py \
  --n-min 5 --n-max 511 \
  --methods popcount,approx_block3 \
  --logic-reprs xag,aig,mig
```

### Paper-style strict baseline comparison

```bash
./scripts/run_cancel_tree_eval.py \
  --n-min 5 --n-max 511 \
  --methods popcount,popcount_strict,approx_block3 \
  --logic-reprs xag,aig,mig
```

### One-shot strict-vs-approx sweep (n=5..61)

Generates separate `XAG/`, `AIG/`, `MIG/` run folders and a single merged CSV that includes:
- original network/device stats from `master_results.csv`
- ABC optimized stats (`dc2`)
- ABC mapped stats (`if -K 6`)

```bash
./scripts/run_strict_vs_approx_sweep.sh 5 61
```

Useful env overrides:
- `RUN_VIVADO=1` (default) or `RUN_VIVADO=0`
- `FLOW_MODES=synth_only` (default) or `FLOW_MODES=synth_only,post_route`
- `LOGIC_REPRS=xag,aig,mig`
- `METHODS=popcount_strict,approx_block3`
- `OUT_PARENT=results/strict_vs_approx_eval`

### With Vivado (representative subset)

```bash
./scripts/run_cancel_tree_eval.py \
  --n-min 5 --n-max 511 \
  --methods popcount,approx_block3 \
  --logic-reprs xag,aig,mig \
  --run-vivado \
  --vivado-subset 63,127,255,511
```

### With Vivado (full sweep)

```bash
./scripts/run_cancel_tree_eval.py \
  --n-min 5 --n-max 511 \
  --methods popcount,approx_block3 \
  --logic-reprs xag,aig,mig \
  --run-vivado \
  --full-sweep-vivado
```

### Accuracy sweep across input bias p

```bash
./scripts/run_cancel_tree_eval.py \
  --n-min 5 --n-max 511 \
  --methods popcount,approx_block3 \
  --logic-reprs xag,aig,mig \
  --run-accuracy \
  --accuracy-ps 0.10,0.25,0.50,0.55,0.75,0.90 \
  --accuracy-trials 4000
```

Outputs are written under:

- `results/cancel_tree_eval/run_<timestamp>/master_results.csv`
- `results/cancel_tree_eval/run_<timestamp>/figures/*.png,*.pdf`
- `results/cancel_tree_eval/run_<timestamp>/tables_cancel_tree.csv`
- `results/cancel_tree_eval/run_<timestamp>/tables_cancel_tree.tex`
- `results/cancel_tree_eval/run_<timestamp>/tables_half_resource.csv`
- `results/cancel_tree_eval/run_<timestamp>/tables_half_resource.tex`
- `results/cancel_tree_eval/run_<timestamp>/accuracy_results.csv` (when `--run-accuracy`)
- `results/cancel_tree_eval/run_<timestamp>/tables_accuracy.csv` (when `--run-accuracy`)
- `results/cancel_tree_eval/run_<timestamp>/tables_accuracy.tex` (when `--run-accuracy`)
- `results/cancel_tree_eval/run_<timestamp>/report_cancel_tree.md` (+ PDF when pandoc is available)

Vivado part defaults to `xc7a200tfbg484-1` and flow modes are `synth_only` and `post_route`.
