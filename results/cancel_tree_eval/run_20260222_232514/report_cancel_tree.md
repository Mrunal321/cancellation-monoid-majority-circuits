# Majority Generator Evaluation

## Method
This run compares exact and/or approximate odd-`n` majority generators using a shared flow (Mockturtle network generation, optional Vivado OOC synth/post-route, and uniform CSV/plot/reporting outputs).

Approximation used here: `approx_block3[_popcount]`, which applies one compression stage of 3-input majority (`maj3`) on disjoint triples, forwards tail bits unchanged, and then applies exact majority on the compressed vector.

## Process
- Sweep range: odd n values in requested interval.
- Methods: `approx_block3`, `popcount_strict`
- Logic representations: `aig`, `mig`, `xag`
- Baseline method for ratios/tables: `popcount_strict`
- Area ratio target tracked: method_area / baseline_area <= 0.5

## Reproducibility
- Date: 2026-02-22T23:25:47.401512
- Repo: /home/mrunal/Desktop/Cancellation-Monoid Tree-Project
- Git commit: N/A (no git repo)
- FPGA part: `xc7a200tfbg484-1`
- Timing constraint: `set_max_delay 5.0 -from [all_inputs] -to [all_outputs]`
- majgen: majgen 0.1.0
- vivado: vivado v2024.2 (64-bit)
- yosys: Yosys 0.33 (git sha1 2584903a060)
- abc: UC Berkeley, ABC 1.01 (compiled Nov 27 2025 00:21:43)

## Outputs
- `master_results.csv` (generation + timing/area/power rows)
- `figures/lut_vs_n.(png|pdf)`
- `figures/delay_vs_n.(png|pdf)`
- `figures/power_vs_n.(png|pdf)`
- `figures/depth_vs_n.(png|pdf)`
- `figures/area_delay_scatter.(png|pdf)`
- `figures/area_ratio_vs_baseline.(png|pdf)`
- `tables_cancel_tree.csv`, `tables_cancel_tree.tex`
- `tables_half_resource.csv`, `tables_half_resource.tex`
- `accuracy_results.csv` (if `--run-accuracy`)
- `tables_accuracy.csv`, `tables_accuracy.tex` (if `--run-accuracy`)
- `figures/accuracy_vs_n_pXX.(png|pdf)` and `figures/accuracy_vs_p.(png|pdf)` (if `--run-accuracy`)

## Observations
- `approx_block3[xag]` half-resource hit rate (area ratio <= 0.5): 30.31% over 254 n-points.
- `approx_block3[aig]` half-resource hit rate (area ratio <= 0.5): 69.69% over 254 n-points.
- `approx_block3[mig]` half-resource hit rate (area ratio <= 0.5): 42.52% over 254 n-points.
- Accuracy at p=0.50 for `approx_block3`: mean=0.8372, min=0.8163, max=0.9413.
- Accuracy at p=0.50 for `popcount_strict`: mean=1.0000, min=1.0000, max=1.0000.

## TODO
- Carry-chain friendly compare/sub implementations.
- Optional multi-stage approximation space exploration.
- Optional pipelining for Fmax-focused experiments.

## Cancel-Tree Background
The `cancel_tree` method builds a balanced merge tree of Boyer-Moore style summaries `(cand, k)`, where `cand` is a 1-bit candidate and `k` is the surplus count after pair cancellation. Each merge follows the exact monoid combine rules with empty-summary priority and deterministic zero-candidate handling for `k=0`. Widths are assigned as `w(sz)=ceil(log2(sz+1))` and both children are zero-extended to the merge width before add/sub/compare.

