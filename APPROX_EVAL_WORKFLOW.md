# Approximate Majority Evaluation Workflow

## Objective
Evaluate whether approximate majority can reduce hardware area to roughly half of exact `popcount+threshold` while keeping usable output accuracy.

## Implemented Approximation
Method name: `approx_block3` (alias `approx_block3_popcount`)

Exact baseline alternative now available:
- `popcount_strict` (alias `baseline_strict`): scaffolded paper-style baseline (`N=2^p-1`) using CSA reduction plus strict ripple threshold carry decision.

Process for odd `n` inputs:
1. Partition input bits into full groups of 3, plus remainder bits.
2. Replace each full 3-bit group with `maj3(a,b,c)`.
3. Pass leftover 1-2 bits unchanged.
4. Run exact majority (`popcount+threshold`) on the compressed vector.

For odd `n`, compressed size is still odd:
- `n = 3g + r`, `r in {0,1,2}`
- compressed size = `g + r` (odd when `n` is odd)

## Tool Flow
1. Generate networks with `majgen`.
2. Sweep logic representations (`xag`, `aig`, `mig`) and collect stats per representation.
3. Run Vivado OOC (`synth_only` and `post_route`) on subset or full sweep.
4. Compute area ratio vs baseline:
   - prefer `lut_count` when available
   - fallback by representation otherwise:
     - `xag`: `xag_and_count + xag_xor_count`
     - `aig`: `aig_and_count`
     - `mig`: `mig_maj_count`
5. Mark half-resource success when `ratio <= 0.5`.
6. Run Monte Carlo accuracy sweeps across `p` values.

## Accuracy Sweep Model
For each method, odd `n`, and input Bernoulli probability `p`:
1. Sample random vectors (`trials` configurable).
2. Compute exact majority label.
3. Compute method output:
   - exact methods: always exact
   - `approx_block3`: compressed pipeline model above
4. Report accuracy and error rate.

Recommended p set for stress + biased regimes:
- `0.10, 0.25, 0.50, 0.55, 0.75, 0.90`

## Primary Command
```bash
./scripts/run_cancel_tree_eval.py \
  --n-min 5 --n-max 511 \
  --methods popcount,popcount_strict,approx_block3 \
  --logic-reprs xag,aig,mig \
  --run-vivado \
  --vivado-subset 63,127,255,511 \
  --run-accuracy \
  --accuracy-ps 0.10,0.25,0.50,0.55,0.75,0.90 \
  --accuracy-trials 4000
```

## Artifacts to Use in Analysis
- `master_results.csv`
- `tables_half_resource.csv` and `.tex`
- `accuracy_results.csv`
- `tables_accuracy.csv` and `.tex`
- `figures/area_ratio_vs_baseline.*`
- `figures/accuracy_vs_n_pXX.*`
- `figures/accuracy_vs_p.*`
- `report_cancel_tree.md`
