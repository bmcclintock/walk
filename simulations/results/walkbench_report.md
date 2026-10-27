# walk devel vs modified (speed-up): benchmark report

Generated 2026-10-27 08:15 by `bench_walk_devel_vs_modified.R`. Run time 10.9 min.

## Setup

- **devel:** https://github.com/dsjohnson/walk.git @ `devel` (b39f00b).
- **modified:** https://github.com/bmcclintock/walk.git @ speed-up (55e35ff).
- **Versions loaded:** walk 0.0.1.9005 (devel) and 0.0.1.9005 (modified; OpenMP TRUE). Each ran in its own R process.
- **Simulated data:** langevinSSM 0.0.2, `simLangevin` (underdamped, beta = (-4, 6, 5, -0.1), sigma = 5, gamma = 0.5) on `exampleCovs`.
  - 4 animals, 1200 fixes, subsampled to irregular times (median dt 0.070).
  - Argos-type error ellipses (median smaj 2.38).
- **Grids:** `exampleCovs` aggregated by 4, 2, 1, i.e. 10,000, 40,000, 160,000 cells.
- **Model:** movement `~cov1 + cov2 + cov3`, log links, `norm = FALSE`. −2ll evaluated at 3 parameter vectors.
- **Timing:** median of 3 run(s) after one warm-up. Machine: Darwin, 12 cores, R 4.5.2.

## Checks

- **L, devel vs modified `buffer` (abs tol 1e-10): PASS.**
  - grid4: 15097 vs 15097 non-zeros, 0/0 cells only in one version, max abs diff 4.44e-16.
  - grid2: 42657 vs 42657 non-zeros, 0/0 cells only in one version, max abs diff 9.16e-16.
  - grid1: 138651 vs 138651 non-zeros, 0/0 cells only in one version, max abs diff 3.15e-15.
- **−2ll, devel vs modified (`buffer` L, and the modified likelihood on devel's own L), rel tol 1e-08: PASS.**
  - Max relative difference 4.22e-16 over 36 comparisons.
- **Joint call vs sum of per-animal values (modified):** max abs diff 0.
- **Full fit of animal 1** (10,000 cells, nlminb, same start):
  - −2ll: devel 897.3738143, modified 897.3738143 (difference -3.41e-13).
  - Max abs parameter difference 2.07e-08; convergence codes 0/0.
- **`disc` vs `buffer` (not expected to be identical):**
  - Extra cells in `disc` that the polygonised buffer misses: 10,000 cells: 15; 40,000 cells: 68; 160,000 cells: 277.
  - Its −2ll differs from devel by 1.7e-07 to 1.1e-06 (relative).

## Speed

Speed-up = devel time / modified time. All times are for all animals together.

### `proc_Langevin`

| cells | devel | buffer, per animal, 1 thr | buffer, one call, 2 thr | disc, per animal, 1 thr | disc, one call, 2 thr |
|---|---|---|---|---|---|
| 10,000 | 2.34 s | 0.132 s (18×) | 0.103 s (23×) | 0.018 s (130×) | 0.007 s (334×) |
| 40,000 | 7.21 s | 0.225 s (32×) | 0.194 s (37×) | 0.028 s (257×) | 0.013 s (554×) |
| 160,000 | 25.67 s | 0.601 s (43×) | 0.522 s (49×) | 0.063 s (407×) | 0.032 s (802×) |

### −2 log-likelihood (one evaluation; range over parameter vectors)

| cells | devel | per animal, devel's L | per animal, buffer | joint, 1 thr | joint, 2 thr |
|---|---|---|---|---|---|
| 10,000 | 0.65–0.98 s | 0.029–0.048 s (20–22×) | 0.030–0.047 s (21–23×) | 0.010–0.026 s (38–65×) | 0.009–0.018 s (55–72×) |
| 40,000 | 2.41–4.00 s | 0.124–0.283 s (14–20×) | 0.130–0.282 s (14–20×) | 0.037–0.195 s (21–65×) | 0.034–0.120 s (33–72×) |
| 160,000 | 8.97–18.85 s | 0.532–3.101 s (6–18×) | 0.571–3.059 s (6–18×) | 0.169–2.647 s (7–53×) | 0.155–1.620 s (12–61×) |

**Full fit, animal 1:** devel 15.2 s → modified 0.71 s (21×).

## Notes

- **Devel has no joint likelihood,** so the devel figure is `ctmc_n2ll` per animal, summed.
- **"per animal, devel's L"** runs the modified likelihood on exactly the inputs devel used. It isolates the likelihood code from the `proc_Langevin` step.
- **`disc`** uses every cell that touches the exact error circle, while `buffer` uses `terra::buffer` + `terra::cells(touches = TRUE)`. Only `buffer` is expected to match devel.
- **Files:** `walkbench_summary.txt` (all tables), `walkbench_*.csv`, `walkbench_results.rds`.
