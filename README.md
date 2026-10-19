# Scenario-based conformal prediction experiments

The experiments compare four prediction sets from `papers/main.tex` on simulated
car trajectories at horizons **T = 2, 30, 50**:

- **Cascade removal:** six fits, returning the final disk.
- **Residual thresholding:** one fit on all training points, followed by a rank-18 threshold.
- **Removal with intersection:** thresholded disks intersected across six batches
  of size 3, or two batches of size 9.

All methods target the same **Beta(183, 18)** coverage reference, using **200
training trajectories**, **5,000 independent test trajectories**, and **100 runs**.
The Beta law is exact for the proposed predictors under the paper's assumptions;
for cascade removal, it provides a stochastic coverage bound.

## Run

Install Python 3.14+ and [uv](https://docs.astral.sh/uv/), then run from this folder:

```powershell
uv sync --locked
uv run python paper_experiments.py
```

For a two-run smoke check:

```powershell
uv run python paper_experiments.py --quick
```

Change the repetition count or list all options:

```powershell
uv run python paper_experiments.py --runs 1000
uv run python paper_experiments.py --help
```

## Results

Regenerate plots from saved results without rerunning the experiments:

```powershell
uv run python paper_experiments.py --replot
# Or use another results folder:
uv run python paper_experiments.py --replot --outdir PATH
```

Results are written directly to `results/paper/`, replacing the previous outputs.
Use `--outdir PATH` to choose another folder. Outputs include:

- `sets.png`: first-run prediction sets and training points, with one panel per horizon.
- `coverage.png`: empirical coverage CDFs compared with the theoretical Beta CDF.
- `volume.png`: prediction-set area across runs.
- `runtime.png`: fitting time across runs.
- `report.md` and CSV files: summaries and paired comparisons between methods.

The experiments show the tradeoff between coverage, prediction-set size, and
computation. In the default 100-run study, cascade removal gives higher coverage
than the proposed methods, while the six-stage intersection gives the smallest
mean area at each horizon. One-fit residual thresholding is fastest. Empirical
coverage uses a finite test sample, so it includes validation noise and need not
follow the theoretical true-coverage CDF exactly.
