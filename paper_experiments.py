"""Matched Beta-reference experiments for hat-, check-, and bar-Gamma in main.tex."""

import argparse
import json
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

from car_data import generate_car_residuals
from experiment_utils import stream, summarize
from methods import (
    SUPPORT_BOUND,
    Disk,
    DiskIntersection,
    fit_cascade_predictor,
    fit_intersection_predictor,
    fit_threshold_predictor,
)
from paper_plots import paired_summary, plot_distributions, plot_grouped_sets, report_tables


def replot(directory):
    """Regenerate figures from saved results without simulation or fitting."""
    config = json.loads((directory / "config.json").read_text(encoding="utf-8"))
    frame = pd.read_csv(directory / "runs.csv")
    first_runs = []
    for horizon in config["horizons"]:
        with np.load(directory / f"first_run_T{horizon}.npz") as saved:
            arrays = {key: saved[key].copy() for key in saved.files}
        models, removed = {}, {}
        for key in arrays:
            if not key.endswith("_centers"):
                continue
            name = key.removesuffix("_centers")
            disks = [
                Disk(center, float(radius))
                for center, radius in zip(arrays[key], arrays[name + "_radii_squared"])
            ]
            models[name] = DiskIntersection(disks) if name.startswith("bar_") else disks[0]
            removed[name] = arrays[name + "_removed"]
        first_runs.append((horizon, arrays["calibration"], arrays["test"], models, removed))
    plot_grouped_sets(directory, first_runs)
    plot_distributions(directory, frame, config["n_cal"], config["rank"], config["n_test"])
    print(f"Regenerated plots: {directory.resolve()}")
    return directory


def schedules(rank):
    """The two fixed intersection schedules specified in the numerical section."""
    if rank % 2 or rank // 2 < SUPPORT_BOUND:
        raise ValueError("The two equal batches require even rank >= 2*d")
    return [(SUPPORT_BOUND,) * (rank // SUPPORT_BOUND), (rank // 2, rank // 2)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rank", type=int, default=18, help="Target Beta(N-r+1,r)")
    parser.add_argument("--n-cal", type=int, default=200)
    parser.add_argument("--n-test", type=int, default=5000)
    parser.add_argument("--runs", type=int, default=100)
    parser.add_argument("--horizons", nargs="+", type=int, default=[2, 30, 50])
    parser.add_argument("--seed", type=int, default=20261019)
    parser.add_argument("--quick", action="store_true", help="Two runs with the paper sample sizes")
    parser.add_argument("--outdir", type=Path, default=Path("results/default"))
    parser.add_argument(
        "--replot", action="store_true", help="Only redraw saved results in --outdir"
    )
    args = parser.parse_args(argv)
    if args.replot:
        return replot(args.outdir)
    if args.quick:
        args.runs = min(args.runs, 2)
    n, r = args.n_cal, args.rank
    if not SUPPORT_BOUND <= r < n or r % SUPPORT_BOUND:
        raise ValueError("Matching Romao requires d <= r < N and r divisible by d=3")
    if min(args.n_test, args.runs) < 1 or args.seed < 0:
        raise ValueError("Positive runs/test size and nonnegative seed required")
    if len(set(args.horizons)) != len(args.horizons) or any(h < 2 or h > 50 for h in args.horizons):
        raise ValueError("Use distinct horizons between 2 and 50")
    batches = schedules(r)
    directory = args.outdir
    directory.mkdir(parents=True, exist_ok=True)
    config = dict(
        vars(args),
        outdir=str(directory),
        support_bound=SUPPORT_BOUND,
        beta_parameters=[n - r + 1, r],
        alpha=r / (n + 1),
        alpha_intersection=float(np.nextafter(r / (n + 1), 1.0)),
        schedules=batches,
        area="deterministic vertical-slice quadrature",
        reference="hat: coverage CDF upper bound; check/bar: exact under paper assumptions",
    )
    (directory / "config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    rows = []
    first_runs = []
    with (directory / "stages.jsonl").open("w", encoding="utf-8") as audit:
        for run in range(args.runs):
            calibration = generate_car_residuals(n, stream(args.seed, run, 0))
            testing = generate_car_residuals(args.n_test, stream(args.seed, run, 1))
            for horizon in args.horizons:
                cal, test = calibration[:, horizon - 1], testing[:, horizon - 1]
                variants = [("hat", None), ("check", None)] + [
                    (f"bar_L{len(s)}", s) for s in batches
                ]
                models, removed_indices = {}, {}
                for name, schedule in variants:
                    row = dict(run=run, horizon=horizon, method=name, valid=False, error="")
                    try:
                        start = perf_counter()
                        if name == "hat":
                            result = fit_cascade_predictor(cal, r // SUPPORT_BOUND - 1)
                            model, stages = result.disk, result.stages
                            fits, removed = r // SUPPORT_BOUND, r - SUPPORT_BOUND
                        elif name == "check":
                            model, _ = fit_threshold_predictor(cal, r)
                            stages, fits, removed = [], 1, 0
                        else:
                            model, stages = fit_intersection_predictor(cal, schedule)
                            fits, removed = len(schedule), sum(schedule[:-1])
                        fit_seconds = perf_counter() - start
                        if run == 0:
                            models[name] = model
                            removed_indices[name] = (
                                result.removed.tolist()
                                if name == "hat"
                                else [i for stage in stages[:-1] for i in stage["removed"]]
                            )
                        start = perf_counter()
                        hits = int(np.count_nonzero(model.contains(test)))
                        query_seconds = perf_counter() - start
                        start = perf_counter()
                        area = model.area()
                        row.update(
                            valid=True,
                            coverage=hits / args.n_test,
                            hits=hits,
                            n_test=args.n_test,
                            area=area,
                            fits=fits,
                            removed=removed,
                            fit_seconds=fit_seconds,
                            query_seconds=query_seconds,
                            area_seconds=perf_counter() - start,
                        )
                        audit.write(
                            json.dumps(dict(run=run, horizon=horizon, method=name, stages=stages))
                            + "\n"
                        )
                    except (
                        ValueError,
                        RuntimeError,
                        FloatingPointError,
                        np.linalg.LinAlgError,
                    ) as exc:
                        row["error"] = str(exc)
                    rows.append(row)
                if run == 0:
                    first_runs.append((horizon, cal, test, models, removed_indices))
            pd.DataFrame(rows).to_csv(directory / "runs.csv", index=False)
            print(f"Completed {run + 1}/{args.runs}", flush=True)
    plot_grouped_sets(directory, first_runs)
    frame = pd.DataFrame(rows)
    summary = summarize(frame)
    summary.to_csv(directory / "summary.csv", index=False)
    paired = paired_summary(frame)
    paired.to_csv(directory / "paired.csv", index=False)
    plot_distributions(directory, frame, n, r, args.n_test)
    report = [
        "# Matched paper experiments",
        "",
        f"Target coverage reference: Beta({n - r + 1}, {r}); alpha={r / (n + 1):.8g}.",
        "hat: Romao final disk, ell=r/d. check: one-fit residual threshold.",
        "bar_L6: batches (3,3,3,3,3,3); bar_L2: batches (9,9) at the default rank 18.",
        "All schedules are fixed before observing data. Intersection alpha is the next float above r/(N+1).",
        "The hat coverage CDF is bounded above by the Beta CDF; equality is not assumed.",
        "The check/bar laws are exact subject to the assumptions in main.tex.",
        "Disk problems can have fewer than d supports: bar uses residual ranking, hat uses label padding.",
        "Fit times include standalone solves; area quadrature is timed separately.",
        "The plotted reference is the theoretical true-coverage Beta CDF. Empirical coverage includes finite-validation noise.",
        "See summary.csv for run-level 95% mean intervals and paired.csv for all pairwise differences.",
        f"Failures: {int((~frame.valid).sum())}; inspect runs.csv. Successful-only summaries are conditional.",
    ]
    report.append(report_tables(summary, paired))
    if args.runs < 30:
        report.append("Smoke check only: too few runs for distributional conclusions.")
    (directory / "report.md").write_text("\n\n".join(report) + "\n", encoding="utf-8")
    print(
        summary[["horizon", "method", "coverage", "area", "fit_seconds", "failures"]].to_string(
            index=False
        )
    )
    print(f"Results: {directory.resolve()}")
    if not frame.valid.all():
        raise RuntimeError("Experiment incomplete; inspect runs.csv")
    return directory


if __name__ == "__main__":
    main()
