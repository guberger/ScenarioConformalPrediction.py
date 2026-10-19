"""Random streams and run-level summaries shared by the paper experiment."""

import math

import numpy as np
import pandas as pd
from scipy.stats import binom, t

from methods import SUPPORT_BOUND


def stream(seed, run, purpose):
    return np.random.default_rng(np.random.SeedSequence([seed, run, purpose]))


def mean_interval(values):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return (float("nan"),) * 3
    mean = float(values.mean())
    half = (
        float(t.ppf(0.975, len(values) - 1) * values.std(ddof=1) / np.sqrt(len(values)))
        if len(values) > 1
        else float("nan")
    )
    return mean, mean - half, mean + half


def summarize(frame):
    rows = []
    for (horizon, method), group in frame.groupby(["horizon", "method"], sort=True):
        good = group[group.valid]
        row = dict(
            horizon=horizon,
            method=method,
            runs=len(group),
            successes=len(good),
            failures=len(group) - len(good),
            status="complete" if group.valid.all() else "incomplete",
        )
        for metric in ("coverage", "area", "fit_seconds", "query_seconds"):
            values = good[metric] if metric in good else []
            mean, low, high = mean_interval(values)
            row.update({metric: mean, metric + "_ci_low": low, metric + "_ci_high": high})
        row["estimated_coverage_q05"] = (
            float(good.coverage.quantile(0.05)) if len(good) else float("nan")
        )
        rows.append(row)
    return pd.DataFrame(rows)


def select_matched_rank(n, eps=0.1, delta=None, rounds=None):
    """Return (rank, removal rounds) matching Beta(n-rank+1, rank).

    eps alone targets marginal miscoverage. With delta, select the largest
    admissible multiple of d whose binomial tail is <= delta.
    """
    if not isinstance(n, (int, np.integer)) or n <= SUPPORT_BOUND or not 0 < eps < 1:
        raise ValueError("Require integer n > d and 0 < eps < 1")
    if delta is not None and not 0 < delta < 1:
        raise ValueError("Require 0 < delta < 1")
    if rounds is not None:
        if delta is not None:
            raise ValueError("Choose either explicit rounds or confidence calibration")
        if not isinstance(rounds, (int, np.integer)) or rounds < 0:
            raise ValueError("rounds must be a nonnegative integer")
        rank = SUPPORT_BOUND * (rounds + 1)
    elif delta is None:
        rank = SUPPORT_BOUND * math.floor(eps * (n + 1) / SUPPORT_BOUND)
    else:
        ranks = np.arange(SUPPORT_BOUND, n, SUPPORT_BOUND)
        eligible = ranks[binom.cdf(ranks - 1, n, eps) <= delta]
        rank = int(eligible[-1]) if len(eligible) else 0
    if not SUPPORT_BOUND <= rank < n:
        raise ValueError("No admissible matched rank; increase n or relax the target")
    return int(rank), int(rank // SUPPORT_BOUND - 1)
