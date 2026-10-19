"""Plot and report helpers for the matched paper experiments."""

from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import beta

from experiment_utils import mean_interval
from methods import minimum_enclosing_disk

plt.rcParams.update(
    {
        "figure.constrained_layout.h_pad": 0.10,
        "font.size": 14,
        "axes.labelsize": 15,
        "axes.titlesize": 17,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "legend.fontsize": 17,
        "figure.titlesize": 19,
    }
)


def shared_legend(fig, axes):
    """Collect unique labeled artists across panels into one figure legend."""
    entries = {}
    for ax in axes:
        handles, labels = ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            entries.setdefault(label, handle)
        if ax.get_legend() is not None:
            ax.get_legend().remove()
    if entries:
        fig.legend(
            entries.values(),
            entries.keys(),
            loc="outside lower center",
            ncol=5,
        )


def method_label(name):
    """Paper notation for the method identifiers stored in result tables."""
    return {
        "hat": r"$\hat\Gamma_{6}$",
        "check": r"$\check\Gamma_{\alpha}$",
        "bar_L6": r"$\bar\Gamma_{\mathbf{s},\alpha},\;\mathbf{s}=(3,3,3,3,3,3)$",
        "bar_L2": r"$\bar\Gamma_{\mathbf{s},\alpha},\;\mathbf{s}=(9,9)$",
    }.get(name, name)


def plot_distributions(directory, frame, n, rank, n_test):
    """Figures X2--X4: empirical coverage CDFs, areas, and fit times."""
    groups = list(frame.groupby("horizon", sort=False))
    if not groups:
        return
    figures = {}
    for filename in ("coverage", "volume", "runtime"):
        fig, axes = plt.subplots(
            1, len(groups), figsize=(7 * len(groups), 6), squeeze=False, constrained_layout=True
        )
        figures[filename] = (fig, axes[0])
    colors = dict(zip(dict.fromkeys(frame.method), plt.get_cmap("tab10").colors))
    for column, (horizon, group) in enumerate(groups):
        good = group[group.valid]
        if good.empty:
            continue
        names = list(dict.fromkeys(group.method))
        names = [name for name in names if (good.method == name).any()]
        a, b = n - rank + 1, rank
        low = max(0.0, min(float(good.coverage.min()) - 0.02, float(beta.ppf(0.001, a, b))))
        fig, axes = figures["coverage"]
        ax = axes[column]
        for name in names:
            values = np.sort(good.loc[good.method == name, "coverage"].to_numpy())
            # Extend to the plot limits so the CDF starts at zero and ends at one.
            ax.step(
                np.r_[low, values, 1.0],
                np.r_[0.0, np.arange(1, len(values) + 1) / len(values), 1.0],
                where="post",
                linewidth=2,
                color=colors[name],
                label=method_label(name),
            )
        x = np.linspace(low, 1, 1000)
        ax.plot(
            x,
            beta.cdf(x, a, b),
            color="black",
            linestyle="--",
            label=r"Theoretical true-coverage CDF ($\mathrm{Beta}$)",
        )
        ax.set(
            xlabel=r"Empirical coverage $\widehat C$",
            ylabel=r"Cumulative probability $F(c)$",
            ylim=(0, 1),
            xlim=(low, 1),
            title=rf"$T={horizon}$",
        )
        status = " | INCOMPLETE: successful runs only" if not group.valid.all() else ""
        ax.set_title(rf"$T={horizon}$" + status)
        n_run = group["run"].nunique()
        fig.suptitle(
            rf"Coverage CDFs | $N={n},\;N_{{\mathrm{{test}}}}={n_test},"
            rf"\;n_{{\mathrm{{run}}}}={n_run}$"
        )
        for metric, filename, ylabel, title in (
            (
                "area",
                "volume",
                r"Prediction-set area",
                "Prediction-set volume",
            ),
            (
                "fit_seconds",
                "runtime",
                r"Standalone fit time $t\;[\mathrm{s}]$",
                "Computation time",
            ),
        ):
            fig, axes = figures[filename]
            ax = axes[column]
            artists = ax.boxplot(
                [good.loc[good.method == name, metric] for name in names],
                patch_artist=True,
                showfliers=True,
            )
            for patch, name in zip(artists["boxes"], names):
                patch.set_facecolor(colors[name])
                patch.set_alpha(0.5)
            ax.set(ylabel=ylabel, title=rf"$T={horizon}$" + status)
            ax.set_xticks([])
            fig.suptitle(rf"{title} | $N={n},\;n_{{\mathrm{{run}}}}={n_run}$")
    for filename, (fig, axes) in figures.items():
        for ax in axes[1:]:
            ax.set_ylabel("")
        if filename == "coverage":
            shared_legend(fig, axes)
        else:
            handles = [
                Patch(facecolor=color, alpha=0.5, label=method_label(name))
                for name, color in colors.items()
            ]
            fig.legend(handles=handles, loc="outside lower center", ncol=5)
        fig.savefig(directory / f"{filename}.png", dpi=180)
        plt.close(fig)


def paired_summary(frame):
    rows = []
    for horizon, group in frame.groupby("horizon"):
        for metric in ("coverage", "area", "fit_seconds", "query_seconds"):
            if metric not in group:
                continue
            pivot = group[group.valid].pivot(index="run", columns="method", values=metric)
            for reference, method in combinations(pivot.columns, 2):
                values = (pivot[method] - pivot[reference]).dropna()
                mean, low, high = mean_interval(values)
                rows.append(
                    dict(
                        horizon=horizon,
                        metric=metric,
                        method=method,
                        reference=reference,
                        pairs=len(values),
                        mean=mean,
                        ci_low=low,
                        ci_high=high,
                    )
                )
    return pd.DataFrame(
        rows,
        columns=["horizon", "metric", "method", "reference", "pairs", "mean", "ci_low", "ci_high"],
    )


def plot_sets(directory, horizon, cal, test, models, removed, ax=None):
    """Show actual intersection boundary arcs, not the union of component disks."""
    base = minimum_enclosing_disk(cal)
    count = len(models)
    if not count:
        return
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    ax.scatter(cal[:, 0], cal[:, 1], s=14, c="#334155", alpha=0.6, label="training")
    colors = dict(zip(models, plt.get_cmap("tab10").colors))
    arrays = dict(
        calibration=cal, test=test, base_center=base.center, base_radius_squared=base.radius_squared
    )
    angle = np.linspace(0, 2 * np.pi, 4097)
    direction = np.column_stack([np.cos(angle), np.sin(angle)])
    for name, model in models.items():
        disks = model.disks if hasattr(model, "disks") else [model]
        for index, disk in enumerate(disks):
            boundary = disk.center + np.sqrt(disk.radius_squared) * direction
            active = np.ones(len(boundary), dtype=bool)
            for other in disks:
                tolerance = 1e-10 * max(1.0, other.radius_squared)
                active &= other.scores(boundary) <= tolerance
            boundary[~active] = np.nan
            ax.plot(
                boundary[:, 0],
                boundary[:, 1],
                color=colors[name],
                linewidth=2,
                label=method_label(name) if index == 0 else None,
            )
        indices = np.asarray(removed.get(name, []), dtype=int)
        arrays[name + "_centers"] = np.array([d.center for d in disks])
        arrays[name + "_radii_squared"] = np.array([d.radius_squared for d in disks])
        arrays[name + "_removed"] = indices
    ax.set(aspect="equal", xlabel=r"Residual $z_1$", ylabel=r"Residual $z_2$")
    # Square viewing windows let equal-scale geometry fill each panel rather
    # than squeezing narrow residual distributions into thin axes.
    x_limits, y_limits = ax.get_xlim(), ax.get_ylim()
    span = max(x_limits[1] - x_limits[0], y_limits[1] - y_limits[0])
    x_center, y_center = np.mean(x_limits), np.mean(y_limits)
    ax.set_xlim(x_center - span / 2, x_center + span / 2)
    ax.set_ylim(y_center - span / 2, y_center + span / 2)
    ax.set_title(rf"$T={horizon}$")
    if standalone:
        shared_legend(fig, [ax])
        fig.savefig(directory / "sets.png", dpi=180)
        plt.close(fig)
    np.savez_compressed(directory / f"first_run_T{horizon}.npz", **arrays)


def plot_grouped_sets(directory, first_runs):
    """Overlay methods in one panel per horizon, with independent spatial scales."""
    if not first_runs:
        return
    fig, axes = plt.subplots(
        1, len(first_runs), figsize=(7 * len(first_runs), 8), squeeze=False, constrained_layout=True
    )
    for ax, (horizon, cal, test, models, removed) in zip(axes[0], first_runs):
        plot_sets(directory, horizon, cal, test, models, removed, ax=ax)
    fig.suptitle("First fitted prediction sets")
    shared_legend(fig, axes[0])
    fig.savefig(directory / "sets.png", dpi=180)
    plt.close(fig)


def report_tables(summary, paired):
    lines = [
        "",
        "## Results",
        "",
        "| Horizon | Method | Coverage | Area | Fit ms | Query ms | Failures |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary.itertuples():
        lines.append(
            f"| {row.horizon} | {row.method} | {row.coverage:.5f} | {row.area:.6g} "
            f"| {1000 * row.fit_seconds:.3f} | {1000 * row.query_seconds:.3f} | {row.failures} |"
        )
    lines += [
        "",
        "Mean intervals use independent runs. Horizons share trajectories; "
        "do not pool them as independent runs. Coverage quantiles include test-sample noise.",
        "",
        "## Paired differences",
        "",
        "Differences are method minus reference, using runs where both succeeded.",
        "",
        "| Horizon | Method | Reference | Metric | Difference | 95% interval | Pairs |",
        "|---|---|---|---|---:|---|---:|",
    ]
    for row in paired.itertuples():
        lines.append(
            f"| {row.horizon} | {row.method} | {row.reference} | {row.metric} "
            f"| {row.mean:.6g} | [{row.ci_low:.6g}, {row.ci_high:.6g}] | {row.pairs} |"
        )
    return "\n".join(lines)
