"""Prediction sets for the three constructions in papers/main.tex.

All methods solve min t subject to ||z_i - c||^2 <= t in two dimensions (d=3):
* fit_cascade_predictor: hat-Gamma, support removal with fixed-label padding;
* fit_threshold_predictor: check-Gamma, one fit followed by residual thresholding;
* fit_intersection_predictor: bar-Gamma, ranked removal and intersection.

Geometry and numerical support checks are shared by all three constructions.
Observation indices in audit records always refer to the original input sample.
"""

import math
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from scipy.integrate import quad
from scipy.spatial import ConvexHull, QhullError

SUPPORT_BOUND = 3  # Two center coordinates and one squared radius.
FEASIBILITY_TOL = 2e-12  # Absolute tolerance after normalizing the input points.
SUPPORT_TOL = 1e-9  # Relative tolerance for support and reconstruction checks.


@dataclass
class Disk:
    """A closed disk, storing its squared radius to match the scenario objective."""

    center: np.ndarray
    radius_squared: float

    def scores(self, z) -> np.ndarray:
        """Constraint residuals g(q; z); nonpositive values lie in the disk."""
        return np.sum((_validate_points(z) - self.center) ** 2, axis=1) - self.radius_squared

    def contains(self, z):
        return self.scores(z) <= 0

    def area(self):
        return float(np.pi * self.radius_squared)


@dataclass
class DiskIntersection:
    """Closed prediction set formed by intersecting the component disks."""

    disks: list[Disk]

    def contains(self, z):
        return np.logical_and.reduce([disk.contains(z) for disk in self.disks])

    def area(self):
        """Integrate vertical slices of the intersection (absolute area units)."""
        centers = np.array([disk.center for disk in self.disks])
        radii = np.sqrt([disk.radius_squared for disk in self.disks])
        left = float(np.max(centers[:, 0] - radii))
        right = float(np.min(centers[:, 0] + radii))
        if left >= right:
            return 0.0

        def height(x):
            half_heights = np.sqrt(np.maximum(0, radii**2 - (x - centers[:, 0]) ** 2))
            return max(
                0.0,
                float(np.min(centers[:, 1] + half_heights) - np.max(centers[:, 1] - half_heights)),
            )

        # Circle crossing abscissae partition the envelope into smooth arcs.
        crossing_x_coordinates = []
        for i, j in combinations(range(len(radii)), 2):
            center_delta = centers[j] - centers[i]
            distance = float(np.linalg.norm(center_delta))
            if distance == 0 or not abs(radii[i] - radii[j]) < distance < radii[i] + radii[j]:
                continue
            distance_along_centers = (radii[i] ** 2 - radii[j] ** 2 + distance**2) / (2 * distance)
            crossing_height = math.sqrt(max(0.0, radii[i] ** 2 - distance_along_centers**2))
            chord_center = centers[i] + distance_along_centers * center_delta / distance
            crossing_x_coordinates.extend(
                [
                    chord_center[0] - crossing_height * center_delta[1] / distance,
                    chord_center[0] + crossing_height * center_delta[1] / distance,
                ]
            )
        integration_limits = sorted(
            set([left, right] + [float(x) for x in crossing_x_coordinates if left < x < right])
        )
        return sum(
            quad(height, a, b, epsabs=1e-9, epsrel=1e-9)[0]
            for a, b in zip(integration_limits[:-1], integration_limits[1:])
        )


@dataclass
class CascadeResult:
    """Final disk and per-stage audit of the support-removal cascade."""

    disk: Disk
    stages: list[dict]
    removed: np.ndarray


def minimum_enclosing_disk(z) -> Disk:
    """Global minimum enclosing disk, up to documented floating-point tolerances.

    Any minimal disk is determined by a diameter pair or three points on the convex
    hull. Isotropic normalization is numerical only and preserves the problem.
    """
    z = _validate_points(z)
    if not len(z):
        raise ValueError("Cannot fit an empty sample")
    origin = z.mean(axis=0)
    scale = float(np.max(np.linalg.norm(z - origin, axis=1)))
    if scale == 0:
        return Disk(origin.copy(), 0.0)
    hull_points = np.unique((z - origin) / scale, axis=0)
    if len(hull_points) > 2:
        try:
            hull_points = hull_points[ConvexHull(hull_points).vertices]
        except QhullError:
            # Collinear inputs are covered by the diameter-pair candidates.
            pass
    centers, squared_radii = _boundary_disk_candidates(hull_points)
    # Stop at the smallest feasible candidate, testing all hull points.
    for index in np.argsort(squared_radii, kind="stable"):
        distances = np.sum((hull_points - centers[index]) ** 2, axis=1)
        if distances.max() <= squared_radii[index] + FEASIBILITY_TOL:
            center = origin + scale * centers[index]
            # Include numerical boundary points without modifying query tolerance.
            radius_squared = float(np.max(np.sum((z - center) ** 2, axis=1)))
            return Disk(center, radius_squared)
    raise RuntimeError("No feasible enclosing disk found")


def find_support_indices(z, disk: Disk) -> np.ndarray:
    """Return local support indices, checking deletion and reconstruction.

    Degenerate samples (e.g. four cocircular points) are rejected rather than
    assigning the nondegenerate-cascade certificate to an arbitrary basis.
    """
    z = _validate_points(z)
    if len(z) == 1:
        return np.array([0], dtype=int)
    radius_scale = max(disk.radius_squared, np.finfo(float).tiny)
    active_indices = np.flatnonzero(np.abs(disk.scores(z)) <= SUPPORT_TOL * radius_scale)
    support = []
    for index in active_indices:
        disk_without_point = minimum_enclosing_disk(np.delete(z, index, axis=0))
        radius_change = disk.radius_squared - disk_without_point.radius_squared
        center_change = np.linalg.norm(disk.center - disk_without_point.center)
        # Support means a change of optimizer, including the center. Near a
        # right triangle, a real center shift can have a tiny objective change.
        if radius_change > SUPPORT_TOL * radius_scale or center_change > SUPPORT_TOL * math.sqrt(
            radius_scale
        ):
            support.append(int(index))
    if not 1 <= len(support) <= SUPPORT_BOUND:
        raise RuntimeError("Degenerate or numerically ambiguous support set")
    reconstructed = minimum_enclosing_disk(z[support])
    radius_error = abs(reconstructed.radius_squared - disk.radius_squared)
    center_error = np.linalg.norm(reconstructed.center - disk.center)
    if radius_error > SUPPORT_TOL * radius_scale or center_error > SUPPORT_TOL * math.sqrt(
        radius_scale
    ):
        raise RuntimeError("Support reconstruction failed at the specified tolerance")
    return np.asarray(support, dtype=int)


def fit_cascade_predictor(
    z, rounds: int, fitted: Disk | None = None, initial_support=None
) -> CascadeResult:
    """Remove all support points plus smallest original row labels to total d.

    Labels stay attached to observations throughout the cascade. Padding never
    uses residual ranks. No deletion occurs after the final solve.
    There are rounds + 1 fits (ell in the paper). Optional fitted and
    initial_support arguments reuse the first fit and its local support indices.
    """
    z = _validate_points(z)
    if (
        not isinstance(rounds, (int, np.integer))
        or rounds < 0
        or len(z) <= SUPPORT_BOUND * (rounds + 1)
    ):
        raise ValueError("Require nonnegative integer rounds and n > d*(rounds+1)")
    retained_indices = np.arange(len(z))
    stages, removed = [], []
    disk = fitted
    for stage in range(rounds + 1):
        if stage > 0 or disk is None:
            disk = minimum_enclosing_disk(z[retained_indices])
        local_support = (
            initial_support
            if stage == 0 and initial_support is not None
            else find_support_indices(z[retained_indices], disk)
        )
        support = retained_indices[local_support]
        padding = np.array([], dtype=int)
        removal_batch = np.array([], dtype=int)
        if stage < rounds:
            padding = np.setdiff1d(retained_indices, support)[: SUPPORT_BOUND - len(support)]
            removal_batch = np.concatenate([support, padding])
        stages.append(
            dict(
                stage=stage,
                retained=int(len(retained_indices)),
                support=support.tolist(),
                padding=padding.tolist(),
                removed=removal_batch.tolist(),
                center=disk.center.tolist(),
                radius_squared=disk.radius_squared,
            )
        )
        removed.extend(removal_batch.tolist())
        retained_indices = np.setdiff1d(retained_indices, removal_batch)
    return CascadeResult(disk, stages, np.asarray(removed, dtype=int))


def fit_threshold_predictor(z, rank: int, fitted: Disk | None = None) -> tuple[Disk, float]:
    """Return the thresholded disk and the rank-th largest fitted residual.

    The fitted center stays fixed. An optional precomputed disk avoids a solve.
    The threshold is a residual value, not the new squared radius.
    """
    z = _validate_points(z)
    if not isinstance(rank, (int, np.integer)) or not SUPPORT_BOUND <= rank <= len(z):
        raise ValueError("SBCP requires integer d <= rank <= n")
    fitted = minimum_enclosing_disk(z) if fitted is None else fitted
    # NumPy partitions in ascending order; rank counts from the largest score.
    order_index = len(z) - rank
    scores = fitted.scores(z)
    threshold = float(np.partition(scores, order_index)[order_index])
    # Compute the equivalent radius directly to avoid cancellation of t + threshold.
    distances = np.sum((z - fitted.center) ** 2, axis=1)
    radius_squared = float(np.partition(distances, order_index)[order_index])
    return Disk(fitted.center.copy(), radius_squared), threshold


def fit_intersection_predictor(z, batches) -> tuple[DiskIntersection, list[dict]]:
    """Paper's bar-Gamma: residual-ranked removal and intersection of all sets.

    Unlike Romao padding, each batch uses residual ranks, not original labels.
    Reject boundary ties that enlarge a batch rather than silently breaking them.
    Each entry of batches is a prescribed rank s_k >= d. Their sum is the
    target rank r. Return the intersection and an audit record for every fit.
    The final batch is recorded but requires no subsequent fit.
    """
    z = _validate_points(z)
    batches = tuple(batches)
    if (
        not batches
        or any(
            not isinstance(batch_size, (int, np.integer)) or batch_size < SUPPORT_BOUND
            for batch_size in batches
        )
        or sum(batches) > len(z)
    ):
        raise ValueError("Require batches >= d with total <= n")
    retained_indices = np.arange(len(z))
    disks, stages = [], []
    for batch_size in batches:
        fitted = minimum_enclosing_disk(z[retained_indices])
        find_support_indices(z[retained_indices], fitted)
        disk, threshold = fit_threshold_predictor(z[retained_indices], batch_size, fitted)
        distances = np.sum((z[retained_indices] - fitted.center) ** 2, axis=1)
        removed = retained_indices[distances >= disk.radius_squared]
        if len(removed) != batch_size:
            raise RuntimeError("Residual ties enlarge the prescribed batch")
        disks.append(disk)
        stages.append(
            dict(
                retained=len(retained_indices),
                batch_size=int(batch_size),
                removed=removed.tolist(),
                threshold=threshold,
                center=fitted.center.tolist(),
                radius_squared=disk.radius_squared,
            )
        )
        retained_indices = np.setdiff1d(retained_indices, removed)
    return DiskIntersection(disks), stages


def _boundary_disk_candidates(hull_points):
    """Enumerate diameter disks and circumdisks of noncollinear hull triples."""
    pairs = np.asarray(list(combinations(range(len(hull_points)), 2)), dtype=int)
    centers = (hull_points[pairs[:, 0]] + hull_points[pairs[:, 1]]) / 2
    squared_radii = np.sum((hull_points[pairs[:, 0]] - centers) ** 2, axis=1)
    if len(hull_points) >= 3:
        triples = np.asarray(list(combinations(range(len(hull_points)), 3)), dtype=int)
        anchors = hull_points[triples[:, 0]]
        first_edges = hull_points[triples[:, 1]] - anchors
        second_edges = hull_points[triples[:, 2]] - anchors
        determinant = (
            first_edges[:, 0] * second_edges[:, 1] - first_edges[:, 1] * second_edges[:, 0]
        )
        noncollinear = np.abs(determinant) > 1e-14
        anchors, first_edges, second_edges, determinant = (
            anchors[noncollinear],
            first_edges[noncollinear],
            second_edges[noncollinear],
            determinant[noncollinear],
        )
        first_half_squared, second_half_squared = (
            np.sum(first_edges * first_edges, axis=1) / 2,
            np.sum(second_edges * second_edges, axis=1) / 2,
        )
        offset = np.column_stack(
            [
                (first_half_squared * second_edges[:, 1] - second_half_squared * first_edges[:, 1])
                / determinant,
                (first_edges[:, 0] * second_half_squared - second_edges[:, 0] * first_half_squared)
                / determinant,
            ]
        )
        centers = np.vstack([centers, anchors + offset])
        squared_radii = np.concatenate([squared_radii, np.sum(offset * offset, axis=1)])
    return centers, squared_radii


def _validate_points(value) -> np.ndarray:
    """Validate and convert observations to a finite array with shape (n, 2)."""
    value = np.asarray(value, dtype=float)
    if value.ndim != 2 or value.shape[1] != 2 or not np.isfinite(value).all():
        raise ValueError("Expected a finite (n, 2) array")
    return value
