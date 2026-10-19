"""Vehicle benchmark dynamics and predictor for the Section 4 experiments."""

import numpy as np

# Vehicle benchmark constants adapted from Tumu et al. (2024).
DT, L = 0.1, 1.0
PREFIX_SIZE, SUFFIX_SIZE = 5, 50
INITIAL_CONDITIONS = np.array([0.0, 0.0, 0.0, 0.5])
STD_DEVS = [0.3, 0.15 / 2.5]
ACCELERATE = np.array([0.0, 0.15])
COAST = np.array([0.0, 0.0])
LEFT = np.array([0.3, 0.0])
RIGHT = np.array([-0.3, 0.0])
STRAIGHT_ACCEL = [ACCELERATE] * PREFIX_SIZE + [COAST] * SUFFIX_SIZE
STRAIGHT_ACCEL_LEFT = [ACCELERATE] * PREFIX_SIZE + [1.2 * ACCELERATE + LEFT] * SUFFIX_SIZE
STRAIGHT_ACCEL_RIGHT = [ACCELERATE] * PREFIX_SIZE + [1.2 * ACCELERATE + RIGHT] * SUFFIX_SIZE


def predict_trajectories(prefix, steps=50):
    """Same constant speed / average heading increment predictor as the notebook."""
    state = np.asarray(prefix[:, -1], dtype=float).copy()
    turn = (prefix[:, -1, 2] - prefix[:, 0, 2]) / (prefix.shape[1] - 1)
    prediction = np.empty((len(prefix), steps, 4))
    for t in range(steps):
        state[:, 0] += DT * state[:, 3] * np.cos(state[:, 2])
        state[:, 1] += DT * state[:, 3] * np.sin(state[:, 2])
        state[:, 2] += turn
        prediction[:, t] = state
    return prediction


def generate_car_trajectories(n, rng):
    """IID equal-probability left/right/straight policies with independent noise.

    Unlike the reference generator, never reset a global seed between policies
    or datasets: doing so duplicates noise across rollouts used for evaluation.
    """
    if n < 1:
        raise ValueError("n must be positive")
    policies = np.asarray([STRAIGHT_ACCEL_LEFT, STRAIGHT_ACCEL_RIGHT, STRAIGHT_ACCEL])
    controls = policies[rng.integers(3, size=n)].copy()
    controls += rng.normal(size=controls.shape) * np.asarray(STD_DEVS)
    state = np.tile(INITIAL_CONDITIONS, (n, 1))
    trajectory = np.empty((n, controls.shape[1], 4))
    for t in range(controls.shape[1]):
        x, y, heading, speed = state.T.copy()
        state[:, 0] = x + DT * speed * np.cos(heading)
        state[:, 1] = y + DT * speed * np.sin(heading)
        state[:, 2] = heading + DT * speed * np.tan(controls[:, t, 0]) / L
        state[:, 3] = speed + DT * controls[:, t, 1]
        trajectory[:, t] = state
    return trajectory[:, :PREFIX_SIZE], trajectory[:, PREFIX_SIZE:]


def generate_car_residuals(n, rng):
    """Return observed minus predicted XY positions, shape (n, 50, 2)."""
    prefix, future = generate_car_trajectories(n, rng)
    return future[:, :, :2] - predict_trajectories(prefix, future.shape[1])[:, :, :2]
