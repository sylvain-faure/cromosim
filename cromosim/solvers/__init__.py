"""Compiled numerical solvers shipped with :mod:`cromosim`."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .._optimization import Minimizer

__all__ = ["Minimizer", "minimize"]


def minimize(
    xyrv: ArrayLike,
    contacts: ArrayLike,
    desired_velocity: ArrayLike,
    mass: ArrayLike | None = None,
    *,
    max_iterations: int = 100_000,
    rho: float = 0.1,
    tolerance: float = 0.01,
    time_step: float,
    minimum_people_distance: float = 0.0,
    minimum_wall_distance: float = 0.0,
    nesterov: bool = False,
    thread_count: int = 0,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], int]:
    """Solve one velocity projection and return NumPy arrays.

    ``xyrv``, ``contacts`` and ``desired_velocity`` respectively have shapes
    ``(N, 4)``, ``(M, 5)`` and ``(N, 2)``. ``mass`` has shape ``(N,)`` and
    defaults to one for every person. A ``thread_count`` of zero lets the
    solver select the number of worker threads.
    """
    people = _matrix(xyrv, "xyrv", 4)
    contact_data = _matrix(contacts, "contacts", 5)
    velocity = _matrix(desired_velocity, "desired_velocity", 2)
    if velocity.shape[0] != people.shape[0]:
        raise ValueError("desired_velocity must have one row per person")

    if mass is None:
        masses = np.ones(people.shape[0], dtype=np.float64)
    else:
        masses = np.asarray(mass, dtype=np.float64)
        if masses.ndim != 1 or masses.shape[0] != people.shape[0]:
            raise ValueError("mass must have shape (N,)")

    minimizer = Minimizer(
        max_iterations,
        rho,
        tolerance,
        time_step,
        minimum_people_distance,
        minimum_wall_distance,
        nesterov,
    )
    minimizer.set_thread_count(thread_count)
    result = minimizer.run(
        people.ravel().tolist(),
        contact_data.ravel().tolist(),
        velocity.ravel().tolist(),
        masses.tolist(),
    )
    projected, pressure, multipliers, info = result
    return (
        np.asarray(projected, dtype=np.float64).reshape((-1, 2)),
        np.asarray(pressure, dtype=np.float64),
        np.asarray(multipliers, dtype=np.float64),
        info,
    )


def _matrix(values: ArrayLike | Sequence[float], name: str, columns: int) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 1 and array.size == 0:
        return array.reshape((0, columns))
    if array.ndim != 2 or array.shape[1] != columns:
        raise ValueError(f"{name} must have shape (N, {columns})")
    return np.ascontiguousarray(array)
