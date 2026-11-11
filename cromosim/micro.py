# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
# License: GPL

"""Microscopic crowd-simulation utilities used by Cromosim.

The module contains contact detection, force evaluation, velocity projection,
pedestrian initialization and transfer helpers, sensor processing, and plotting
utilities for the microscopic crowd model.

"""

import logging
import sys

import numpy as np
import matplotlib.pyplot as plt

# from scipy.spatial import KDTree
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse import spdiags
from matplotlib.patches import Arrow
from matplotlib.lines import Line2D
from matplotlib.collections import EllipseCollection, LineCollection


_LOGGER = logging.getLogger(__name__)

PEOPLE_ARRAY_KEYS = (
    "xyrv",
    "Vd",
    "U",
    "Uold",
    "destinations",
    "frustration",
    "with_frustration",
    "beta",
    "viewing_angle",
    "mass",
    "id",
    "frustration_last_pos",
    "frustration_memory",
    "frustration_memory_times",
)
PEOPLE_LIST_KEYS = (
)
PEOPLE_DICT_KEYS = (
    "paths",
)


def compute_contacts(dom, xyrv, dmax):
    """Compute active pedestrian-pedestrian and pedestrian-wall contacts.

    Pedestrian pairs are identified with a :class:`scipy.spatial.cKDTree`.
    Wall contacts are obtained from the domain wall-distance field.

    Parameters
    ----------
    dom : Domain
        Domain object providing wall distances and wall-distance gradients.
    xyrv : numpy.ndarray, shape (n_people, 4)
        Pedestrian state array with columns ``x``, ``y``, radius ``r`` and
        velocity coefficient ``v``.
    dmax : float
        Maximum gap distance for an active contact. A contact is retained when
        ``dij < dmax``.

    Returns
    -------
    numpy.ndarray, shape (n_contacts, 5)
        Contact rows ``(i, j, dij, eij_x, eij_y)``. Pedestrian-pedestrian
        contacts satisfy ``i < j``. Wall contacts use ``j = -1``.
    """
    if xyrv.size>0:
        # lf: the number of points at which the algorithm
        # switches over to brute-force. Has to be positive.
        lf = 100
        if (lf > sys.getrecursionlimit()):
            sys.setrecursionlimit(lf)
        kd = cKDTree(xyrv[:, :2], leafsize=lf)
        # Find all pairs of points whose distance is at most dmax+2*rmax
        rmax = xyrv[:, 2].max()
        neighbors = kd.query_ball_tree(kd, dmax+2*rmax)
        # Create the contact array: i,j,dij,eij_x,eij_y
        first_elements = np.arange(xyrv.shape[0])  # i
        other_elements = list(map(lambda x: x[1:], neighbors))  # all j for each i
        lengths = list(map(len, other_elements))
        tt = np.stack([first_elements, lengths], axis=1)
        II = np.concatenate(list(map(lambda x: np.full((x[1],), x[0]), tt))).astype(int)
        JJ = np.concatenate(other_elements).astype(int)
        ind = np.where(II < JJ)[0]
        II = II[ind]
        JJ = JJ[ind]
        DP = xyrv[JJ, :2]-xyrv[II, :2]
        Norm = np.linalg.norm(DP, axis=1, ord=2)
        Dij = Norm - xyrv[II, 2]-xyrv[JJ, 2]
        ind = np.where(Dij < dmax)[0]
        Dij = Dij[ind]
        II = II[ind]
        JJ = JJ[ind]
        Norm = Norm[ind]
        DP = DP[ind]
        contacts = np.stack([II, JJ, Dij, DP[:, 0]/Norm, DP[:, 1]/Norm], axis=1)
        # Add contacts with the walls
        II, JJ, DD = dom.people_wall_distance(xyrv)
        ind = np.where(DD < dmax)[0]
        wall_contacts = np.stack([ind, -1*np.ones(ind.shape), DD[ind],
                                dom.wall_grad_X[II[ind], JJ[ind]],
                                dom.wall_grad_Y[II[ind], JJ[ind]]], axis=1)
        contacts = np.vstack([contacts, wall_contacts])
        return np.array(contacts)
    else:
        return []


def compute_forces(F, Fwall, xyrv, contacts, U, Vd, lambda_, delta, k, eta):
    """Compute social, overlap and friction forces for all pedestrians.

    The social-force contribution includes the directional weighting associated
    with the pedestrian viewing direction. Overlap and tangential friction
    forces are added for active contacts.

    Parameters
    ----------
    F : float
        Pedestrian-pedestrian social-force amplitude.
    Fwall : float
        Pedestrian-wall social-force amplitude.
    xyrv : numpy.ndarray, shape (n_people, 4)
        Pedestrian state array with columns ``x``, ``y``, radius ``r`` and
        velocity coefficient ``v``.
    contacts : numpy.ndarray, shape (n_contacts, 5)
        Contact rows ``(i, j, dij, eij_x, eij_y)``.
    U : numpy.ndarray, shape (n_people, 2)
        Current pedestrian velocities.
    Vd : numpy.ndarray, shape (n_people, 2)
        Desired pedestrian velocities.
    lambda_ : float
        Directional-dependence coefficient in ``[0, 1]``. The value ``1``
        gives an isotropic interaction.
    delta : float
        Characteristic interaction distance in the exponential social force.
    k : float
        Contact stiffness used when pedestrians overlap.
    eta : float
        Tangential friction coefficient.

    Returns
    -------
    numpy.ndarray, shape (n_people, 2)
        Total force applied to each pedestrian.
    """
    Np = xyrv.shape[0]
    Nc = contacts.shape[0]
    Forces = np.zeros((Np, 2))
    # Social forces, friction,...
    for ic in np.arange(Nc):
        i = int(contacts[ic, 0])
        j = int(contacts[ic, 1])
        dij = contacts[ic, 2]
        dij_moins = min(dij, 0.0)
        eij_x = contacts[ic, 3]
        eij_y = contacts[ic, 4]
        if (j > -1):  # contact person/person
            # Angular dependence
            norm_Vdi = np.sqrt(Vd[i, 0]**2+Vd[i, 1]**2)
            if (norm_Vdi > 0):
                theta_ij = np.arccos( max(-1.0, min(1.0, ( (Vd[i,0]*eij_x+Vd[i,1]*eij_y)/norm_Vdi ))))
                omega_ij = lambda_+(1-lambda_)*(1+np.cos(theta_ij))/2
            else:
                omega_ij = 1
            norm_Vdj = np.sqrt(Vd[j, 0]**2+Vd[j, 1]**2)
            if (norm_Vdj > 0):
                theta_ji = np.arccos( max(-1.0, min(1.0, ( -(Vd[j,0]*eij_x+Vd[j,1]*eij_y)/norm_Vdj ))))
                omega_ji = lambda_+(1-lambda_)*(1+np.cos(theta_ji))/2
            else:
                omega_ji = 1
            # Social force + force to handle overlapping
            fij = -omega_ij*F*np.exp(-dij/delta) + k*dij_moins
            fji = -omega_ji*F*np.exp(-dij/delta) + k*dij_moins
            Forces[i, 0] += fij*eij_x
            Forces[i, 1] += fij*eij_y
            Forces[j, 0] -= fji*eij_x
            Forces[j, 1] -= fji*eij_y
            # Friction
            fij_friction = eta*dij_moins*(-(U[i, 0]-U[j, 0])*eij_y +
                                           (U[i, 1]-U[j, 1])*eij_x)
            Forces[i, 0] -= fij_friction*eij_y
            Forces[i, 1] += fij_friction*eij_x
            Forces[j, 0] += fij_friction*eij_y
            Forces[j, 1] -= fij_friction*eij_x
        else:  # contact person/walls
            fij = -Fwall*np.exp(-dij/delta) + k*dij_moins
            Forces[i, 0] -= fij*eij_x
            Forces[i, 1] -= fij*eij_y
    return Forces


def _build_contact_matrix(contacts, n_people):
    """Build the sparse linearized non-overlap constraint matrix."""
    ii = contacts[:, 0].astype(int)
    jj = contacts[:, 1].astype(int)
    person_contacts = np.where(jj >= 0)[0]
    wall_contacts = np.where(jj < 0)[0]

    row = np.concatenate([
        person_contacts,
        person_contacts,
        person_contacts,
        person_contacts,
        wall_contacts,
        wall_contacts,
    ])
    col = np.concatenate([
        2 * ii[person_contacts],
        2 * ii[person_contacts] + 1,
        2 * jj[person_contacts],
        2 * jj[person_contacts] + 1,
        2 * ii[wall_contacts],
        2 * ii[wall_contacts] + 1,
    ])
    data = np.concatenate([
        contacts[person_contacts, 3],
        contacts[person_contacts, 4],
        -contacts[person_contacts, 3],
        -contacts[person_contacts, 4],
        -contacts[wall_contacts, 3],
        -contacts[wall_contacts, 4],
    ])
    matrix = csr_matrix((data, (row, col)), shape=(contacts.shape[0], 2 * n_people))
    return matrix, ii, jj, person_contacts, wall_contacts, row, col, data


def _run_uzawa_projection(
    dt,
    xyrv,
    contacts,
    desired_velocity,
    mass,
    constraint_matrix,
    dmin_people,
    dmin_walls,
    nb_iter_max,
    rho,
    tol,
):
    """Run the original Uzawa iteration used by :func:`projection`."""
    n_people = xyrv.shape[0]
    n_contacts = contacts.shape[0]

    multipliers = np.zeros((n_contacts,))
    residual = 99 * np.ones((n_contacts,))
    velocity = np.zeros((2 * n_people,))
    target_velocity = np.zeros((2 * n_people,))
    distances = contacts[:, 2]
    target_velocity[::2] = desired_velocity[:, 0]
    target_velocity[1::2] = desired_velocity[:, 1]

    iteration = 0
    wall_contacts = np.where(contacts[:, 1] == -1)[0]
    minimum_distance = np.ones(contacts.shape[0]) * dmin_people
    minimum_distance[wall_contacts] = dmin_walls

    mass_vector = np.zeros((2 * n_people,))
    mass_vector[::2] = mass
    mass_vector[1::2] = mass

    while (
        dt * residual.max() > tol * 2 * xyrv[:, 2].min()
        and iteration < nb_iter_max
    ):
        velocity[:] = (
            target_velocity[:]
            - constraint_matrix.transpose() @ multipliers[:] / mass_vector
        )
        residual[:] = (
            constraint_matrix @ velocity[:]
            - (distances[:] - minimum_distance) / dt
        )
        multipliers[:] = np.maximum(multipliers[:] + rho * residual[:], 0)
        iteration += 1

    return velocity, multipliers, residual, iteration


def projection(dt, xyrv, contacts, Vd, Mass, dmin_people=0, dmin_walls=0,
               nb_iter_max=100000, rho=0.1, tol=0.01, log=False, method="cvxopt"):
    """Project desired velocities onto the non-overlapping feasible set.

    The projection minimizes the mass-weighted distance to the desired velocity
    under linearized non-overlap constraints. Several numerical solvers are
    supported, including CVXOPT, MOSEK, Uzawa and the C++ Cromosim solvers.

    Parameters
    ----------
    dt : float
        Time step.
    xyrv : numpy.ndarray, shape (n_people, 4)
        Pedestrian state array with columns ``x``, ``y``, radius ``r`` and
        velocity coefficient ``v``.
    contacts : numpy.ndarray, shape (n_contacts, 5)
        Contact rows ``(i, j, dij, eij_x, eij_y)``.
    Vd : numpy.ndarray, shape (n_people, 2)
        Desired pedestrian velocities.
    Mass : numpy.ndarray, shape (n_people,)
        Pedestrian masses.
    dmin_people : float, optional
        Minimum pedestrian-pedestrian distance.
    dmin_walls : float, optional
        Minimum pedestrian-wall distance.
    nb_iter_max : int, optional
        Maximum number of iterations for iterative solvers.
    rho : float, optional
        Step parameter used by the Uzawa method.
    tol : float, optional
        Solver tolerance.
    log : bool, optional
        If ``True``, emit solver diagnostics.
    method : str, optional
        Solver name. Supported values are ``"cvxopt"``, ``"mosek"``,
        ``"uzawa"``, ``"cppuzawa"``, ``"nesterov"`` and
        ``"cppnesterov"``.

    Returns
    -------
    info : int
        Solver information, typically the number of iterations. ``-1`` denotes
        an iterative solver that reached its iteration limit.
    B : scipy.sparse.csr_matrix or None
        Contact-constraint matrix.
    U : numpy.ndarray, shape (n_people, 2)
        Projected pedestrian velocities.
    L : numpy.ndarray or None
        Contact Lagrange multipliers when available.
    P : numpy.ndarray or None
        Per-pedestrian pressure estimate when multipliers are available.
    frustrations : numpy.ndarray, shape (n_people,)
        Frustration measure in ``[0, 1]`` based on the loss of velocity in the
        desired direction.
    """
    Np = xyrv.shape[0]
    Nc = contacts.shape[0]
    info = 0
    L = None
    P = None
    B = None

    if Nc == 0:
        info = 1
        frustrations = np.zeros((Np,))
        return info, None, Vd, L, P, frustrations

    B, II, JJ, Jpos, Jneg, row, col, data = _build_contact_matrix(contacts, Np)

    if method in ("cvxopt", "mosek"):
        import cvxopt

        cvxopt.solvers.options["show_progress"] = False
        cvxopt.solvers.options["maxiters"] = 1000
        cvxopt.solvers.options["abstol"] = 1e-8
        cvxopt.solvers.options["reltol"] = 1e-7

        V = np.zeros((2 * Np,))
        wall_contacts = np.where(contacts[:, 1] == -1)[0]
        DMIN = np.ones(contacts.shape[0]) * dmin_people
        DMIN[wall_contacts] = dmin_walls
        Z = (contacts[:, 2] - DMIN) / dt
        V[::2] = Vd[:, 0]
        V[1::2] = Vd[:, 1]

        V_cvx = cvxopt.matrix(V)
        Z_cvx = cvxopt.matrix(Z, (Nc, 1))
        M_cvx = cvxopt.spdiag(cvxopt.matrix(np.repeat(Mass, 2)))
        cvxoptB = cvxopt.spmatrix(
            np.array(data), np.array(row), np.array(col), size=(Nc, 2 * Np)
        )

        if method == "mosek":
            from mosek import iparam

            cvxopt.solvers.options["mosek"] = {iparam.log: 0}
            solution = cvxopt.solvers.qp(
                M_cvx, -M_cvx * V_cvx, cvxoptB, Z_cvx, solver="mosek"
            )
        else:
            solution = cvxopt.solvers.qp(M_cvx, -M_cvx * V_cvx, cvxoptB, Z_cvx)
            info = solution["iterations"]

        if solution["status"] == "optimal":
            U = np.asarray(solution["x"], dtype=np.float64)
            if log:
                constraint_margin = np.asarray(Z_cvx - B @ U)
                if method == "mosek":
                    _LOGGER.info(
                        "projection[%s]: contacts=%d, min(Z - B@U)=%.6e",
                        method,
                        Nc,
                        float(constraint_margin.min()),
                    )
                else:
                    _LOGGER.info(
                        "projection[%s]: contacts=%d, iterations=%s, status=%s, "
                        "min(Z - B@U)=%.6e",
                        method,
                        Nc,
                        solution["iterations"],
                        solution["status"],
                        float(constraint_margin.min()),
                    )
        else:
            _LOGGER.warning(
                "projection[%s] failed with status=%s; falling back to Uzawa",
                method,
                solution["status"],
            )
            info = 0
            U, L, R, k = _run_uzawa_projection(
                dt,
                xyrv,
                contacts,
                Vd,
                Mass,
                B,
                dmin_people,
                dmin_walls,
                nb_iter_max,
                rho,
                tol,
            )
            if log:
                _LOGGER.info(
                    "projection[uzawa fallback]: contacts=%d, iterations=%d, "
                    "residual_min=%.6e, residual_max=%.6e, tol=%.6e",
                    Nc,
                    k,
                    float(R.min()),
                    float(R.max()),
                    tol,
                )
            if k == nb_iter_max:
                _LOGGER.warning(
                    "projection[uzawa fallback] reached the maximum number of "
                    "iterations; some constraints remain unsatisfied"
                )
                info = -1
            else:
                info = k

    elif method == "uzawa":
        U, L, R, k = _run_uzawa_projection(
            dt,
            xyrv,
            contacts,
            Vd,
            Mass,
            B,
            dmin_people,
            dmin_walls,
            nb_iter_max,
            rho,
            tol,
        )
        if log:
            _LOGGER.info(
                "projection[uzawa]: contacts=%d, iterations=%d, "
                "residual_min=%.6e, residual_max=%.6e, tol=%.6e",
                Nc,
                k,
                float(R.min()),
                float(R.max()),
                tol,
            )
        if k == nb_iter_max:
            _LOGGER.warning(
                "projection[uzawa] reached the maximum number of iterations; "
                "some constraints remain unsatisfied"
            )
            info = -1
        else:
            info = k

    # elif method == "cppuzawa":
    #     from cromosim.solvers import Minimizer

    #     minimizer = Minimizer(
    #         max_iterations=nb_iter_max,
    #         tolerance=tol,
    #         rho=rho,
    #         time_step=dt,
    #         minimum_people_distance=dmin_people,
    #         minimum_wall_distance=dmin_walls,
    #         nesterov=False,
    #     )
    #     threads = 4
    #     minimizer.set_thread_count(threads)
    #     U, P, L, info = minimizer.run(
    #         xyrv.ravel().tolist(),
    #         contacts.ravel().tolist(),
    #         Vd.ravel().tolist(),
    #         Mass.tolist(),
    #     )
    #     U = np.asarray(U, dtype=np.float64)
    #     L = np.asarray(L, dtype=np.float64)
    #     P = np.asarray(P, dtype=np.float64)

    elif method == "nesterov":
        B = -dt * B
        L = np.zeros((Nc,))
        R = 99 * np.ones((Nc,))
        U = np.zeros((2 * Np,))
        V = np.zeros((2 * Np,))
        V[::2] = Vd[:, 0]
        V[1::2] = Vd[:, 1]
        D = contacts[:, 2]
        k = 0
        wall_contacts = np.where(contacts[:, 1] == -1)[0]
        DMIN = np.ones((Nc,)) * dmin_people
        DMIN[wall_contacts] = dmin_walls

        mass_vector = np.zeros((2 * Np,))
        mass_vector[::2] = Mass
        mass_vector[1::2] = Mass
        M = spdiags(
            np.array([mass_vector]), np.array([0]), mass_vector.shape[0], mass_vector.shape[0]
        ).tocsc()
        Minv = spdiags(
            np.array([1 / mass_vector]),
            np.array([0]),
            mass_vector.shape[0],
            mass_vector.shape[0],
        ).tocsc()

        # g(L) = -0.5 * L.T @ A @ L - L.T @ C,
        # where A = B @ Minv @ B.T and C = B @ V + (D - DMIN).
        # The code below minimizes the opposite convex quadratic over L >= 0.
        A = B @ Minv @ B.T
        C = B @ V + (D - DMIN)

        theta_old = 1
        rho = 0.05
        Lold = L.copy()
        Y = Lold.copy()
        Lchap = np.ones((Lold.shape[0],))
        lip = np.linalg.norm(A @ (Lold - Lchap)) / np.linalg.norm(Lold - Lchap)
        rho = 1 / lip

        while k < nb_iter_max:
            DG = A @ Y + C
            L = np.maximum(Y - rho * DG, 0.0)
            cccc = 0
            lip_itermax = 10
            while (
                cccc < lip_itermax
                and 0.5 * L.T @ A @ L + L.T @ C
                >= 0.5 * Y.T @ A @ Y
                + Y.T @ C
                + DG.T @ (L - Y)
                + 0.5 * lip * np.sum((L - Y) ** 2)
            ):
                lip = 2 * lip
                rho = 1 / lip
                L = np.maximum(Y - rho * DG, 0)
                cccc += 1

            delta = theta_old**2 * (theta_old**2 + 4)
            theta = 0.5 * (-theta_old**2 + np.sqrt(delta))
            beta = theta_old * (1 - theta_old) / (theta_old**2 + theta)
            Y = L + beta * (L - Lold)
            if DG.T @ (L - Lold) > 0:
                Y = L.copy()
                theta = 1.0
            lip = 0.9 * lip
            rho = 1 / lip
            k += 1
            Lold = L.copy()
            theta_old = theta

            R = np.array(-A @ L - C)
            if R.max() <= tol * 2 * xyrv[:, 2].min():
                break

        U = V + Minv @ B.T @ L
        info = k
        if log:
            _LOGGER.info(
                "projection[nesterov]: contacts=%d, iterations=%d, "
                "residual_min=%.6e, residual_max=%.6e, tol=%.6e",
                Nc,
                k,
                float(R.min()),
                float(R.max()),
                tol,
            )
        if k == nb_iter_max:
            _LOGGER.warning(
                "projection[nesterov] reached the maximum number of iterations; "
                "some constraints remain unsatisfied"
            )
            info = -1
        else:
            info = k

    elif method == "cppnesterov":
        from cromosim.solvers import Minimizer

        minimizer = Minimizer(
            max_iterations=nb_iter_max,
            tolerance=tol,
            rho=rho,
            time_step=dt,
            minimum_people_distance=dmin_people,
            minimum_wall_distance=dmin_walls,
            nesterov=True,
        )
        threads = 4
        minimizer.set_thread_count(threads)
        U, P, L, info = minimizer.run(
            xyrv.ravel().tolist(),
            contacts.ravel().tolist(),
            Vd.ravel().tolist(),
            Mass.tolist(),
        )
        U = np.asarray(U, dtype=np.float64)
        L = np.asarray(L, dtype=np.float64)
        P = np.asarray(P, dtype=np.float64)

    UU = U.reshape((Np, 2))

    # Frustration is 1 - <U_i, Vd_i> / <Vd_i, Vd_i>, clipped to [0, 1].
    satisfactions = np.minimum(
        np.maximum(
            (Vd[:, 0] * UU[:, 0] + Vd[:, 1] * UU[:, 1])
            / np.maximum(Vd[:, 0] ** 2 + Vd[:, 1] ** 2, 1.0e-8),
            0.0,
        ),
        1.0,
    )
    frustrations = 1 - satisfactions

    if L is not None:
        # Pressure is estimated from one half of the stress trace in 2D.
        P = np.zeros(Np)
        P[II[Jpos]] += 1 / (2 * np.pi * xyrv[II[Jpos], 2]) * L[Jpos]
        P[JJ[Jpos]] += 1 / (2 * np.pi * xyrv[JJ[Jpos], 2]) * L[Jpos]
        P[II[Jneg]] += 1 / (2 * np.pi * xyrv[II[Jneg], 2]) * L[Jneg]

    return info, B, UU, L, P, frustrations


def move_people(time, dt, people, sensors):
    """Advance pedestrian positions and update crossing sensors.

    Parameters
    ----------
    time : float
        Current simulation time.
    dt : float
        Time step.
    people : dict
        Pedestrian data. The dictionary must contain at least ``"xyrv"``,
        ``"U"`` and ``"id"``.
    sensors : list of dict
        Sensor definitions and accumulated sensor data.

    Returns
    -------
    people : dict
        Updated pedestrian data after advancing positions by ``dt``.
    sensors : list of dict
        Updated sensors. Crossing identifiers, times, intersection points and
        directions are appended to the corresponding sensor records.
    """
    if (len(sensors) > 0):
        xyrv_old = people["xyrv"].copy()
    # Update people coordinates
    people["xyrv"][:, 0] += dt*people["U"][:, 0]
    people["xyrv"][:, 1] += dt*people["U"][:, 1]
    for il, s in enumerate(sensors):
        ll = s["line"]
        # line [x0,y0,x1,y1],
        dist_in = 5.0
        dist_out = 5.0
        id, pt, dir, times, entries, exits, id_people_in, id_people_out = sensor(ll, xyrv_old[:, :2],
                                                    people["xyrv"][:, :2], 
                                                    time, time+dt,
                                                    sensors[il]["dist_in"], sensors[il]["dist_out"])
        sensors[il]["id"]    = sensors[il]["id"] + people["id"][id].tolist()
        sensors[il]["times"] = sensors[il]["times"] + times.tolist()
        sensors[il]["xy"]    = sensors[il]["xy"] + pt.tolist()
        sensors[il]["dir"]   = sensors[il]["dir"] + dir.tolist()
        sensors[il]["id_people_in"]  = id_people_in
        sensors[il]["id_people_out"] = id_people_out
    return people, sensors


def create_people_in_box(nn, box, dest_name, radius_distribution,
                         velocity_distribution, viewing_angle_distribution,
                         mass_distribution, dom, rng, log=True):
    """Create pedestrians at random admissible positions inside a box.

    Radii, velocity coefficients, viewing angles and masses are sampled from
    the requested distributions. Positions are resampled until they belong to a
    region where the selected destination is admissible. Overlaps between
    pedestrians are not removed by this function.

    Parameters
    ----------
    nn : int
        Number of pedestrians to create.
    box : sequence of float
        Bounding box ``[xmin, xmax, ymin, ymax]``.
    dest_name : str
        Initial destination assigned to every created pedestrian.
    radius_distribution : sequence
        Distribution specification ``["uniform", min, max]`` or
        ``["normal", mean, std]`` for pedestrian radii.
    velocity_distribution : sequence
        Distribution specification for velocity coefficients.
    viewing_angle_distribution : sequence
        Distribution specification for half viewing angles, in degrees.
    mass_distribution : sequence
        Distribution specification for pedestrian masses.
    dom : Domain
        Domain object used to test destination admissibility.
    rng : numpy.random.RandomState
        Random-number generator.
    log : bool, optional
        If ``True``, emit initialization diagnostics.

    Returns
    -------
    xyrv : numpy.ndarray, shape (nn, 4)
        New pedestrian states ``(x, y, r, v)``.
    dest_names : list of str
        Destination name for each pedestrian.
    view : numpy.ndarray, shape (nn,)
        Half viewing angle for each pedestrian, clipped to ``[0, 90]`` degrees.
    mass : numpy.ndarray, shape (nn,)
        Pedestrian masses.
    """
    if log:
        _LOGGER.info(
            "create_people_in_box: count=%d, box=%s; initial overlaps may occur",
            nn,
            box,
        )

    xyrv = np.zeros((nn, 4))  # columns: x, y, radius, velocity coefficient
    dest_names = [dest_name] * nn

    if radius_distribution[0] == "uniform":
        rmin = radius_distribution[1]
        rmax = radius_distribution[2]
        xyrv[:, 2] = rng.uniform(rmin, rmax, nn)
    elif radius_distribution[0] == "normal":
        mean = radius_distribution[1]
        sigma = radius_distribution[2]
        xyrv[:, 2] = rng.normal(mean, sigma, nn)
    else:
        _LOGGER.warning(
            "create_people_in_box: unknown radius distribution %r",
            radius_distribution[0],
        )
    p_rmax = xyrv[:, 2].max()

    if velocity_distribution[0] == "uniform":
        vmin = velocity_distribution[1]
        vmax = velocity_distribution[2]
        xyrv[:, 3] = rng.uniform(vmin, vmax, nn)
    elif velocity_distribution[0] == "normal":
        mean = velocity_distribution[1]
        sigma = velocity_distribution[2]
        xyrv[:, 3] = rng.normal(mean, sigma, nn)
    else:
        _LOGGER.warning(
            "create_people_in_box: unknown velocity distribution %r",
            velocity_distribution[0],
        )

    if viewing_angle_distribution[0] == "uniform":
        wmin = viewing_angle_distribution[1]
        wmax = viewing_angle_distribution[2]
        if wmax > 90:
            _LOGGER.warning(
                "create_people_in_box: viewing angle exceeds 90 degrees; "
                "values are clipped to the valid half-angle range [0, 90]"
            )
        view = np.maximum(np.minimum(rng.uniform(wmin, wmax, nn), 90.0), 0.0)
    elif viewing_angle_distribution[0] == "normal":
        mean = viewing_angle_distribution[1]
        sigma = viewing_angle_distribution[2]
        view = rng.normal(mean, sigma, nn)
        if view.max() > 90:
            _LOGGER.warning(
                "create_people_in_box: viewing angle exceeds 90 degrees; "
                "values are clipped to the valid half-angle range [0, 90]"
            )
        view = np.maximum(np.minimum(view, 90.0), 0.0)
    else:
        _LOGGER.warning(
            "create_people_in_box: unknown viewing-angle distribution %r",
            viewing_angle_distribution[0],
        )

    if mass_distribution[0] == "uniform":
        mass_min = mass_distribution[1]
        mass_max = mass_distribution[2]
        mass = rng.uniform(mass_min, mass_max, nn)
    elif mass_distribution[0] == "normal":
        mass_mean = mass_distribution[1]
        mass_sigma = mass_distribution[2]
        mass = rng.normal(mass_mean, mass_sigma, nn)
    else:
        _LOGGER.warning(
            "create_people_in_box: unknown mass distribution %r",
            mass_distribution[0],
        )

    xmin, xmax, ymin, ymax = box
    xyrv[:, 0] = rng.uniform(xmin + p_rmax, xmax - p_rmax, nn)
    xyrv[:, 1] = rng.uniform(ymin + p_rmax, ymax - p_rmax, nn)

    # Resample positions located in destination-specific exclusion zones.
    while True:
        II, JJ, _ = dom.people_desired_velocity(xyrv, dest_names)
        invalid = np.where(dom.destinations[dest_name].excluded_id[II, JJ] > 0)[0]
        if invalid.shape[0] == 0:
            break
        xyrv[invalid, 0] = rng.uniform(
            xmin + p_rmax, xmax - p_rmax, invalid.shape[0]
        )
        xyrv[invalid, 1] = rng.uniform(
            ymin + p_rmax, ymax - p_rmax, invalid.shape[0]
        )

    return xyrv, dest_names, view, mass


def create_people_at_positions(nn, rgb_color, dest_name, radius_distribution,
                               velocity_distribution, viewing_angle_distribution,
                               mass_distribution, dom, rng, log=True):
    """Create pedestrians at positions encoded by an RGB map color.

    Pixels matching ``rgb_color`` define the available positions. If more
    positions than pedestrians are available, a random subset is selected.

    Parameters
    ----------
    nn : int
        Number of pedestrians to create.
    rgb_color : sequence of int
        RGB triplet identifying admissible pixels in the domain image.
    dest_name : str
        Initial destination assigned to every created pedestrian.
    radius_distribution : sequence
        Distribution specification ``["uniform", min, max]`` or
        ``["normal", mean, std]`` for pedestrian radii.
    velocity_distribution : sequence
        Distribution specification for velocity coefficients.
    viewing_angle_distribution : sequence
        Distribution specification for half viewing angles, in degrees.
    mass_distribution : sequence
        Distribution specification for pedestrian masses.
    dom : Domain
        Domain containing the RGB image and geometric metadata.
    rng : numpy.random.RandomState
        Random-number generator.
    log : bool, optional
        If ``True``, emit initialization diagnostics.

    Returns
    -------
    xyrv : numpy.ndarray, shape (nn, 4)
        New pedestrian states ``(x, y, r, v)``.
    dest_names : list of str
        Destination name for each pedestrian.
    view : numpy.ndarray, shape (nn,)
        Half viewing angle for each pedestrian, clipped to ``[0, 90]`` degrees.
    mass : numpy.ndarray, shape (nn,)
        Pedestrian masses.
    """
    if log:
        _LOGGER.info(
            "create_people_at_positions: count=%d, rgb_color=%s",
            nn,
            rgb_color,
        )

    xyrv = np.zeros((nn, 4))  # columns: x, y, radius, velocity coefficient
    dest_names = [dest_name] * nn

    if radius_distribution[0] == "uniform":
        rmin = radius_distribution[1]
        rmax = radius_distribution[2]
        xyrv[:, 2] = rng.uniform(rmin, rmax, nn)
    elif radius_distribution[0] == "normal":
        mean = radius_distribution[1]
        sigma = radius_distribution[2]
        xyrv[:, 2] = rng.normal(mean, sigma, nn)
    else:
        _LOGGER.warning(
            "create_people_at_positions: unknown radius distribution %r",
            radius_distribution[0],
        )

    if velocity_distribution[0] == "uniform":
        vmin = velocity_distribution[1]
        vmax = velocity_distribution[2]
        xyrv[:, 3] = rng.uniform(vmin, vmax, nn)
    elif velocity_distribution[0] == "normal":
        mean = velocity_distribution[1]
        sigma = velocity_distribution[2]
        xyrv[:, 3] = rng.normal(mean, sigma, nn)
    else:
        _LOGGER.warning(
            "create_people_at_positions: unknown velocity distribution %r",
            velocity_distribution[0],
        )

    if viewing_angle_distribution[0] == "uniform":
        wmin = viewing_angle_distribution[1]
        wmax = viewing_angle_distribution[2]
        if wmax > 90:
            _LOGGER.warning(
                "create_people_at_positions: viewing angle exceeds 90 degrees; "
                "values are clipped to the valid half-angle range [0, 90]"
            )
        view = np.maximum(np.minimum(rng.uniform(wmin, wmax, nn), 90.0), 0.0)
    elif viewing_angle_distribution[0] == "normal":
        mean = viewing_angle_distribution[1]
        sigma = viewing_angle_distribution[2]
        view = rng.normal(mean, sigma, nn)
        if view.max() > 90:
            _LOGGER.warning(
                "create_people_at_positions: viewing angle exceeds 90 degrees; "
                "values are clipped to the valid half-angle range [0, 90]"
            )
        view = np.maximum(np.minimum(view, 90.0), 0.0)
    else:
        _LOGGER.warning(
            "create_people_at_positions: unknown viewing-angle distribution %r",
            viewing_angle_distribution[0],
        )

    if mass_distribution[0] == "uniform":
        mass_min = mass_distribution[1]
        mass_max = mass_distribution[2]
        mass = rng.uniform(mass_min, mass_max, nn)
    elif mass_distribution[0] == "normal":
        mass_mean = mass_distribution[1]
        mass_sigma = mass_distribution[2]
        mass = rng.normal(mass_mean, mass_sigma, nn)
    else:
        _LOGGER.warning(
            "create_people_at_positions: unknown mass distribution %r",
            mass_distribution[0],
        )

    positions_mask = (
        (dom.image_red == rgb_color[0])
        * (dom.image_green == rgb_color[1])
        * (dom.image_blue == rgb_color[2])
    )
    mask_id = np.where(positions_mask >= 1)
    if nn > mask_id[0].size:
        raise ValueError(
            "create_people_at_positions: too many pedestrians for the number "
            "of available colored pixels"
        )

    if nn != mask_id[0].size:
        indices = np.arange(mask_id[0].size)
        rng.shuffle(indices)
        mask_id = (mask_id[0][indices[:nn]], mask_id[1][indices[:nn]])

    xyrv[:, 0] = dom.ymin + dom.pixel_size * mask_id[1]
    xyrv[:, 1] = dom.xmin + dom.pixel_size * mask_id[0]

    return xyrv, dest_names, view, mass


def check_people_in_box(dom, box, xyrv, dest_names, rng, log=True):
    """Relocate pedestrians that are outside a box or an admissible region.

    The function repeatedly resamples invalid positions until every pedestrian
    lies inside the computational domain, inside ``box`` and outside destination
    exclusion zones.

    Parameters
    ----------
    dom : Domain
        Domain object used for geometry and destination checks.
    box : sequence of float
        Bounding box ``[xmin, xmax, ymin, ymax]``.
    xyrv : numpy.ndarray, shape (n_people, 4)
        Pedestrian state array.
    dest_names : sequence of str
        Destination name for each pedestrian.
    rng : numpy.random.RandomState
        Random-number generator.
    log : bool, optional
        If ``True``, emit relocation diagnostics.

    Returns
    -------
    info : bool
        ``True`` if at least one pedestrian had to be relocated, otherwise
        ``False``.
    xyrv : numpy.ndarray
        Updated pedestrian states. The input array is modified in place.
    """
    if log:
        _LOGGER.info(
            "check_people_in_box: checking %d pedestrians in box %s",
            xyrv.shape[0],
            box,
        )

    p_rmax = xyrv[:, 2].max()
    xmin, xmax, ymin, ymax = box
    info = False

    while True:
        # Test 1: positions must map to valid domain pixels.
        II = np.floor(
            (xyrv[:, 1] - dom.ymin - 0.5 * dom.pixel_size) / dom.pixel_size
        ).astype(int)
        JJ = np.floor(
            (xyrv[:, 0] - dom.xmin - 0.5 * dom.pixel_size) / dom.pixel_size
        ).astype(int)
        test1 = (II >= 0) * (II < dom.height) * (JJ >= 0) * (JJ < dom.width)
        ind1 = np.where(test1 == 0)[0]
        if ind1.shape[0] > 0:
            if log:
                _LOGGER.info(
                    "check_people_in_box: relocating %d pedestrians outside the domain",
                    ind1.shape[0],
                )
            info = True
            xyrv[ind1, 0] = rng.uniform(
                xmin + p_rmax, xmax - p_rmax, ind1.shape[0]
            )
            xyrv[ind1, 1] = rng.uniform(
                ymin + p_rmax, ymax - p_rmax, ind1.shape[0]
            )
            continue

        # Test 2: positions must remain inside the requested box and have a
        # non-zero desired velocity.
        II, JJ, Vd = dom.people_desired_velocity(xyrv, dest_names)
        normVd = Vd[:, 0] ** 2 + Vd[:, 1] ** 2
        test2 = (
            (xyrv[:, 0] > xmin + p_rmax)
            * (xyrv[:, 0] < xmax - p_rmax)
            * (xyrv[:, 1] > ymin + p_rmax)
            * (xyrv[:, 1] < ymax - p_rmax)
            * (normVd > 0)
        )
        ind2 = np.where(test2 == 0)[0]
        if ind2.shape[0] > 0:
            if log:
                _LOGGER.info(
                    "check_people_in_box: relocating %d pedestrians with an "
                    "undefined desired velocity or outside the box",
                    ind2.shape[0],
                )
            info = True
            xyrv[ind2, 0] = rng.uniform(
                xmin + p_rmax, xmax - p_rmax, ind2.shape[0]
            )
            xyrv[ind2, 1] = rng.uniform(
                ymin + p_rmax, ymax - p_rmax, ind2.shape[0]
            )
            continue

        # Test 3: destination-specific exclusion zones must be avoided.
        II, JJ, _ = dom.people_desired_velocity(xyrv, dest_names)
        ind3 = np.array([], dtype=int)
        destination_array = np.asarray(dest_names)
        for dest in np.unique(dest_names):
            ind_p = np.where(destination_array == dest)[0]
            ind_p_pb = np.where(
                dom.destinations[dest].excluded_id[II[ind_p], JJ[ind_p]] > 0
            )[0]
            ind3 = np.concatenate((ind3, ind_p[ind_p_pb]))

        if ind3.shape[0] > 0:
            if log:
                _LOGGER.info(
                    "check_people_in_box: relocating %d pedestrians located in "
                    "destination exclusion zones",
                    ind3.shape[0],
                )
            info = True
            xyrv[ind3, 0] = rng.uniform(
                xmin + p_rmax, xmax - p_rmax, ind3.shape[0]
            )
            xyrv[ind3, 1] = rng.uniform(
                ymin + p_rmax, ymax - p_rmax, ind3.shape[0]
            )
            continue

        if log:
            _LOGGER.info("check_people_in_box: all positions are admissible")
        break

    return info, xyrv


def remove_overlaps_in_box(dom, box, xyrv, dest_names, mass, dt, rng,
                           projection_method="cvxopt", dmin_people=0,
                           dmin_walls=0, itermax=10, log=True):
    """Iteratively remove pedestrian overlaps inside a bounding box.

    Repeated velocity projections are used to separate overlapping pedestrians.
    After each projection, positions are checked and invalid positions are
    resampled inside the requested box.

    Parameters
    ----------
    dom : Domain
        Domain object.
    box : sequence of float
        Bounding box ``[xmin, xmax, ymin, ymax]``.
    xyrv : numpy.ndarray, shape (n_people, 4)
        Pedestrian states. Positions are modified in place.
    dest_names : sequence of str
        Destination name for each pedestrian.
    mass : numpy.ndarray, shape (n_people,)
        Pedestrian masses.
    dt : float
        Time step used by the projection.
    rng : numpy.random.RandomState
        Random-number generator.
    projection_method : str, optional
        Projection solver passed to :func:`projection`.
    dmin_people : float, optional
        Minimum pedestrian-pedestrian distance.
    dmin_walls : float, optional
        Minimum pedestrian-wall distance.
    itermax : int, optional
        Maximum number of outer overlap-removal iterations.
    log : bool, optional
        If ``True``, emit progress diagnostics.

    Returns
    -------
    numpy.ndarray
        Updated pedestrian states.
    """
    if log:
        _LOGGER.info(
            "remove_overlaps_in_box: starting overlap removal for %d pedestrians",
            xyrv.shape[0],
        )

    dmax_init = xyrv[:, 2].max()
    iteration = 0

    while True:
        if log:
            _LOGGER.info(
                "remove_overlaps_in_box: box=%s, iteration=%d/%d",
                box,
                iteration,
                itermax,
            )

        contacts = compute_contacts(dom, xyrv, dmax_init)
        _, _, desired_velocity = dom.people_desired_velocity(xyrv, dest_names)
        info_projection, _, velocity, _, _, _ = projection(
            dt,
            xyrv,
            contacts,
            desired_velocity,
            mass,
            method=projection_method,
            dmin_people=dmin_people,
            dmin_walls=dmin_walls,
            nb_iter_max=10000,
            log=False,
        )

        xyrv[:, 0] += dt * velocity[:, 0]
        xyrv[:, 1] += dt * velocity[:, 1]

        iteration += 1
        positions_changed, xyrv = check_people_in_box(
            dom, box, xyrv, dest_names, rng, log=log
        )

        if info_projection >= 0 and positions_changed is False:
            if log:
                _LOGGER.info(
                    "remove_overlaps_in_box: completed successfully after %d iterations",
                    iteration,
                )
            break

        if iteration >= itermax:
            _LOGGER.warning(
                "remove_overlaps_in_box: initialization stopped after %d iterations; "
                "some overlaps or invalid positions may remain",
                itermax,
            )
            break

    return xyrv


def _create_people_for_group(
    nn,
    group,
    dom,
    dt,
    rng,
    dmin_people,
    dmin_walls,
    itermax,
    projection_method,
    log,
):
    """Create one population group using its box or RGB-position definition."""
    box = group["box"]
    color = group["color_for_positions"]

    if len(box) == 4 and len(color) == 0:
        if log:
            _LOGGER.info(
                "population initialization: count=%d, placement=box, box=%s",
                nn,
                box,
            )
        people, destinations, viewing_angles, masses = create_people_in_box(
            nn,
            box,
            group["destination"],
            group["radius_distribution"],
            group["velocity_distribution"],
            group["viewing_angle_distribution"],
            group["mass_distribution"],
            dom,
            rng,
            log=log,
        )
        people = remove_overlaps_in_box(
            dom,
            box,
            people,
            destinations,
            masses,
            dt,
            rng,
            dmin_people=dmin_people,
            dmin_walls=dmin_walls,
            itermax=itermax,
            projection_method=projection_method,
            log=log,
        )
        return people, destinations, viewing_angles, masses

    if len(box) == 0 and len(color) == 3:
        if log:
            _LOGGER.info(
                "population initialization: count=%d, placement=rgb, color=%s",
                nn,
                color,
            )
        return create_people_at_positions(
            nn,
            color,
            group["destination"],
            group["radius_distribution"],
            group["velocity_distribution"],
            group["viewing_angle_distribution"],
            group["mass_distribution"],
            dom,
            rng,
            log=log,
        )

    raise ValueError(
        "A population group must define either a four-value 'box' or a "
        "three-value 'color_for_positions', but not both."
    )


def add_people_flows(dom, groups, dt, dmin_people=0, dmin_walls=0, rng=None,
                     itermax=10, projection_method="cvxopt",
                     flow_duration=0.0, seed=0, frustration_memory_size=10,
                     log=True):
    """Create pedestrians generated by time-dependent population flows.

    Each group marked with ``"flux_per_minute"`` accumulates a fractional flow.
    Whenever the accumulated value exceeds one pedestrian, the corresponding
    integer number of pedestrians is created and the fractional remainder is
    kept in ``groups[...]["flux"]`` for the next call.

    Parameters
    ----------
    dom : Domain
        Domain in which pedestrians are created.
    groups : sequence of dict
        Population-group definitions.
    dt : float
        Simulation time step.
    dmin_people : float, optional
        Minimum pedestrian-pedestrian distance during initialization.
    dmin_walls : float, optional
        Minimum pedestrian-wall distance during initialization.
    rng : numpy.random.RandomState or None, optional
        Random-number generator. A new generator is created when ``None``.
    itermax : int, optional
        Maximum number of overlap-removal iterations.
    projection_method : str, optional
        Projection solver used while removing overlaps.
    flow_duration : float, optional
        Duration, in seconds, represented by the current flow update.
    seed : int, optional
        Seed used only when a new random generator is created. ``0`` requests a
        non-deterministic seed.
    frustration_memory_size : int, optional
        Number of entries allocated for frustration history.
    log : bool, optional
        If ``True``, emit flow-initialization diagnostics.

    Returns
    -------
    dict or None
        Dictionary containing newly created pedestrians, or ``None`` when no
        pedestrian must be added during this call.
    """
    if log:
        _LOGGER.info("add_people_flows: updating population flows")

    if rng is None:
        if log:
            _LOGGER.info("add_people_flows: creating a new RandomState")
        if seed == 0:
            rng = np.random.RandomState()
        else:
            rng = np.random.RandomState(seed)
        if log:
            _LOGGER.info(
                "add_people_flows: RandomState seed=%d",
                rng.get_state()[1][0],
            )

    # Convert accumulated flow rates into integer numbers of pedestrians while
    # retaining the fractional remainder for the next call.
    group_counts = []
    for group in groups:
        if group["flux_per_minute"] is True:
            if "flux" not in group:
                group["flux"] = group["nb"] * flow_duration / 60.0
            else:
                group["flux"] += group["nb"] * flow_duration / 60.0

            if group["flux"] >= 1:
                group_counts.append(int(group["flux"]))
                group["flux"] -= group_counts[-1]
            else:
                group_counts.append(0)
        else:
            group["flux"] = 0.0
            group_counts.append(0)

    if log:
        _LOGGER.info(
            "add_people_flows: group_counts=%s, remaining_flows=%s",
            group_counts,
            [group["flux"] for group in groups],
        )

    n_people = sum(group_counts)
    if n_people == 0:
        return None

    xyrv = np.zeros((n_people, 4))
    destinations = []
    viewing_angles = np.zeros((n_people,))
    masses = np.zeros((n_people,))
    with_frustration = np.zeros((n_people,))
    pos = 0

    for group, nn in zip(groups, group_counts):
        if nn <= 0:
            continue
        try:
            pp, pp_dest, pp_view, pp_mass = _create_people_for_group(
                nn,
                group,
                dom,
                dt,
                rng,
                dmin_people,
                dmin_walls,
                itermax,
                projection_method,
                log,
            )
        except Exception as exc:
            raise ValueError(f"The population group {group!r} cannot be created") from exc

        xyrv[pos:pos + nn, :] = pp
        destinations += pp_dest
        viewing_angles[pos:pos + nn] = pp_view
        masses[pos:pos + nn] = pp_mass
        if group["with_frustration"] is True:
            with_frustration[pos:pos + nn] = 1
        pos += nn

    people_to_add = {
        "xyrv": xyrv,
        "Vd": np.zeros((n_people, 2)),
        "U": np.zeros((n_people, 2)),
        "Uold": np.zeros((n_people, 2)),
        "viewing_angle": viewing_angles,
        "mass": masses,
        "destinations": np.array(destinations, dtype="U20"),
        "with_frustration": with_frustration,
        "frustration": np.zeros((n_people,)),
        "frustration_last_pos": -np.ones((n_people,), dtype=int),
        "frustration_memory": -np.ones((n_people,frustration_memory_size)),
        "frustration_memory_times": -np.ones((n_people,frustration_memory_size)),
        "beta": np.ones((n_people,)),
        "paths": {},
        "rng": rng,
    }
    people_to_add["id"] = np.char.add(
        [dom.name + "_"] * n_people,
        np.arange(dom.people_counter, dom.people_counter + n_people).astype("<U6"),
    )
    for ip, person_id in enumerate(people_to_add["id"]):
        people_to_add["paths"][person_id] = people_to_add["xyrv"][ip, :2]

    dom.people_counter += n_people
    return people_to_add


def people_initialization(dom, groups, dt, dmin_people=0, dmin_walls=0, rng=None,
                          seed=0, itermax=10, frustration_memory_size=10,
                          flow_duration=0.0, projection_method="cvxopt", log=True):
    """Initialize the pedestrian population for one domain.

    Flow-driven groups are processed first, followed by groups with a fixed
    initial population. Pedestrians are assigned identifiers and auxiliary
    arrays used by the microscopic model.

    Parameters
    ----------
    dom : Domain
        Domain in which pedestrians are initialized.
    groups : sequence of dict
        Population-group definitions. Each group may define ``"nb"``,
        ``"flux_per_minute"``, distribution specifications, ``"box"``,
        ``"color_for_positions"``, ``"destination"`` and
        ``"with_frustration"``.
    dt : float
        Simulation time step.
    dmin_people : float, optional
        Minimum pedestrian-pedestrian distance.
    dmin_walls : float, optional
        Minimum pedestrian-wall distance.
    rng : numpy.random.RandomState or None, optional
        Random-number generator.
    seed : int, optional
        Seed used only when ``rng`` is ``None``. ``0`` requests a
        non-deterministic seed.
    itermax : int, optional
        Maximum number of overlap-removal iterations.
    frustration_memory_size : int, optional
        Number of entries allocated for frustration history.
    flow_duration : float, optional
        Duration, in seconds, represented by the current flow update.
    projection_method : str, optional
        Projection solver used while removing overlaps.
    log : bool, optional
        If ``True``, emit initialization diagnostics.

    Returns
    -------
    dict or None
        Initialized pedestrian dictionary, or ``None`` if no pedestrian is
        created.
    """
    if log:
        _LOGGER.info("people_initialization: starting population initialization")

    if rng is None:
        if log:
            _LOGGER.info("people_initialization: creating a new RandomState")
        if seed == 0:
            rng = np.random.RandomState()
        else:
            rng = np.random.RandomState(seed)
        if log:
            _LOGGER.info(
                "people_initialization: RandomState seed=%d",
                rng.get_state()[1][0],
            )

    # Flow groups are initialized first. The original solver/default iteration
    # settings used by add_people_flows are intentionally preserved here.
    people = add_people_flows(
        dom,
        groups,
        dt,
        dmin_people=dmin_people,
        dmin_walls=dmin_walls,
        rng=rng,
        flow_duration=flow_duration,
        frustration_memory_size=frustration_memory_size,
        log=log,
    )
    n_flow_people = 0 if people is None else people["xyrv"].shape[0]
    if log:
        _LOGGER.info(
            "people_initialization: created %d pedestrians from flows",
            n_flow_people,
        )

    no_flow_groups = [group for group in groups if group["flux_per_minute"] is False]
    n_initial = sum(group["nb"] for group in no_flow_groups)

    if n_initial > 0:
        if people is not None:
            pos = people["xyrv"].shape[0]
            xyrv = np.concatenate((people["xyrv"], np.zeros((n_initial, 4))))
            destinations = list(people["destinations"])
            ids = people["id"]
            viewing_angles = np.concatenate(
                (people["viewing_angle"], np.zeros((n_initial,)))
            )
            masses = np.concatenate((people["mass"], np.zeros((n_initial,))))
            with_frustration = np.concatenate(
                (people["with_frustration"], np.zeros((n_initial,)))
            )
        else:
            xyrv = np.zeros((n_initial, 4))
            destinations = []
            ids = np.array([])
            viewing_angles = np.zeros((n_initial,))
            masses = np.zeros((n_initial,))
            with_frustration = np.zeros((n_initial,))
            pos = 0

        for group in no_flow_groups:
            nn = group["nb"]
            if nn <= 0:
                continue
            if log:
                _LOGGER.info(
                    "people_initialization: creating fixed group with count=%d",
                    nn,
                )
            try:
                pp, pp_dest, pp_view, pp_mass = _create_people_for_group(
                    nn,
                    group,
                    dom,
                    dt,
                    rng,
                    dmin_people,
                    dmin_walls,
                    itermax,
                    projection_method,
                    log,
                )
            except Exception as exc:
                raise ValueError(
                    f"The population group {group!r} cannot be created"
                ) from exc

            xyrv[pos:pos + nn, :] = pp
            destinations += pp_dest
            ids = np.concatenate((ids,
                np.char.add(
                    [dom.name + "_"] * nn,
                    np.arange(dom.people_counter, dom.people_counter + nn).astype("<U6"),
                )
            ))
            dom.people_counter += nn
            viewing_angles[pos:pos + nn] = pp_view
            masses[pos:pos + nn] = pp_mass
            if group["with_frustration"] is True:
                with_frustration[pos:pos + nn] = 1
            pos += nn

        if log:
            _LOGGER.info(
                "people_initialization: removing overlaps across the full domain"
            )
        xyrv = remove_overlaps_in_box(
            dom,
            [dom.xmin, dom.xmax, dom.ymin, dom.ymax],
            xyrv,
            destinations,
            masses,
            dt,
            rng,
            dmin_people=dmin_people,
            dmin_walls=dmin_walls,
            projection_method=projection_method,
            itermax=itermax,
            log=log,
        )

        if log:
            _LOGGER.info(
                "people_initialization: created %d pedestrians from fixed groups",
                n_initial,
            )

        people = {
            "xyrv": xyrv,
            "Vd": np.zeros((xyrv.shape[0], 2)),
            "U": np.zeros((xyrv.shape[0], 2)),
            "Uold": np.zeros((xyrv.shape[0], 2)),
            "viewing_angle": viewing_angles,
            "mass": masses,
            "destinations": np.array(destinations, dtype="U20"),
            "with_frustration": with_frustration,
            "frustration": np.zeros((xyrv.shape[0],)),
            "frustration_last_pos": -np.ones((xyrv.shape[0],), dtype=int),
            "frustration_memory": -np.ones((xyrv.shape[0],frustration_memory_size)),
            "frustration_memory_times": -np.ones((xyrv.shape[0],frustration_memory_size)),
            "beta": np.ones((xyrv.shape[0],)),
            "paths": {},
            "id": ids,
            "rng": rng,
        }
        for ip, person_id in enumerate(people["id"]):
            people["paths"][person_id] = people["xyrv"][ip, :2]

    if people is not None and log:
        _LOGGER.info(
            "people_initialization: total pedestrians=%d",
            people["xyrv"].shape[0],
        )

    return people


def find_duplicate_people(all_people, domains):
    """Build virtual pedestrians for overlapping transit regions.

    Pedestrians located inside a destination transit box are copied to the
    connected domain using the configured spatial shift. The reverse copy is
    also performed for pedestrians located in the shifted transit box.

    Parameters
    ----------
    all_people : dict
        Mapping from domain names to pedestrian dictionaries.
    domains : dict
        Mapping from domain names to :class:`Domain` objects.

    Returns
    -------
    dict
        Virtual pedestrian dictionaries for every domain.
    """
    virtual_people = {}
    for name in domains:
        virtual_people[name] = {}
        for key in PEOPLE_ARRAY_KEYS:
            virtual_people[name][key] = None
        for key in PEOPLE_LIST_KEYS:
            virtual_people[name][key] = []
        for key in PEOPLE_DICT_KEYS:        
            virtual_people[name][key] = None
        
    for name in domains:
        # print("===> Look at virtual people from domain ",name)
        dom = domains[name]
        for dest_name in dom.destinations:
            # _LOGGER.info(
            #     f"find_duplicate_people: domain={name}, "
            #     f"current destination={dest_name}, "
            #     f"all upcoming available destinations={dom.destinations[dest_name].next_destination}"
            # )
            for inextdest in range(len(dom.destinations[dest_name].next_destination)):
                box = dom.destinations[dest_name].next_transit_box[inextdest]
                box_xshift = dom.destinations[dest_name].next_transit_box_xshift[inextdest]
                box_yshift = dom.destinations[dest_name].next_transit_box_yshift[inextdest]
                if box is not None:
                    box = np.array(box)
                    box = np.vstack((box,box[0,:]))
                    
                    next_domain_name = dom.destinations[dest_name].next_domain[inextdest]

                    ## People in this domain "name" that must be replicated in
                    ## the domain "next_domain_name"
                    p = all_people[name]["xyrv"]
                    if (p.shape[0]>0):
                        #ind = np.where(all_people[name]["destinations"]==dest_name)[0]
                        #test = np.zeros(ind.shape[0])
                        test = np.zeros(p.shape[0])
                        for i in np.arange(box.shape[0]-1):
                            # vx = box[i,0]-p[ind,0]
                            # vy = box[i,1]-p[ind,1]
                            # wx = box[i+1,0]-p[ind,0]
                            # wy = box[i+1,1]-p[ind,1]
                            vx = box[i,0]-p[:,0]
                            vy = box[i,1]-p[:,1]
                            wx = box[i+1,0]-p[:,0]
                            wy = box[i+1,1]-p[:,1]
                            test += np.arccos(np.maximum(-1.0, np.minimum(1.0, ( (vx*wx+vy*wy)/np.sqrt(vx*vx+vy*vy)/np.sqrt(wx*wx+wy*wy)))))
                        #inside = ind[np.where(np.abs(test-2*np.pi)<0.001)[0]]
                        inside = np.where(np.abs(test-2*np.pi)<0.001)[0]
                        if (inside.shape[0]>0):
                            # _LOGGER.info(
                            #     f"find_duplicate_people: domain={name}, "
                            #     f"current destination={dest_name}, "
                            #     f"next destination={dom.destinations[dest_name].next_destination[inextdest]}, "
                            #     f"next domain={dom.destinations[dest_name].next_domain[inextdest]}, "
                            #     f"{inside.shape[0]} people inside transit box {box}, "
                            #     f"inside = {inside}, "
                            #     f"apply xhift {box_xshift} and yshift {box_yshift} "
                            #     f"=> put them in virtual_people"
                            # )
                            xyrv_inside = all_people[name]["xyrv"][inside].copy()
                            xyrv_inside[:,0] += box_xshift
                            xyrv_inside[:,1] += box_yshift
                    
                            if (virtual_people[next_domain_name]["xyrv"] is not None) \
                                and (virtual_people[next_domain_name]["xyrv"].shape[0]>0):
                                    
                                virtual_people[next_domain_name]["xyrv"] = np.concatenate((
                                    virtual_people[next_domain_name]["xyrv"],
                                    xyrv_inside))
                                
                                for key in [k for k in PEOPLE_ARRAY_KEYS if k != "xyrv"]:
                                    virtual_people[next_domain_name][key] = np.concatenate((
                                        virtual_people[next_domain_name][key],
                                        all_people[name][key][inside])
                                    )          
                                                   
                                for key in [k for k in PEOPLE_LIST_KEYS]:
                                    virtual_people[next_domain_name][key] = \
                                        virtual_people[next_domain_name][key] + \
                                        all_people[name][key][inside]
                            
                            else:
                                
                                virtual_people[next_domain_name]["xyrv"] \
                                    = xyrv_inside
                                    
                                for key in PEOPLE_ARRAY_KEYS:
                                    if key=="xyrv":
                                        continue
                                    else:
                                        virtual_people[next_domain_name][key] \
                                        = all_people[name][key][inside]
                                    
                                for key in [k for k in PEOPLE_LIST_KEYS]:
                                    virtual_people[next_domain_name][key] = \
                                        all_people[name][key][inside]

                            # _LOGGER.info(
                            #     f"find_duplicate_people: virtual_people[{next_domain_name}]="
                            #     f"{virtual_people[next_domain_name]}"
                            # )
                        ## To check if an individual appears two times...
                        if (inside.shape[0]>0):
                            virtual_people[next_domain_name]["xyrv"], indices \
                                = np.unique(virtual_people[next_domain_name]["xyrv"],
                                            return_index=True, axis=0)
                            for key in PEOPLE_ARRAY_KEYS:
                                if key=="xyrv":
                                    continue
                                else:
                                    virtual_people[next_domain_name][key] = \
                                        virtual_people[next_domain_name][key][indices]
                            
                    ## People in the domain "next_domain_name" that must be replicated
                    ## in this domain "name"
                    p = all_people[next_domain_name]["xyrv"]
                    if (p.shape[0]>0):
                        test = np.zeros(p.shape[0])
                        next_domain_box = box.copy()
                        next_domain_box[:,0] += box_xshift
                        next_domain_box[:,1] += box_yshift
                        # print("next_domain_box =", next_domain_box)

                        for i in np.arange(next_domain_box.shape[0]-1):
                            vx = next_domain_box[i,0]-p[:,0]
                            vy = next_domain_box[i,1]-p[:,1]
                            wx = next_domain_box[i+1,0]-p[:,0]
                            wy = next_domain_box[i+1,1]-p[:,1]
                            test += np.arccos( np.maximum(-1.0, np.minimum(1.0, ( (vx*wx+vy*wy)/np.sqrt(vx*vx+vy*vy)/np.sqrt(wx*wx+wy*wy)))))
                        #print("test=",test)
                        inside = np.where(np.abs(test-2*np.pi)<0.001)[0]
                        # print("people inside: ",inside,
                        #       " must be copy in domain ",name)
                        next_domain_xyrv_inside = all_people[next_domain_name]["xyrv"][inside].copy()
                        next_domain_xyrv_inside[:,0] -= box_xshift
                        next_domain_xyrv_inside[:,1] -= box_yshift
                        # print("next_domain box_xshift =",box_xshift," box_yshift =",box_yshift," next_domain_xyrv_inside =",next_domain_xyrv_inside)
                        try:
                            virtual_people[name]["xyrv"] = np.concatenate((
                                virtual_people[name]["xyrv"],
                                next_domain_xyrv_inside))
                        except Exception:
                            virtual_people[name]["xyrv"] \
                                = next_domain_xyrv_inside
                            for key in PEOPLE_ARRAY_KEYS:
                                if key=="xyrv":
                                    continue
                                else:
                                    try:
                                        virtual_people[name][key] = np.concatenate((
                                            virtual_people[name][key],
                                            all_people[next_domain_name][key][inside]))
                                    except Exception:
                                        virtual_people[name][key] \
                                            = all_people[next_domain_name][key][inside]

                        ## To check if an individual appears two times...
                        if (inside.shape[0]>0):
                            virtual_people[name]["xyrv"], indices \
                                = np.unique(virtual_people[name]["xyrv"],
                                            return_index=True, axis=0)
                            for key in PEOPLE_ARRAY_KEYS:
                                if key=="xyrv":
                                    continue
                                else:
                                    virtual_people[name][key] \
                                        = virtual_people[name][key][indices]

    return virtual_people


def append_people(people, people_to_add):
    """Append a newly generated pedestrian batch to an existing population."""
    for key in PEOPLE_ARRAY_KEYS:
        people[key] = np.concatenate(
            (people[key], people_to_add[key])
        )
    
    for key in PEOPLE_LIST_KEYS:
        people[key] = people[key] + people_to_add[key]
        
    for key in PEOPLE_DICT_KEYS:
        people[key] = people[key] | people_to_add[key]

    
def people_update_destination(all_people, domains, thld, box=None):
    """Update pedestrian destinations, domain transfers and removals.

    Pedestrians outside the allowed geometry are removed. Pedestrians close to
    their current destination are either transferred to a next domain, shifted
    through a transit box, assigned a new destination, or removed when no next
    destination is defined.

    Parameters
    ----------
    all_people : dict
        Mapping from domain names to pedestrian dictionaries.
    domains : dict
        Mapping from domain names to :class:`Domain` objects.
    thld : float
        Distance threshold used to decide whether a destination has been
        reached.
    box : sequence of float or None, optional
        Optional removal box ``[xmin, xmax, ymin, ymax]``. If ``None``, each
        domain bounding box is used.

    Returns
    -------
    dict
        Updated ``all_people`` mapping.
    """

    for name in domains:
        p = all_people[name]["xyrv"]
        if (p.shape[0]>0):
            dom = domains[name]
            # Test if there are people outside the domain or outside the given box
            # If it is the case, we remove these persons...
            if box is None:
                # Remove people who are outside the domain
                S = (p[:, 0]-p[:, 2] <= dom.xmin+dom.pixel_size) +\
                    (p[:, 0]-p[:, 2] >= dom.xmax-dom.pixel_size) +\
                    (p[:, 1]-p[:, 2] <= dom.ymin+dom.pixel_size) +\
                    (p[:, 1]-p[:, 2] >= dom.ymax-dom.pixel_size)
            else:
                # Remove people who are outside the given box
                S = (p[:, 0]-p[:, 2] <= box[0]+dom.pixel_size) +\
                    (p[:, 0]-p[:, 2] >= box[1]-dom.pixel_size) +\
                    (p[:, 1]-p[:, 2] <= box[2]+dom.pixel_size) +\
                    (p[:, 1]-p[:, 2] >= box[3]-dom.pixel_size)
            out = np.where(S)[0]
            if (out.shape[0] > 0):
                _LOGGER.info(
                    f"people_update_destination: removing {out.shape[0]} pedestrians "
                    f"outside domain={name} xy = {p[out,:2]}"
                )
                for ik, key in enumerate(all_people[name]):
                    # print(ik, key)
                    if (key != "paths"):
                        try:
                            if (all_people[name][key].shape[0] == p.shape[0]):
                                all_people[name][key] = \
                                    np.delete(all_people[name][key], out, 0)
                        except Exception:
                            pass
            
    for idom,name in enumerate(domains):
        p = all_people[name]["xyrv"]
        if (p.shape[0]>0):
            dom = domains[name]
            destinations = all_people[name]["destinations"]
            # Test if people are close to their destinations
            II = np.floor((p[:, 1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
            JJ = np.floor((p[:, 0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
            all_people[name]["II"] = II
            all_people[name]["JJ"] = JJ

            # Remove people with II or JJ outside
            testI = (II >= 0)*(II < dom.height)
            testJ = (JJ >= 0)*(JJ < dom.width)
            outI = np.where(testI == False)[0]
            outJ = np.where(testJ == False)[0]
            out = np.unique(np.concatenate((outI, outJ)))
            if (out.shape[0] > 0):
                _LOGGER.info(
                    "people_update_destination: deleting %d pedestrians outside "
                    "domain=%s; xy=%s, II=%s, JJ=%s. Check whether the "
                    "domain is overcrowded.",
                    out.shape[0],
                    name,
                    all_people[name]["xyrv"][out, :2],
                    II[out],
                    JJ[out],
                )
                for key in all_people[name]:
                    try:
                        if (all_people[name][key].shape[0] == p.shape[0]):
                            all_people[name][key] = np.delete(all_people[name][key], out, 0)
                    except Exception:
                        pass
                destinations = all_people[name]["destinations"]
                p = all_people[name]["xyrv"]
                II = np.floor((p[:, 1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
                JJ = np.floor((p[:, 0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
                all_people[name]["II"] = II
                all_people[name]["JJ"] = JJ

            # Remove people inside the walls
            in_wall = np.where(dom.wall_mask[II, JJ] > 0)[0]
            if (in_wall.shape[0] > 0):
                _LOGGER.info(
                    "people_update_destination: deleting %d pedestrians inside walls "
                    "in domain=%s; positions=%s, II=%s, JJ=%s. Check whether the "
                    "domain is overcrowded.",
                    in_wall.shape[0],
                    name,
                    all_people[name]["xyrv"][in_wall, :2],
                    II[in_wall],
                    JJ[in_wall],
                )
                for key in all_people[name]:
                    try:
                        if (all_people[name][key].shape[0] == p.shape[0]):
                            all_people[name][key] = \
                                np.delete(all_people[name][key], in_wall, 0)
                    except Exception:
                        pass
                destinations = all_people[name]["destinations"]
                p = all_people[name]["xyrv"]
                II = np.floor((p[:, 1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
                JJ = np.floor((p[:, 0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
                all_people[name]["II"] = II
                all_people[name]["JJ"] = JJ

            # Keep all people paths
            for ip, pid in enumerate(all_people[name]["id"]):
                try:
                    all_people[name]["paths"][pid] = \
                        np.hstack((all_people[name]["paths"][pid], p[ip, :2]))
                except Exception:
                    all_people[name]["paths"][pid] = p[ip, :2]
            # print("path: ",all_people[name]["paths"])

            # Change destination or remove people who are at their final destination
            id_rm = None
            for id,dest_name in enumerate(np.unique(destinations)):

                ind = np.where(np.array(destinations) == dest_name)[0]
                
                D = dom.destinations[dest_name].distance[II[ind], JJ[ind]]-p[ind, 2]
                # Test if people are close to their destination
                ind2 = np.where(D <= thld)[0]
                # print("     ",ind2," are close to their destination ",dest_name)
                if (ind2.shape[0] > 0):
                    # weight renormalization
                    weight_sum = np.sum(dom.destinations[dest_name].next_destination_weights)
                    if weight_sum<=0:
                        raise ValueError(f"[micro people_update_destination] : !!! ERROR !!! \
                            weights are not correct ! current destination = {dest_name}, \
                            next destinations = {dom.destinations[dest_name].next_destination}, \
                            next_destination_weights = {dom.destinations[dest_name].next_destination_weights} \
                            ")
                    else:
                        next_destination_normalized_weights = dom.destinations[dest_name].next_destination_weights/weight_sum
                        ndest = len(dom.destinations[dest_name].next_destination)
                        rng = all_people[name]["rng"]
                        # if rng is None:
                        #     rng = np.random.RandomState()
                        #     print("** WARNING: NEW RandomState with seed = 0 **")
                        next_destination_id = rng.choice(
                            np.arange(ndest), 
                            size=ind2.shape[0], 
                            replace=True, p=next_destination_normalized_weights)
                        next_destinations = [dom.destinations[dest_name].next_destination[i] for i in next_destination_id]
                        next_domains = [dom.destinations[dest_name].next_domain[i] for i in next_destination_id]
                        
                        for inextdest in range(ndest):
                            ind3 = np.where(next_destination_id==inextdest)[0]
                            nextdest = dom.destinations[dest_name].next_destination[inextdest]
                            if (nextdest is not None):
                                next_domain = dom.destinations[dest_name].next_domain[inextdest]
                                if (next_domain != name):
                                    # In this case we have to copy these persons in their new
                                    # domain, and then remove them
                                    # update destination
                                    all_people[name]["destinations"][ind[ind2[ind3]]] \
                                        = nextdest
                                    # move
                                    if (all_people[next_domain]["xyrv"].shape[0]>0):
                                        for key in PEOPLE_ARRAY_KEYS:
                                            all_people[next_domain][key] = np.concatenate(
                                                (all_people[next_domain][key],
                                                all_people[name][key][ind[ind2[ind3]]]))
                                    else:
                                        for key in PEOPLE_ARRAY_KEYS:
                                            all_people[next_domain][key] = all_people[name][key][ind[ind2[ind3]]]
                                            
                                    # remove
                                    try:
                                        id_rm = np.concatenate((id_rm,ind[ind2[ind3]]))
                                    except Exception:
                                        id_rm = ind[ind2[ind3]]

                                elif (dom.destinations[dest_name].next_transit_box[inextdest] is not None):
                                    # next_domain == domain BUT a transit box + a shift are defined
                                    # print("people before shift : ", all_people[name]["xyrv"][ind[ind2]])
                                    all_people[name]["xyrv"][ind[ind2[ind3]],0] += dom.destinations[dest_name].next_transit_box_xshift[inextdest]
                                    all_people[name]["xyrv"][ind[ind2[ind3]],1] += dom.destinations[dest_name].next_transit_box_yshift[inextdest]
                                    # print("people after shift : ", all_people[name]["xyrv"][ind[ind2]])
                                    all_people[name]["destinations"][ind[ind2[ind3]]] \
                                        = nextdest       
                                else:
                                    # People do not change domains, they juste change their
                                    # destination
                                    #print("     => change destination...")
                                    _LOGGER.info(
                                        "people_update_destination: people who don't change domains, "
                                        f"just destinations => ids={all_people[name]["id"][ind[ind2[ind3]]]}, "
                                        f"xy={all_people[name]["id"][ind[ind2[ind3]]]}"
                                    )
                                    all_people[name]["destinations"][ind[ind2[ind3]]] \
                                        = nextdest    
                            else: ## in this case, we will only remove these persons...
                                #print("     => remove...")
                                try:
                                    id_rm = np.concatenate((id_rm,ind[ind2[ind3]]))
                                except Exception:
                                    id_rm = ind[ind2[ind3]]

            if (id_rm is not None):
                for key in all_people[name]:
                    try:
                        if (all_people[name][key].shape[0] == p.shape[0]):
                            all_people[name][key] \
                                = np.delete(all_people[name][key], id_rm, 0)
                    except Exception as e:
                        # print(f"!! warning !! Failed to reduce '{key}', error = {e}")
                        pass

    return all_people


def sensor(door, xy0, xy1, t0, t1, dist_in, dist_out):
    """Detect pedestrian crossings of a line segment sensor.

    The pedestrian trajectory during ``[t0, t1]`` is approximated by a segment.
    The function computes intersections with the sensor segment, crossing
    directions, approximate crossing times and pedestrians currently located
    near either side of the sensor.

    Parameters
    ----------
    door : array-like, shape (4,)
        Sensor segment coordinates ``[x0, y0, x1, y1]``.
    xy0 : numpy.ndarray, shape (n_people, 2)
        Pedestrian positions at time ``t0``.
    xy1 : numpy.ndarray, shape (n_people, 2)
        Pedestrian positions at time ``t1``.
    t0 : float
        Initial time.
    t1 : float
        Final time.
    dist_in : float
        Maximum distance used to identify pedestrians on the entry side.
    dist_out : float
        Maximum distance used to identify pedestrians on the exit side.

    Returns
    -------
    id : numpy.ndarray
        Indices of pedestrians crossing the sensor.
    p : numpy.ndarray
        Intersection points for crossing pedestrians.
    io : numpy.ndarray
        Crossing directions: ``1`` for exits and ``-1`` for entries.
    times : numpy.ndarray
        Estimated crossing times.
    entries : int
        Number of entries.
    exits : int
        Number of exits.
    id_people_in : numpy.ndarray
        Indices of pedestrians near the entry side.
    id_people_out : numpy.ndarray
        Indices of pedestrians near the exit side.
    """
    #
    #   trajectories:
    #               xy0
    #                |       id_people_in    
    #                 |
    #   door:     d0--p--d1
    #                   |    id_people_out 
    #                   xy1
    d0 = np.empty(xy0.shape)
    d0[:, 0] = door[0]
    d0[:, 1] = door[1]
    d1 = np.empty(xy1.shape)
    d1[:, 0] = door[2]
    d1[:, 1] = door[3]
    T = np.array([[0, -1], [1, 0]])
    vdoor = np.atleast_2d(d1 - d0)
    vtraj = np.atleast_2d(xy1 - xy0)
    v0 = np.atleast_2d(d0 - xy0)
    dot_vdoor_T = np.dot(vdoor, T)
    denom = np.sum(dot_vdoor_T * vtraj, axis=1)
    num = np.sum(dot_vdoor_T * v0, axis=1)
    # Intersection points
    # can be inf or nan if parallel lines...
    p = np.atleast_2d(num / denom).T * vtraj + xy0
    # Test if the intersection point is on the door segment
    vp0 = np.atleast_2d(p - d0)
    norm_vdoor_2 = np.sum(vdoor*vdoor, axis=1)
    dot_vdoor_vp0 = np.sum(vdoor*vp0, axis=1)
    is_p_in_door = (dot_vdoor_vp0 >= 0)*(dot_vdoor_vp0 <= norm_vdoor_2)
    # Test if the intersection point is on the person trajectory
    vpxy0 = np.atleast_2d(p - xy0)
    norm_vtraj_2 = np.sum(vtraj*vtraj, axis=1)
    dot_vtraj_vpxy0 = np.sum(vtraj*vpxy0, axis=1)
    is_p_in_traj = (dot_vtraj_vpxy0 >= 0)*(dot_vtraj_vpxy0 <= norm_vtraj_2)
    # Keep only points on the door and on the trajectory
    is_p_intersect = is_p_in_door*is_p_in_traj
    id = np.where(is_p_intersect)[0]
    # Test if the direction is the output normal: (d1-d0)_y , (d1-d0)_x
    vn = np.empty(vdoor.shape)
    vn[:, 0] = vdoor[:, 1]
    vn[:, 1] = -vdoor[:, 0]
    dot_vn_vtraj = np.sum(vn*vtraj, axis=1)
    is_normal_dir = (dot_vn_vtraj > 0)
    io = (is_normal_dir[id])*1+(is_normal_dir[id] is False)*(-1)
    exits = np.sum(io == 1)
    entries = np.sum(io == -1)
    # Find people close to the sensor in the direction of the inputs (before the sensor)
    # and find people close to the sensor in the direction of the outputs (after the sensor)
    # print("door :",door)
    vtraj_to_sensor = np.atleast_2d(0.5*(d0+d1)-xy1)
    # print("vtraj_to_sensor :", vtraj_to_sensor)
    dist_to_sensor = np.linalg.norm(vtraj_to_sensor, axis=1)
    # print("dist_to_sensor :", dist_to_sensor)
    dot_vn_vtraj_to_sensor = np.sum(vn*vtraj_to_sensor, axis=1)
    # print("dot_vn_vtraj_to_sensor :", dot_vn_vtraj_to_sensor)
    all_io = (dot_vn_vtraj_to_sensor>0)*1+(dot_vn_vtraj_to_sensor<=0)*(-1)
    # print("all_io :", all_io)
    id_people_in  = np.where((all_io ==  1)*(dist_to_sensor<=dist_in))[0]
    id_people_out = np.where((all_io == -1)*(dist_to_sensor<=dist_out))[0]
    # Compute the distance from the intersection point p to xy0
    norm_vpxy0_2 = np.sqrt(np.sum(vpxy0 * vpxy0, axis=1))
    # Compute the distance from the intersection point p to xy1
    vpxy1 = np.atleast_2d(p - xy1)
    norm_vpxy1_2 = np.sqrt(np.sum(vpxy1 * vpxy1, axis=1))
    # Compute the intersection time
    norm_vtraj = np.sqrt(norm_vtraj_2)
    dt = t1-t0
    times = t0 + (is_normal_dir)*(norm_vpxy0_2*dt/norm_vtraj) + \
        (is_normal_dir is False)*(norm_vpxy1_2*dt/norm_vtraj)
    return id, p[id, :], io, times[id], entries, exits, id_people_in, id_people_out

