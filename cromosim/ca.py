# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
# License: GPL

import numpy as np
import random
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from matplotlib import colors
import matplotlib.cm as cm


def plot_people_according_to_current_door_distance(ifig, people, domain,
                                                   axis=None,
                                                   savefig=False,
                                                   filename='fig.png'):
    """
    To draw occupied cells with colors depending on current door distances

    Parameters
    ----------

    ifig: int
        figure number
    people: numpy masked arrays
        equal to 1 if the cell (i,j) is occupied, 0 elsewhere
    domain: Domain
        contains everything for managing the domain
    axis: list
        matplotlib axis
    savefig: boolean
        writes the figure as a png file if true
    filename: string
        png filename used to write the figure
    """
    fig = plt.figure(ifig)
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.imshow(domain.image, interpolation='nearest',
               extent=[domain.xmin, domain.xmax, domain.ymin, domain.ymax],
               origin='lower')
    ax1.imshow(domain.destinations["door"].distance * (people == 1), interpolation='nearest',
               extent=[domain.xmin, domain.xmax, domain.ymin, domain.ymax],
               origin='lower', alpha=1.0, cmap='Greys')
    if (axis):
        ax1.set_xlim(axis[0], axis[1])
        ax1.set_ylim(axis[2], axis[3])
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.axis('off')
    # ax1.set_xlabel('x',color='white')
    # ax1.set_ylabel('y',color='white')
    # ax1.set_title('t',color='white')
    # ax1.spines['left'].set_color('white')
    # ax1.spines['right'].set_color('white')
    # ax1.spines['bottom'].set_color('white')
    # ax1.spines['top'].set_color('white')
    fig.canvas.draw()
    if (savefig):
        fig.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=0)
    plt.pause(0.01)


def plot_people_according_to_initial_door_distance(ifig, people, domain, results,
                                                   axis=None, savefig=False,
                                                   filename='fig.png'):
    """
    To draw occupied cells with colors depending on initial (time=0) door \
    distances

    Parameters
    ----------

    ifig: int
        figure number
    people: numpy masked arrays
        equal to 1 if the cell (i,j) is occupied, 0 elsewhere
    domain: Domain
        contains everything for managing the domain
    results: numpy array
        contains the results for each iteration in time, used here to \
        determine the initial distance to the door
    axis: list
        matplotlib axis
    savefig: boolean
        writes the figure as a png file if true
    filename: string
        png filename used to write the figure
    """
    fig = plt.figure(ifig)
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.imshow(domain.image, interpolation='nearest',
               extent=[domain.xmin, domain.xmax, domain.ymin, domain.ymax],
               origin='lower')
    I0 = results[:, 0, 0]
    J0 = results[:, 1, 0]
    D0 = domain.destinations["door"].distance[I0, J0]
    II = results[:, 0, -1]
    JJ = results[:, 1, -1]
    tmp = domain.destinations["door"].distance.copy()
    tmp[II, JJ] = D0
    ax1.imshow(tmp * (people == 1), interpolation='nearest',
               extent=[domain.xmin, domain.xmax, domain.ymin, domain.ymax],
               origin='lower', alpha=1.0, cmap='Greys')
    if (axis):
        ax1.set_xlim(axis[0], axis[1])
        ax1.set_ylim(axis[2], axis[3])
    ax1.set_xticks([])
    ax1.set_yticks([])
    # ax1.set_xlabel('x',color='white')
    # ax1.set_ylabel('y',color='white')
    # ax1.set_title('t',color='white')
    # ax1.spines['left'].set_color('white')
    # ax1.spines['right'].set_color('white')
    # ax1.spines['bottom'].set_color('white')
    # ax1.spines['top'].set_color('white')
    ax1.axis('off')
    fig.canvas.draw()
    if (savefig):
        fig.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=0)
    plt.pause(0.01)


def compute_exit_times(dt, results):
    """
    To compute exit times for all the individuals from the results array

    Parameters
    ----------

    dt: float
        time step
    results: numpy array
        positions at each time

    Returns
    -------

    exit_times: all the exit times
    """
    Np = results.shape[0]
    exit_times = -np.ones(Np)
    for id in np.arange(Np):
        exit_times[id] = dt * np.where(results[id, 0, :] == -1)[0].min()
    return exit_times


def plot_people_according_to_exit_times(ifig, dt, people, domain, results,
                                        axis=None, savefig=False,
                                        filename='fig.png'):
    """
    To draw occupied cells with colors depending on the exit times

    Parameters
    ----------

    ifig: int
        figure number
    dt: float
        time step
    people: masked ndarray
        equal to 1 if the cell (i,j) is occupied, 0 elsewhere
    domain: Domain
        contains everything for managing the domain
    results: ndarray
        positions at each time
    axis: list
        matplotlib axis
    savefig: boolean
        writes the figure as a png file if true
    filename: string
        png filename used to write the figure
    """
    fig = plt.figure(ifig)
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.imshow(domain.image, interpolation='nearest',
               extent=[domain.xmin, domain.xmax, domain.ymin, domain.ymax],
               origin='lower')
    I0 = results[:, 0, 0]
    J0 = results[:, 1, 0]
    exit_times = compute_exit_times(dt, results)
    tmp = domain.destinations["door"].distance.copy()
    tmp[I0, J0] = exit_times
    people_init = np.ma.MaskedArray(np.zeros(people.shape, dtype=int),
                                    mask=domain.wall_mask)
    people_init.data[I0, J0] = 1
    ax1.imshow(tmp * (people_init == 1), interpolation='nearest',
               extent=[domain.xmin, domain.xmax, domain.ymin, domain.ymax],
               origin='lower', alpha=1.0, cmap='Greys')
    if (axis):
        ax1.set_xlim(axis[0], axis[1])
        ax1.set_ylim(axis[2], axis[3])
    ax1.set_xticks([])
    ax1.set_yticks([])
    # ax1.set_xlabel('x',color='white')
    # ax1.set_ylabel('y',color='white')
    # ax1.set_title('t',color='white')
    # ax1.spines['left'].set_color('white')
    # ax1.spines['right'].set_color('white')
    # ax1.spines['bottom'].set_color('white')
    # ax1.spines['top'].set_color('white')
    ax1.axis('off')
    fig.canvas.draw()
    if (savefig):
        fig.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=0)
    plt.pause(0.01)


def plot_people_paths(ifig, dt, pixel_size, people, domain, results, axis=None,
                      savefig=False, filename='fig.png'):
    """
    To draw all the individual paths from intial time to final time

    Parameters
    ----------

    ifig: int
        figure number
    dt: float
        time step
    pixel_size: float
        size of one pixel in meters
    people: numpy masked arrays
        equal to 1 if the cell (i,j) is occupied, 0 elsewhere
    domain: Domain
        contains everything for managing the domain
    results: numpy.ndarray
        positions at each time
    axis: list
        matplotlib axis
    savefig: boolean
        writes the figure as a png file if true
    filename: string
        png filename used to write the figure
    """
    Np = results.shape[0]
    fig = plt.figure(ifig)
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.imshow(domain.image, interpolation='nearest',
               extent=[domain.xmin, domain.xmax, domain.ymin, domain.ymax],
               origin='lower')
    I0 = results[:, 0, 0]
    J0 = results[:, 1, 0]
    exit_times = compute_exit_times(dt, results)
    tmp = domain.destinations["door"].distance.copy()
    tmp[I0, J0] = exit_times
    people_init = np.ma.MaskedArray(np.zeros(people.shape, dtype=int),
                                    mask=domain.wall_mask)
    people_init.data[I0, J0] = 1
    vmin = (tmp * (people_init == 1)).min()
    vmax = (tmp * (people_init == 1)).max()
    cmap = cm.Greys
    ax1.imshow(tmp * (people_init == 1), interpolation='nearest',
               extent=[domain.xmin, domain.xmax, domain.ymin, domain.ymax],
               origin='lower', alpha=1.0, cmap=cmap,
               norm=colors.BoundaryNorm(boundaries=np.linspace(vmin, vmax, cmap.N + 1),
                                        ncolors=cmap.N))
    mask = (results == -1)
    traj = np.swapaxes(np.ma.MaskedArray(results, dtype=int, mask=mask), 1, 2)
    traj = traj[:, :, [1, 0]]
    col = np.array([])
    for id in np.arange(Np):
        ns = np.where(results[id, 0, :] != -1)[0].shape[0] - 1
        col = np.hstack((col, np.ones(ns) * exit_times[id]))
    paths = LineCollection(0.5 * pixel_size + pixel_size * traj, linewidths=2,
                           linestyle="solid", cmap=cmap,
                           norm=colors.BoundaryNorm(
                               boundaries=np.linspace(vmin, vmax, cmap.N + 1),
                               ncolors=cmap.N)
                           )
    paths.set_array(exit_times)
    ax1.add_collection(paths)
    if (axis):
        ax1.set_xlim(axis[0], axis[1])
        ax1.set_ylim(axis[2], axis[3])
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.axis('off')
    # ax1.set_xlabel('x',color='white')
    # ax1.set_ylabel('y',color='white')
    # ax1.set_title('t',color='white')
    # ax1.spines['left'].set_color('white')
    # ax1.spines['right'].set_color('white')
    # ax1.spines['bottom'].set_color('white')
    # ax1.spines['top'].set_color('white')
    # ax1.axis('off')
    fig.canvas.draw()
    if (savefig):
        fig.savefig(filename, dpi=150, bbox_inches='tight', pad_inches=0)


def sequential_update(people, people_ij, weight, shuffle=None,
                      randomstate=None):
    """
    To move all individuals sequentially according to the following rule: \
    the update of individual i is determined in a stochastic way by computing \
    transition probabilities on neighboring cells (including the current \
    position of i).

    Parameters
    ----------

    people: numpy masked arrays
        equal to 1 if the cell (i,j) is occupied, 0 elsewhere
    people_ij: numpy array
        (i,j) for each individual
    weight: numpy array
        weights for the probabilities in order to move in such a way like to \
        reach the door
    shuffle: string
        shuffle kind ('random' or 'random_frozen'): if the sequential order \
        changes or not at each time
    randomstate: numpy randomstate
        create a new one or reuse the given random state

    Returns
    -------

    people: numpy masked arrays
        new positions equal to 1 if the cell (i,j) is occupied, 0 otherwise
    people_id: numpy masked arrays
        new people index (i,j)
    """
    Np = people_ij.shape[0]
    height, width = people.shape
    if (randomstate is None):
        randomstate = np.random.RandomState()
    # print("---------- sequential_update: shuffle = ",shuffle)
    # Default: no shuffle update
    order = np.arange(Np)
    if (shuffle == 'random'):
        # Random shuffle update
        randomstate.shuffle(order)
    elif (shuffle == 'random_frozen'):
        # Frozen shuffle update
        local_seed = 0.5
        random.shuffle(order, lambda: local_seed)
    for id in order:
        i = people_ij[id, 0]
        j = people_ij[id, 1]
        # print("--------- sequential_update: id = ",id," i = ",i," j = ",j)
        if ((i != -1) and (j != -1)):
            II = [i]
            JJ = [j]
            if (i >= 1):
                im = i - 1
                if ((people[im, j] is not np.ma.masked) and (people[im, j] != 1)):
                    II.append(im)
                    JJ.append(j)
            if (i <= height - 2):
                ip = i + 1
                if ((people[ip, j] is not np.ma.masked) and (people[ip, j] != 1)):
                    II.append(ip)
                    JJ.append(j)
            if (j >= 1):
                jm = j - 1
                if ((people[i, jm] is not np.ma.masked) and (people[i, jm] != 1)):
                    II.append(i)
                    JJ.append(jm)
            if (j <= width - 2):
                jp = j + 1
                if ((people[i, jp] is not np.ma.masked) and (people[i, jp] != 1)):
                    II.append(i)
                    JJ.append(jp)
            w = weight[II, JJ] / np.sum(weight[II, JJ])
            # print("--------- sequential_update: I = ",I," J = ", J," w = ",w)

            pos = randomstate.choice(len(II), 1, p=w)
            # print("--------- sequential_update: pos = ",pos)
            people.data[i, j] = 0
            people.data[II[pos[0]], JJ[pos[0]]] = 1
            people_ij[id, 0] = II[pos[0]]
            people_ij[id, 1] = JJ[pos[0]]
    # print("===> Move: people_ij = ",people_ij)
    return people, people_ij


def parallel_update(people, people_ij, weight, friction=0, randomstate=None):
    """
    To move all individuals in parallel according to the following rule: \
    first, desired moves are precomputed and then the conflicts (two individuals \
    at the same position) are resolved

    Parameters
    ----------

    people: numpy masked arrays
        equal to 1 if the cell (i,j) is occupied, 0 elsewhere
    people_ij: numpy array
        (i,j) for each individual
    weight: numpy array
        weights for the probabilities in order to move in such a way like to \
        reach the door
    friction: float
        to designate the effect induced by a modified handling of conflicts, \
        friction is the probability that a conflict remains unresolved (no one \
        moves)
    randomstate: numpy randomstate
        create a new one or reuse the given random state

    Returns
    -------

    people: numpy masked arrays
        new positions equal to 1 if the cell (i,j) is occupied, 0 elsewhere
    people_id: numpy masked arrays
        new people index (i,j)
    """
    height, width = people.shape
    nn = people_ij.shape[0]
    if (randomstate is None):
        randomstate = np.random.RandomState()
    new_ij = -np.ones((nn, 3))  # new_i new_j proba
    for id in np.arange(nn):
        i = people_ij[id, 0]
        j = people_ij[id, 1]
        if ((i != -1) and (j != -1)):
            II = [i]
            JJ = [j]
            if (i >= 1):
                im = i - 1
                if ((people[im, j] is not np.ma.masked) and (people[im, j] != 1)):
                    II.append(im)
                    JJ.append(j)
            if (i <= height - 2):
                ip = i + 1
                if ((people[ip, j] is not np.ma.masked) and (people[ip, j] != 1)):
                    II.append(ip)
                    JJ.append(j)
            if (j >= 1):
                jm = j - 1
                if ((people[i, jm] is not np.ma.masked) and (people[i, jm] != 1)):
                    II.append(i)
                    JJ.append(jm)
            if (j <= width - 2):
                jp = j + 1
                if ((people[i, jp] is not np.ma.masked) and (people[i, jp] != 1)):
                    II.append(i)
                    JJ.append(jp)
            w = weight[II, JJ] / np.sum(weight[II, JJ])
            pos = randomstate.choice(len(II), 1, p=w)
            new_ij[id, 0] = II[pos[0]]
            new_ij[id, 1] = JJ[pos[0]]
            new_ij[id, 2] = w[pos[0]]
    ind = np.where(new_ij[:, 0] != -1)[0]
    # Find the conflicts
    unique_ij, index, inverse, counts = np.unique(new_ij[ind, :2], axis=0, return_index=True, return_inverse=True,
                                                  return_counts=True)  # numpy 1.13 only !!
    if (unique_ij.shape[0] < ind.shape[0]):
        for ic, cc in enumerate(counts):
            if (cc > 1):  # Conflict...
                ind_conf = np.where(inverse == ic)[0]  # Positions of people concerned by the Conflicts
                npc = ind_conf.shape[0]  # Number of concerned persons
                print("---------- Conflicts between ", npc, " persons")
                # No agent will move to the empty cell with the probability "friction"
                choice = randomstate.choice(2, 1, p=[friction, 1 - friction])[0]  # 0 = no agent move...
                if (choice == 0):
                    print("---------- => Friction... No agent move...")
                    for id in np.arange(ind_conf.shape[0]):
                        new_ij[ind[ind_conf[id]], :2] = people_ij[ind[ind_conf[id]], :2]
                else:
                    print("---------- => One agent move...")
                    w = new_ij[ind[ind_conf], 2] / np.sum(new_ij[ind[ind_conf], 2])
                    pos = randomstate.choice(ind_conf.shape[0], 1, p=w)  # Winner
                    for id in np.arange(ind_conf.shape[0]):
                        if (id != pos):
                            new_ij[ind[ind_conf[id]], :2] = people_ij[ind[ind_conf[id]], :2]
    # Update
    for id in np.arange(nn):
        i = people_ij[id, 0]
        j = people_ij[id, 1]
        II = int(new_ij[id, 0])
        JJ = int(new_ij[id, 1])
        people.data[i, j] = 0
        people.data[II, JJ] = 1
        people_ij[id, 0] = II
        people_ij[id, 1] = JJ
    # print("===> Move: people_ij = ",people_ij)
    return people, people_ij


def exit(domain, people, people_ij):
    """
    To update people and people_ij arrays in removing individuals who \
    left the domain

    Parameters
    ----------

    domain: Domain
        contains everything for managing the domain
    people: numpy masked arrays
        equal to 1 if the cell (i,j) is occupied, 0 elsewhere
    people_ij: numpy array
        (i,j) for each individual

    Returns
    -------

    people: numpy masked arrays
        new positions equal to 1 if the cell (i,j) is occupied, 0 elsewhere
    people_id: numpy masked arrays
        new people index (i,j)
    """
    Sexit = 0.2
    ind_in = np.where(people_ij[:, 0] != -1)[0]
    d = domain.destinations["door"].distance[people_ij[ind_in, 0], people_ij[ind_in, 1]]
    dind = np.where(d < Sexit)[0]
    if (dind.shape[0] > 0):
        people.data[people_ij[ind_in[dind], 0], people_ij[ind_in[dind], 1]] = 0
        people_ij[ind_in[dind], :] = [-1, -1]
        # print("===> Exit: people_ij = ",people_ij)
    return people, people_ij
