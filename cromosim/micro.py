# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
# License: GPL

import numpy as np
import scipy as sp
import sys
import random
import matplotlib
import matplotlib.pyplot as plt

#from scipy.spatial import KDTree
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from matplotlib.patches import Ellipse, Circle, Rectangle, Polygon, Arrow
from matplotlib.lines import Line2D
from matplotlib.collections import EllipseCollection, LineCollection


def compute_contacts(dom, xyrv, dmax):
    """This function uses a KDTree method to find the contacts between
    individuals. Moreover the contacts with the walls are also determined
    from the wall distance (obtained by the fast-marching method).

    Parameters
    ----------
    dom: Domain
        contains everything for managing the domain
    xyrv: numpy array
        people coordinates, radius and velocity coefficient: ``x,y,r,v``
    dmax: float
        threshold value used to consider a contact as active i.e. ``dij<dmax``

    Returns
    -------
    contacts: numpy array
        all the contacts ``i,j,dij,eij_x,eij_y`` such that ``dij<dmax`` and
        ``i<j`` (no duplication)
    """
    # lf: the number of points at which the algorithm
    # switches over to brute-force. Has to be positive.
    lf = 100
    if (lf>sys.getrecursionlimit()):
        sys.setrecursionlimit(lf)
    kd = cKDTree(xyrv[:,:2],leafsize = lf)
    ## Find all pairs of points whose distance is at most dmax+2*rmax
    rmax = xyrv[:,2].max()
    neighbors = kd.query_ball_tree(kd,dmax+2*rmax)
    ## Create the contact array: i,j,dij,eij_x,eij_y
    first_elements = np.arange(xyrv.shape[0]) ## i
    other_elements = list(map(lambda x: x[1:], neighbors)) ## all j for each i
    lengths = list(map(len, other_elements))
    tt = np.stack([first_elements,lengths],axis=1)
    I = np.concatenate(list(map(lambda x: np.full((x[1],), x[0]), tt))).astype(int)
    J = np.concatenate(other_elements).astype(int)
    ind = np.where(I<J)[0]
    I = I[ind] ; J = J[ind]
    DP = xyrv[J,:2]-xyrv[I,:2]
    Norm = np.linalg.norm(DP,axis=1,ord=2)
    Dij = Norm - xyrv[I,2]-xyrv[J,2]
    ind = np.where(Dij<dmax)[0]
    Dij = Dij[ind]
    I = I[ind]
    J = J[ind]
    Norm = Norm[ind]
    DP = DP[ind]
    contacts = np.stack([I,J,Dij,DP[:,0]/Norm,DP[:,1]/Norm],axis=1)
    # Add contacts with the walls
    II,JJ,DD = dom.people_wall_distance(xyrv)
    ind = sp.where(DD<dmax)[0]
    wall_contacts = sp.stack([ind,-1*sp.ones(ind.shape),DD[ind],
                             dom.wall_grad_X[II[ind],JJ[ind]],
                             dom.wall_grad_Y[II[ind],JJ[ind]] ],axis=1)
    contacts = np.vstack([contacts,wall_contacts])
    return np.array(contacts)


def compute_forces(F, Fwall, xyrv, contacts, U, Vd, lambda_, delta, k, eta):
    """This function computes all the forces (isentropic interaction and
    friction) and sums them. The correcting pre-factor due to the vision
    angle is also used into the social force term.

    Parameters
    ----------
    F: float
        social trend of an individual to keep apart from another (homogeneous
        to a force)
    Fwall: float
        social trend of an individual to keep apart from a wall (homogeneous
        to a force)
    xyrv: numpy array
        people coordinates, radius and velocity coefficient: ``x,y,r,v``
    contacts: numpy array
        all the contacts: ``i,j,dij,eij_x,eij_y``
    U: numpy array
        people velocities
    Vd: numpy array
        people desired velocities
    lambda_: float
        quantifies the directional dependence when the vision angle is
        considered (between ``[0,1]``, if equal to 1 the fully isotropic case is
        recovered)
    delta: float
        maintains a certain distance between neighbors
    k: float
        used when there is overlapping, k is a stiffness constant of
        individuals seen as deformable bodies
    eta: float
        friction coefficient

    Returns
    -------
    Forces: numpy array
        sum of all forces for each individual
    """
    Np = xyrv.shape[0]
    Nc = contacts.shape[0]
    Forces = np.zeros((Np,2))
    ## Social forces, friction,...
    for ic in np.arange(Nc):
        i = int(contacts[ic,0])
        j = int(contacts[ic,1])
        dij = contacts[ic,2]
        dij_moins = min(dij,0.0)
        eij_x = contacts[ic,3]
        eij_y = contacts[ic,4]
        if (j>-1): ## contact person/person
            # Angular dependence
            norm_Vdi = np.sqrt( Vd[i,0]**2+Vd[i,1]**2 )
            if (norm_Vdi > 0):
                theta_ij = np.arccos(  (Vd[i,0]*eij_x+Vd[i,1]*eij_y)/norm_Vdi )
                omega_ij = lambda_+(1-lambda_)*(1+np.cos(theta_ij))/2
            else:
                omega_ij= 1
            norm_Vdj = np.sqrt( Vd[j,0]**2+Vd[j,1]**2 )
            if (norm_Vdj > 0):
                theta_ji = np.arccos( -(Vd[j,0]*eij_x+Vd[j,1]*eij_y)/norm_Vdj )
                omega_ji = lambda_+(1-lambda_)*(1+np.cos(theta_ji))/2
            else:
                omega_ji = 1
            # Social force + force to handle overlapping
            fij = -omega_ij*F*np.exp(-dij/delta) + k*dij_moins
            fji = -omega_ji*F*np.exp(-dij/delta) + k*dij_moins
            Forces[i,0] += fij*eij_x
            Forces[i,1] += fij*eij_y
            Forces[j,0] -= fji*eij_x
            Forces[j,1] -= fji*eij_y
            # Friction
            fij_friction = eta*dij_moins*( -(U[i,0]-U[j,0])*eij_y
                + (U[i,1]-U[j,1])*eij_x )
            Forces[i,0] -= fij_friction*eij_y
            Forces[i,1] += fij_friction*eij_x
            Forces[j,0] += fij_friction*eij_y
            Forces[j,1] -= fij_friction*eij_x
        else: ## contact person/walls
            fij = -Fwall*np.exp(-dij/delta) + k*dij_moins
            Forces[i,0] -= fij*eij_x
            Forces[i,1] -= fij*eij_y
    return Forces


def projection(dt, xyrv, contacts, Vd, dmin_people=0, dmin_walls=0,
        nb_iter_max = 100000, rho=0.1, tol = 0.01, log=False, method="cvxopt"):
    """From the desired velocities Vd, this projection step consists of
    computing the global velocity field defined as the closest velocity to the
    desired one among all the feasible fields (i.e. fields which do not lead 
    to disks overlapping).

    Parameters
    ----------
    dt: float
        time step
    xyrv: numpy array
        people coordinates, radius and velocity coefficient: ``x,y,r,v``
    contacts: numpy array
        all the contacts: ``i,j,dij,eij_x,eij_y``
    Vd: numpy array
        people desired velocities
    dmin_people: float
        minimum distance guaranteed between individuals
    dmin_walls: float
        minimum distance guaranteed between an individual and a wall
    nb_iter_max: integer
        maximum number of iterations allowed
    rho: float
        parameter of the Uzawa method
    tol: float
        tolerance wished
    log: boolean
        to print the final accuracy, number of iterations,...
    method: string
        optimization algorithm: ``"cvxopt"`` (default) or ``"uzawa"`` (or
        ``"mosek"`` if installed with a license file).

    Returns
    -------
    info: integer
        number of iterations needed
    B: numpy array
        constraint matrix
    U: numpy array
        new people velocities ensuring that there is no overlap between
        individuals
    L: numpy array
        Lagrange multipliers (only when ``method="uzawa"``, None otherwise)
    P: numpy array
        pressure on each individual (only when ``method="uzawa"``, None
        otherwise)
    """
    Np = xyrv.shape[0]
    Nc = contacts.shape[0]
    info = 0
    if (Nc == 0):
        info = 1
        return info, None, Vd, None, None
    else:

        if (method == "cvxopt") or (method == "mosek"):

            import cvxopt
            cvxopt.solvers.options['show_progress'] = False
            cvxopt.solvers.maxiters = 1000
            cvxopt.solvers.abstol = 1e-8
            cvxopt.solvers.reltol = 1e-7
            L = None
            P = None
            U = np.zeros((2*Np,))
            V = np.zeros((2*Np,))
            ind_contacts_walls = np.where(contacts[:,1] ==-1)[0]
            DMIN = np.ones(contacts.shape[0])*dmin_people
            DMIN[ind_contacts_walls] = dmin_walls
            Z = (contacts[:,2]-DMIN)/dt ## ie Dij/dt
            V[::2] = Vd[:,0]; V[1::2] = Vd[:,1] ## A priori velocity
            V = cvxopt.matrix(V)
            Z = cvxopt.matrix(Z, (Nc,1))
            Id = cvxopt.spdiag([1]*(U.shape[0]))
            if (Nc>0):
                II = contacts[:,0].astype(int)
                JJ = contacts[:,1].astype(int)
                Jpos = np.where(JJ>=0)[0]
                Jneg = np.where(JJ<0)[0]
                row = np.concatenate([Jpos, Jpos, Jpos, Jpos, Jneg, Jneg])
                col = np.concatenate([2*II[Jpos], 2*II[Jpos]+1,
                                      2*JJ[Jpos], 2*JJ[Jpos]+1,
                                      2*II[Jneg], 2*II[Jneg]+1])
                data = np.concatenate([contacts[Jpos,3], contacts[Jpos,4],
                                       -contacts[Jpos,3], -contacts[Jpos,4],
                                       -contacts[Jneg,3], -contacts[Jneg,4]])
                B = csr_matrix((data, (row, col)), shape=(Nc, 2*Np))#.toarray()
                cvxoptB = cvxopt.spmatrix(np.array(data),np.array(row),
                                          np.array(col),size=(Nc, 2*Np))
                if (method == "mosek"):
                    from mosek import iparam
                    cvxopt.solvers.options['mosek'] = {iparam.log: 0}
                    solution = cvxopt.solvers.qp(Id, -V, cvxoptB, Z, solver='mosek')
                else:
                    solution = cvxopt.solvers.qp(Id, -V, cvxoptB, Z)
                    info = solution["iterations"]
                if log:
                    C = Z - B@U
                    if (method == "mosek"):
                        print("    projection (mosek): nb of contacts = ",Nc,
                              ", contrainte (Z-B@U).min() = ",C.min())
                    else:
                        print("    projection (cvxopt): nb of contacts = ",Nc,
                              ", nb of iterations = ",solution["iterations"],
                              ", status = ",solution["status"],
                              ", contrainte (Z-B@U).min() = ",C.min())
                if (solution["status"] == "optimal"):
                    U = solution['x']
                else:
                    print("    projection (mosek or cvxopt): optimization failed")
                    print("    ---> try with uzawa method...")
                    info = 0
                    II = contacts[:,0].astype(int)
                    JJ = contacts[:,1].astype(int)
                    Jpos = np.where(JJ>=0)[0]
                    Jneg = np.where(JJ<0)[0]
                    row = np.concatenate([Jpos, Jpos, Jpos, Jpos, Jneg, Jneg])
                    col = np.concatenate([2*II[Jpos], 2*II[Jpos]+1,
                                          2*JJ[Jpos], 2*JJ[Jpos]+1,
                                          2*II[Jneg], 2*II[Jneg]+1])
                    data = np.concatenate([contacts[Jpos,3], contacts[Jpos,4],
                                           -contacts[Jpos,3], -contacts[Jpos,4],
                                           -contacts[Jneg,3], -contacts[Jneg,4]])
                    B = csr_matrix((data, (row, col)), shape=(Nc, 2*Np))
                    L = np.zeros((Nc,))
                    R = 99*np.ones((Nc,))
                    U = np.zeros((2*Np,))
                    V = np.zeros((2*Np,))
                    D = contacts[:,2]
                    V[::2] = Vd[:,0]; V[1::2] = Vd[:,1]
                    k = 0
                    ind_contacts_walls = np.where(contacts[:,1] ==-1)[0]
                    DMIN = np.ones(contacts.shape[0])*dmin_people
                    DMIN[ind_contacts_walls] = dmin_walls
                    while (( dt*R.max()>tol*2*xyrv[:,2].min()) and (k<nb_iter_max)):
                        U[:] = V[:] - B.transpose()@L[:]
                        R[:] = B@U[:] - (D[:]-DMIN)/dt
                        L[:] = np.maximum(L[:] + rho*R[:], 0)
                        k += 1
                    if log:
                        print("    projection (uzawa): nb of contacts = ",Nc,
                              ", nb of iterations = ",k,", min = ",R.min(),
                              ", max = ",R.max(),", tol = ",tol)
                    if (k==nb_iter_max):
                        print("** WARNING: Method projection **")
                        print("** WARNING: you have reached the maximum number \
                               of iterations,")
                        print("** WARNING: it remains unsatisfied constraints \
                               !! ")
                        info = -1
                    else:
                        info = k

            U = np.array(U).reshape((Np, 2))

        elif (method=="uzawa"):

            info = 0
            II = contacts[:,0].astype(int)
            JJ = contacts[:,1].astype(int)
            Jpos = np.where(JJ>=0)[0]
            Jneg = np.where(JJ<0)[0]
            row = np.concatenate([Jpos, Jpos, Jpos, Jpos, Jneg, Jneg])
            col = np.concatenate([2*II[Jpos], 2*II[Jpos]+1,
                                  2*JJ[Jpos], 2*JJ[Jpos]+1,
                                  2*II[Jneg], 2*II[Jneg]+1])
            data = np.concatenate([contacts[Jpos,3], contacts[Jpos,4],
                                   -contacts[Jpos,3], -contacts[Jpos,4],
                                   -contacts[Jneg,3], -contacts[Jneg,4]])
            B = csr_matrix((data, (row, col)), shape=(Nc, 2*Np))#.toarray()
            L = np.zeros((Nc,))
            R = 99*np.ones((Nc,))
            U = np.zeros((2*Np,))
            V = np.zeros((2*Np,))
            D = contacts[:,2]
            V[::2] = Vd[:,0]; V[1::2] = Vd[:,1]
            k = 0
            ind_contacts_walls = np.where(contacts[:,1] ==-1)[0]
            DMIN = np.ones(contacts.shape[0])*dmin_people
            DMIN[ind_contacts_walls] = dmin_walls
            while (( dt*R.max()>tol*2*xyrv[:,2].min()) and (k<nb_iter_max)):
                U[:] = V[:] - B.transpose()@L[:]
                R[:] = B@U[:] - (D[:]-DMIN)/dt
                L[:] = np.maximum(L[:] + rho*R[:], 0)
                k += 1
            P = np.zeros(Np) ## Pressure
            P[II[Jpos]] +=  3/(4*np.pi*xyrv[II[Jpos],2]**2)*L[Jpos]
            P[JJ[Jpos]] +=  3/(4*np.pi*xyrv[JJ[Jpos],2]**2)*L[Jpos]
            P[II[Jneg]] +=  3/(4*np.pi*xyrv[II[Jneg],2]**2)*L[Jneg]
            if log:
                print("    projection (uzawa): nb of contacts = ",Nc,
                      ", nb of iterations = ",k,", min = ",R.min(),
                      ", max = ",R.max(),", tol = ",tol)
            if (k==nb_iter_max):
                print("** WARNING: Method projection **")
                print("** WARNING: you have reached the maximum number \
                       of iterations,")
                print("** WARNING: it remains unsatisfied constraints !! ")
                info = -1
            else:
                info = k

        return info, B, U.reshape((Np, 2)), L, P


def move_people(time, dt, people, sensors):
    """
    Updates the people positions according to the new velocities U. \
    If there exists crosslines (i.e. sensors), computes also the id, time, \
    direction and impact points for the individuals who cross the lines.

    Parameters
    ----------
    time: float
        current time
    dt: float
        time step
    people: dict
        dictionary containing everything related to individuals (``"xyrv"``, \
        ``"U"``,...)
    sensors: dict
        dictionary containing everything related to sensors (``"line"``,...)

    Returns
    -------
    people: dict
        dictionary updated with new people positions after moving
    sensors: dict
        dictionary enriched by:
         * ``"id"``: people id who cross the sensor lines
         * ``"times"``: times when people cross the sensor lines
         * ``"xy"``: impact points for people crossing the sensor lines
         * ``"dir"``: directions for people crossing the sensor lines (entry or exit)
    """
    if len(sensors)>0:
        xyrv_old = people["xyrv"].copy()
    ## Update people coordinates
    people["xyrv"][:,0] += dt*people["U"][:,0]
    people["xyrv"][:,1] += dt*people["U"][:,1]
    for il,s in enumerate(sensors):
        l = s["line"]
        # line [x0,y0,x1,y1],
        id, pt, dir, times, entries, exits = sensor(l, xyrv_old[:,:2],
                        people["xyrv"][:,:2], time, time+dt)
        sensors[il]["id"] = sensors[il]["id"] + people["id"][id].tolist()
        sensors[il]["times"] = sensors[il]["times"] + times.tolist()
        sensors[il]["xy"] = sensors[il]["xy"] + pt.tolist()
        sensors[il]["dir"] = sensors[il]["dir"] + dir.tolist()
    return people, sensors


##----------------------------------------------------------------##
## Deprecated function (used only in the version 1.0 of cromosim) ##
##----------------------------------------------------------------##
# def exit_door(sexit, dom, people, people_dest, U, arrays=[]):
#     """
#     Removes individuals who are at a distance less than sexit \
#     to the closest door
#
#     Parameters
#     ----------
#     sexit: float
#         threshold value used to remove individuals close to the exit
#     dom: Domain
#         contains everything for managing the domain
#     people: numpy array
#         people coordinates and radius: x,y,r
#     U: numpy array
#         people velocities
#     arrays: list of numpy array
#         other arrays to resize similarly as people and U arrays
#
#     Returns
#     -------
#     people: numpy array
#         new people positions (outside individuals had been removed)
#     U: numpy array
#         new poeple velocities (outside individuals had been removed)
#     arrays: list of numpy array
#         new array resized
#     """
#     I,J,D = dom.people_target_distance(people,people_dest)
#     ind = sp.where(D>sexit)
#     if (len(arrays)>0):
#         return people[ind[0],:], U[ind[0],:], [ a[ind[0]] for a in arrays]
#     else:
#         return people[ind[0],:], U[ind[0],:]


##----------------------------------------------------------------##
## Deprecated function (used only in the version 1.0 of cromosim) ##
##----------------------------------------------------------------##
# def exit_out_of_domain(dom, people, arrays=[], box=None):
#     """
#     Removes individuals who are outside the domain or outside a given box
#
#     Parameters
#     ----------
#     dom: Domain
#         contains everything for managing the domain
#     people: numpy array
#         people coordinates and radius: x,y,r
#     arrays: list of numpy array
#         other arrays to resize similarly as people and U
#     box: numpy array
#         box coordinates [xmin,xmax,ymin,ymax] which replace the \
#         domain minimum and maximum coordinates
#
#     Returns
#     -------
#     people: numpy array
#         new people array (outside individuals had been removed)
#     arrays: list of numpy array
#         new arrays resized similarly as people array
#     """
#     if box is None:
#         ## Remove people who are outside the domain
#         S = (people[:,0]-people[:,2]<=dom.xmin+dom.pixel_size) + \
#             (people[:,0]-people[:,2]>=dom.xmax-dom.pixel_size) + \
#             (people[:,1]-people[:,2]<=dom.ymin+dom.pixel_size) + \
#             (people[:,1]-people[:,2]>=dom.ymax-dom.pixel_size)
#     else:
#         ## Remove people who are outside the given box
#         S = (people[:,0]-people[:,2]<=box[0]+dom.pixel_size) + \
#             (people[:,0]-people[:,2]>=box[1]-dom.pixel_size) + \
#             (people[:,1]-people[:,2]<=box[2]+dom.pixel_size) + \
#             (people[:,1]-people[:,2]>=box[3]-dom.pixel_size)
#     ind = np.where(S==False)[0]
#     people = people[ind,:]
#     if (len(arrays)>0):
#         for a in arrays:
#             a = a[ind]
#     ## Remove people who are too close to walls or with a masked door distance
#     I,J,Dwall = dom.people_wall_distance(people)
#     I,J,Ddoor = dom.people_door_distance(people, I=I, J=J)
#     indDwall = sp.where(Dwall<=dom.pixel_size)[0]
#     indDdoor = sp.where(Ddoor.mask==True)[0]
#     ind = sp.unique(sp.concatenate((indDwall,indDdoor)))
#     comp_ind = sp.setdiff1d(sp.arange(people.shape[0]), ind)
#     if (len(arrays)>0):
#         return people[comp_ind,:], [ a[comp_ind] for a in arrays]
#     else:
#         return people[comp_ind,:]


##----------------------------------------------------------------##
## Deprecated function (used only in the version 1.0 of cromosim) ##
##----------------------------------------------------------------##
# def exit_box(box, sexit, people, U, arrays=[]):
#     """
#     Removes individuals outside the box according to a threshold value.
#
#     Parameters
#     ----------
#     box: numpy array
#         box vertice coordinates [xmin, xmax, ymin, ymax]
#     sexit: float
#         threshold value used to remove individuals close to the exit
#     people: numpy array
#         people coordinates and radius: x,y,r
#     U: numpy array
#         people velocities
#     arrays: list of numpy array
#         other arrays to resize similarly as people and U
#
#     Returns
#     -------
#     people: numpy array
#         new people coordinates (outside individuals had been removed)
#     U: numpy array
#         new people velocities (outside individuals had been removed)
#     arrays: list of numpy array
#         new resized arrays taking into account deleted individuals
#     """
#     S = (people[:,0]-people[:,2]<=box[0]+sexit) + \
#         (people[:,0]-people[:,2]>=box[1]-sexit) + \
#         (people[:,1]-people[:,2]<=box[2]+sexit) + \
#         (people[:,1]-people[:,2]>=box[3]-sexit)
#     ind = np.where(S==False)
#     if (len(arrays)>0):
#         return people[ind[0],:], U[ind[0],:], [ a[ind[0]] for a in arrays]
#     else:
#         return people[ind[0],:], U[ind[0],:]


##----------------------------------------------------------------##
## Deprecated function (used only in the version 1.0 of cromosim) ##
## BUT this function will be adapted to version 2.0 of cromosim   ##
##----------------------------------------------------------------##
# def periodic_bc_vertical(ymin, ymax, people, U, xmin=None, xmax=None, rng=None):
#     """
#     Does not exactly correspond to periodic boundary conditions (in the
#     mathematical sense): the persons having an y-coordinate y less than ymin
#     (respectively greater than ymax) are reinjected at y+(ymax-ymin) (respectively
#     at y-(ymax-ymin)) with (optionally) random x-coordinates (between xmin and xmax,
#     and with velocities equal to 0.
#
#     Parameters
#     ----------
#     ymin: float
#         minimal ordinate of the box
#     ymax: float
#         maximal ordinate of the box
#     people: numpy array
#         people coordinates and radius: x,y,r
#     U: numpy array
#         people velocities
#     xmin: float
#         minimal abscissa of the box
#     xmax: float
#         maximal abscissa of the box
#     rng: RandomState
#         scipy random state object (see scipy.random.RandomState)
#
#     Returns
#     -------
#     people: numpy array
#         new people coordinates (outside individuals had been moved)
#     U: numpy array
#         new people velocities
#     """
#     Smin = (people[:,1]<ymin)
#     Smax = (people[:,1]>ymax)
#     indmin = np.where(Smin==True)[0]
#     indmax = np.where(Smax==True)[0]
#     people[indmin,1] += (ymax-ymin)
#     people[indmax,1] -= (ymax-ymin)
#     U[indmin,:]=0
#     U[indmax,:]=0
#     if xmin is not None:
#         if xmax is not None:
#             if rng is None:
#                 rng = np.random.RandomState()
#             x_indmin = rng.uniform(xmin, xmax, indmin.shape[0])
#             people[indmin,0] = x_indmin
#             x_indmax = rng.uniform(xmin, xmax, indmax.shape[0])
#             people[indmax,0] = x_indmax
#     return people, U


def create_people_in_box(nn, box, dest_name, radius_distribution,
        velocity_distribution, dom, rng, verbose=True):
    """
    To create nn persons in a given box. The overlaps are not treated for \
    the moment but one checks that the individuals are located in an area where \
    the desired velocity is well defined (outside inaccessible areas or \
    obstacles). If it is not the case, one changes their coordinates in \
    consequence.

    Parameters
    ----------
    nn: integer
        number of individuals to create
    box: list
        coordinates of the box ``[xmin, xmax, ymin, ymax]``
    dest_name: string
        destination name for all individuals in this group
    radius_distribution: list
        distribution with its parameters
    velocity_distribution: list
        distribution with its parameters
    dom: Domain
        contains everything for managing the domain
    rng: RandomState
        scipy random state object (see ``scipy.random.RandomState``)
    verbose: boolean
        logs for debug

    Returns
    -------
    p: numpy array
        new people coordinates, radius and velocity coefficient ``x y r v``
    dest_names: numpy array
        people destination names
    """
    if (verbose):
        print("------ create_people_in_box --> Create "+str(nn)+ \
          " individuals in the box "+str(box)+", overlaps can occur...")
    xyrv = np.zeros((nn,4))  # x y r v
    dest_names = [dest_name]*nn
    if (radius_distribution[0]=="uniform"):
        rmin = radius_distribution[1]
        rmax = radius_distribution[2]
        xyrv[:,2] = rng.uniform(rmin, rmax, nn)
    elif (radius_distribution[0]=="normal"):
        mean = radius_distribution[1]
        sigma = radius_distribution[2] # standard deviation
        xyrv[:,2] = rng.normal(mean, sigma, nn)  # r: radius
    else:
        print("------ create_people_in_box --> unknown distribution for \
               people radius ")
    p_rmax = xyrv[:,2].max()
    if (velocity_distribution[0]=="uniform"):
        vmin = velocity_distribution[1]
        vmax = velocity_distribution[2]
        xyrv[:,3] = rng.uniform(vmin, vmax, nn) # v: velocity coefficient
    elif (velocity_distribution[0]=="normal"):
        mean = velocity_distribution[1]
        sigma = velocity_distribution[2] # standard deviation
        xyrv[:,3] = rng.normal(mean, sigma, nn) # v: velocity coefficient
    else:
        print("------ create_people_in_box --> unknown distribution for \
               people velocity")
    xmin, xmax, ymin, ymax = box
    xyrv[:,0] = rng.uniform(xmin+p_rmax, xmax-p_rmax, nn)
    xyrv[:,1] = rng.uniform(ymin+p_rmax, ymax-p_rmax, nn)
    ## To check if people are well localized in the domain
    ## i.e. where a desired velocity is defined
    while True:
        I, J, Vd = dom.people_desired_velocity(xyrv,dest_names)
        normVd = Vd[:,0]**2+Vd[:,1]**2
        ind = np.where(normVd==0)[0]
        if (ind.shape[0]>0):
            xyrv[ind,0] = rng.uniform(xmin+p_rmax, xmax-p_rmax, ind.shape[0])
            xyrv[ind,1] = rng.uniform(ymin+p_rmax, ymax-p_rmax, ind.shape[0])
        else:
            break
    return xyrv, dest_names


def check_people_in_box(dom, box, xyrv, dest_names, rng, verbose=True):
    """To check that people coordinates are in the given box (test 1) and in an
    usable space i.e. in an area accessible and not concerned by obstacles
    (test 2). On the other hand, one moves the individuals which do not satisfy
    these two tests.

    Parameters
    ----------
    dom: Domain
        contains everything for managing the domain
    box: list
        coordinates of the box [xmin, xmax, ymin, ymax]
    xyrv: numpy array
        people coordinates, radius and velocity coefficient: ``x,y,r,v``
    dest_names: numpy array
        people destination names
    rng: RandomState
        scipy random state object (see scipy.random.RandomState)
    verbose: boolean
        logs for debug

    Returns
    -------
    xyrv: numpy array
        new people coordinates ``x y r``

    """
    if (verbose):
        print("------ check_people_in_box --> To verify that "+str(xyrv.shape[0])+ \
          " individuals are in the domain, in the box and with a defined"+ \
          " desired velocity")
    p_rmax = xyrv[:,2].max()
    xmin, xmax, ymin, ymax = box
    info = False
    while True:
        ## test 1
        I = np.floor((xyrv[:,1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
        J = np.floor((xyrv[:,0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
        test1 = (I>=0)*(I<dom.height)*(J>=0)*(J<dom.width)
        ind1 = np.where(test1==0)[0]
        if (ind1.shape[0]>0):
            if (verbose):
                print("------ check_people_in_box --> "+str(ind1.shape[0])+ \
                  " individuals outside the domain")
            info = True
            xyrv[ind1,0] = rng.uniform(xmin+p_rmax, xmax-p_rmax, ind1.shape[0])
            xyrv[ind1,1] = rng.uniform(ymin+p_rmax, ymax-p_rmax, ind1.shape[0])
        else:
            ## test 2
            I, J, Vd = dom.people_desired_velocity(xyrv,dest_names)
            normVd = Vd[:,0]**2+Vd[:,1]**2
            test2 = (xyrv[:,0]>xmin+p_rmax)*(xyrv[:,0]<xmax-p_rmax) \
                  *(xyrv[:,1]>ymin+p_rmax)*(xyrv[:,1]<ymax-p_rmax) \
                  *(normVd>0)
            ind2 = np.where(test2==0)[0]
            if (ind2.shape[0]>0):
                if (verbose):
                    print("------ check_people_in_box --> "+str(ind2.shape[0])+ \
                      " individuals with an undefined desired velocity ")
                info = True
                xyrv[ind2,0] = rng.uniform(xmin+p_rmax, xmax-p_rmax, ind2.shape[0])
                xyrv[ind2,1] = rng.uniform(ymin+p_rmax, ymax-p_rmax, ind2.shape[0])
            else:
                if (verbose):
                    print("------ check_people_in_box --> OK !")
                break
    return info, xyrv


def remove_overlaps_in_box(dom, box, xyrv, dest_names, dt, rng,
        projection_method="cvxopt", dmin_people=0, dmin_walls=0, itermax=10,
        verbose=True):
    """
    To remove the overlaps between individuals (spheres) in the given box. \
    Several projections are used to give a better robustness to this \
    process when the number of overlaps is very high e.g. during the \
    initialization when the individuals have random positions.

    Parameters
    ----------
    dom: Domain
        contains everything for managing the domain
    box: list
        coordinates of the box ``[xmin, xmax, ymin, ymax]``
    xyrv: numpy array
        people coordinates ``x y``, radius ``r`` and velocity coefficient ``v``
    dest_names: numpy array
        people destination names
    dt: float
        time step
    rng: RandomState
        scipy random state object (see ``scipy.random.RandomState``)
    dmin_people: float
        minimal distance allowed between individuals
    dmin_walls: float
        minimal distance allowed between an individual and a wall
    itermax: integer
        maximal number of Uzawa projections (10 by default)
    verbose: boolean
        logs for debug

    Returns
    -------
    xyrv: numpy array
        new people coordinates ``x y``, same radius ``r`` and velocity
        coefficient ``v``
    """
    if (verbose):
        print("------ remove_overlaps_in_box --> Remove overlaps...")
    dmax_init = xyrv[:,2].max() # maximal radius
    it1 = 0
    t = 0
    while (True):
        if (verbose):
            print("------ remove_overlaps_in_box --> Remove overlaps in box "
              + str(box)+": iteration "+str(it1)+" / "+str(itermax))
        c = compute_contacts(dom, xyrv, dmax_init)
        I, J, Vd = dom.people_desired_velocity(xyrv,dest_names)
        info1, B, U, L, P = projection(dt, xyrv, c, Vd,
                                       method=projection_method,
                                       dmin_people=dmin_people,
                                       dmin_walls=dmin_walls,
                                       nb_iter_max = 10000, log = False)
        # Move people
        xyrv[:,0] += dt*U[:,0]
        xyrv[:,1] += dt*U[:,1]

        it1 += 1
        info2, people = check_people_in_box(dom, box, xyrv, dest_names,
                                            rng, verbose=verbose)
        if ((info1>=0) and (info2==False)):
            if (verbose):
                print("------ remove_overlaps_in_box --> Successful... \
                       No overlaps")
            break
        if (it1>=itermax):
            if (verbose):
                print("------ remove_overlaps_in_box --> ** WARNING: \
                       Unsuccessful initialization **")
                print("------ remove_overlaps_in_box --> ** WARNING: Let us \
                       continue...")
            break
    return xyrv


def people_initialization(dom, groups, dt, dmin_people=0, dmin_walls=0,
        seed=0, itermax=10, projection_method="cvxopt", verbose=True):
    """To initialize people.

    Parameters
    ----------
    dom: Domain
        contains everything for managing the domain
    groups: dict
        contains all groups related to this domain. A group contains:
         * ``"nb"``: number of persons
         * ``"radius_distribution"``: ``["uniform",min,max]`` or ``["normal",mean,std_dev]``
         * ``"velocity_distribution"``: ``["uniform",min,max]`` or ``["normal",mean,std_dev]``
         * ``"box"``: ``[xmin,xmax,ymin,ymax]``
         * ``"destination"``: initial destination
    dt: float
        time step
    dmin_people: float
        minimal distance allowed between individuals (0 by default)
    dmin_walls: float
        minimal distance allowed between an individual and a wall
    seed: integer
        seed to initialize the RandomState (0 by default means different seed \
        at each run)
    itermax: integer
        maximal number of Uzawa projections (10 by default)
    projection_method: string
        optimizer name for the projection step: ``"cvxopt"``, ``"uzawa"``
        or ``"mosek"``
    verbose: boolean
        logs for debug

    Returns
    -------
    people: dict
        people dictionary which contains:
         * ``"xyrv"``: coordinates, radius and velocity coefficient for each individual
         * ``"destinations"``: destination for each individual
         * ``"gpid"``: group id

    rng: RandomState
        scipy random state object (see ``scipy.random.RandomState``)
    """
    if (verbose):
        print("\n =================> INITIALIZATION: PEOPLE POSITIONS")
    ## To initialize people positions
    rng = np.random.RandomState()
    if (seed>0):
        rng = np.random.RandomState(seed)
    if (verbose):
        print("=================> INITIALIZATION: SEED = ",rng.get_state()[1][0])
    ## People properties: radius, velocity coefficient and initial random coordinates
    ## overlaps can occur...
    N = 0
    for gp in groups:
        N += gp["nb"]
    xyrv = np.zeros((N,4))  # x y r v
    gpid = np.zeros((N,),dtype=int)
    dest = []
    pos = 0
    ## Loop over all the groups
    for igp,gp in enumerate(groups):
        nn = gp["nb"]
        if (verbose):
            print("=================> INITIALIZATION: "+str(nn)+" IN BOX "
              + str(gp["box"]))
        pp, pp_dest = create_people_in_box(nn, gp["box"],
                                  gp["destination"],
                                  gp["radius_distribution"],
                                  gp["velocity_distribution"],
                                  dom, rng, verbose=verbose)
        pp = remove_overlaps_in_box(dom, gp["box"], pp,
                                    pp_dest, dt, rng,
                                    dmin_people=dmin_people,
                                    dmin_walls=dmin_walls,
                                    itermax=itermax,
                                    projection_method=projection_method,
                                    verbose=verbose)
        xyrv[pos:pos+nn,0] = pp[:,0]
        xyrv[pos:pos+nn,1] = pp[:,1]
        xyrv[pos:pos+nn,2] = pp[:,2]
        xyrv[pos:pos+nn,3] = pp[:,3]
        gpid[pos:pos+nn] = igp
        dest += pp_dest
        pos += nn
    if (verbose):
        print("=================> INITIALIZATION: LAST STEP (remove the overlaps"
          + " between individuals in different boxes...)")
    xyrv = remove_overlaps_in_box(dom, [dom.xmin, dom.xmax, dom.ymin, dom.ymax],
                                    xyrv, dest, dt, rng,
                                    dmin_people=dmin_people,
                                    dmin_walls=dmin_walls,
                                    projection_method=projection_method,
                                    itermax=itermax,verbose=verbose)
    if (verbose):
        print("=================> INITIALIZATION: END ! \n")
    people = {}
    people["xyrv"] = xyrv
    people["Vd"] = np.zeros((xyrv.shape[0],2))
    people["U"] = np.zeros((xyrv.shape[0],2))
    people["Uold"] = np.zeros((xyrv.shape[0],2))
    people["destinations"] = np.array(dest, dtype="U20")
    people["paths"] = {}
    people["rng"] = rng
    people["id"] = np.char.add([dom.name+'_']*xyrv.shape[0],
                                np.arange(xyrv.shape[0]).astype('<U3'))
    people["last_id"] = dom.name+'_'+str(xyrv.shape[0]-1)
    return people


def find_duplicate_people(all_people, domains):
    """
    This function determines people who have to be duplicated
    due to their presence in transit boxes (in Destination object)

    Parameters
    ----------
    all_people: dict
        Has one dictionary ``"people"`` by domain which has at least:
         * ``"xyrv"``: people coordinates and radius
         * ``"destinations"``: destination for each individual
    domains: dict
        Contains all the Domain objets

    Returns
    -------
    virtual_people: dict
        all duplicated people
    """
    virtual_people = {}
    for name in domains:
        virtual_people[name] = {}
        virtual_people[name]["xyrv"] = None
        virtual_people[name]["Vd"] = None
        virtual_people[name]["U"] = None
        virtual_people[name]["Uold"] = None

    for idom,name in enumerate(domains):
        #print("===> Look at virtual people from domain ",name)
        dom = domains[name]
        for dest_name in dom.destinations:
            box = dom.destinations[dest_name].next_transit_box

            if box is not None:
                box = np.array(box)
                box = np.vstack((box,box[0,:]))

                next_domain_name = dom.destinations[dest_name].next_domain

                ## People in this domain "name" that must be replicated in
                ## the domain "next_domain_name"
                p = all_people[name]["xyrv"]
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
                    test += np.arccos((vx*wx+vy*wy)/np.sqrt(vx*vx+vy*vy)/np.sqrt(wx*wx+wy*wy))
                #inside = ind[np.where(np.abs(test-2*np.pi)<0.001)[0]]
                inside = np.where(np.abs(test-2*np.pi)<0.001)[0]
                #print("test=",test)
                # print("people inside: ",inside,
                #      " must be copy in domain ",next_domain_name)
                for key in ["xyrv","Vd","U","Uold"]:
                    try:
                        virtual_people[next_domain_name][key] = np.concatenate((
                            virtual_people[next_domain_name][key],
                            all_people[name][key][inside]))
                    except:
                        virtual_people[next_domain_name][key] \
                                = all_people[name][key][inside]

                ## To check if an individual appears two times...
                if (inside.shape[0]>0):
                    virtual_people[next_domain_name]["xyrv"], indices \
                        = np.unique(virtual_people[next_domain_name]["xyrv"],
                                    return_index=True, axis=0)
                    virtual_people[next_domain_name]["Vd"] = \
                        virtual_people[next_domain_name]["Vd"][indices]
                    virtual_people[next_domain_name]["U"] = \
                        virtual_people[next_domain_name]["U"][indices]
                    virtual_people[next_domain_name]["Uold"] = \
                        virtual_people[next_domain_name]["Uold"][indices]

                ## People in the domain "next_domain_name" that must be replicated
                ## in this domain "name"
                p = all_people[next_domain_name]["xyrv"]
                test = np.zeros(p.shape[0])
                for i in np.arange(box.shape[0]-1):
                    vx = box[i,0]-p[:,0]
                    vy = box[i,1]-p[:,1]
                    wx = box[i+1,0]-p[:,0]
                    wy = box[i+1,1]-p[:,1]
                    test += np.arccos((vx*wx+vy*wy)/np.sqrt(vx*vx+vy*vy)/np.sqrt(wx*wx+wy*wy))
                #print("test=",test)
                inside = np.where(np.abs(test-2*np.pi)<0.001)[0]
                # print("people inside: ",inside,
                #       " must be copy in domain ",name)
                for key in ["xyrv","Vd","U","Uold"]:
                    try:
                        virtual_people[name][key] = np.concatenate((
                            virtual_people[name][key],
                            all_people[next_domain_name][key][inside]))
                    except:
                        virtual_people[name][key] \
                            = all_people[next_domain_name][key][inside]

                ## To check if an individual appears two times...
                if (inside.shape[0]>0):
                    virtual_people[name]["xyrv"], indices \
                        = np.unique(virtual_people[name]["xyrv"],
                                    return_index=True, axis=0)
                    virtual_people[name]["Vd"] \
                        = virtual_people[name]["Vd"][indices]
                    virtual_people[name]["U"] \
                        = virtual_people[name]["U"][indices]
                    virtual_people[name]["Uold"] \
                        = virtual_people[name]["Uold"][indices]

    return virtual_people

def people_update_destination(all_people, domains, thld, box=None):

    """
    This function updates people destinations and domains

    Parameters
    ----------
    all_people: dict
        Has one dictionary ``"people"`` by domain which has at least:
         * ``"xyrv"``: people coordinates and radius
         * ``"destinations"``: destination for each individual
    domains: dict
        Contains all the Domain objets
    thld: float
        threshold value used to decide if an individual has reached its destination
    box: numpy array
        box used to remove people outside

    Returns
    -------
    all_people: dict
        all people in all domains updated
    """

    for idom,name in enumerate(domains):
        print("===> Remove people outside its domain ",name," or the given box")
        p = all_people[name]["xyrv"]
        dom = domains[name]
        # Test if there are people outside the domain or outside the given box
        # If it is the case, we remove these persons...
        if box is None:
            ## Remove people who are outside the domain
            S = (p[:,0]-p[:,2]<=dom.xmin+dom.pixel_size) \
                + (p[:,0]-p[:,2]>=dom.xmax-dom.pixel_size) \
                + (p[:,1]-p[:,2]<=dom.ymin+dom.pixel_size) \
                + (p[:,1]-p[:,2]>=dom.ymax-dom.pixel_size)
        else:
            ## Remove people who are outside the given box
            S = (p[:,0]-p[:,2]<=box[0]+dom.pixel_size) \
                + (p[:,0]-p[:,2]>=box[1]-dom.pixel_size) \
                + (p[:,1]-p[:,2]<=box[2]+dom.pixel_size) \
                + (p[:,1]-p[:,2]>=box[3]-dom.pixel_size)
        out = np.where(S==True)[0]
        if (out.shape[0]>0):
            for ik,key in enumerate(all_people[name]):
                print(ik,key)
                if (key is not "paths"):
                    try:
                        if (all_people[name][key].shape[0]==p.shape[0]):
                            all_people[name][key] = \
                                np.delete(all_people[name][key],out,0)
                    except:
                        pass
    for idom,name in enumerate(domains):
        print("===> Update people destinations for domain ",name)
        #print("    => people at ",name," must be updated: ",all_people[name])
        p = all_people[name]["xyrv"]
        dom = domains[name]
        destinations = all_people[name]["destinations"]
        # Test if people are close to their destinations
        I = sp.floor((p[:,1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
        J = sp.floor((p[:,0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
        all_people[name]["I"] = I
        all_people[name]["J"] = J

        ## Remove people with I or J outside
        testI = (I>=0)*(I<dom.height)
        testJ = (J>=0)*(J<dom.width)
        outI = np.where( testI==False )[0]
        outJ = np.where( testJ==False )[0]
        out = np.unique(np.concatenate((outI,outJ)))
        if (out.shape[0]>0):
            print("** WARNING: people_update_destination **")
            print("** WARNING: "+str(out.shape[0])+" persons are deleted **")
            print("** WARNING: because outside the domain... **")
            print("** WARNING: Are there too many people in the domain ? **")
            for ik,key in enumerate(all_people[name]):
                try:
                    if (all_people[name][key].shape[0]==p.shape[0]):
                        all_people[name][key] = np.delete(all_people[name][key],out,0)
                except:
                    pass
            destinations = all_people[name]["destinations"]
            p = all_people[name]["xyrv"]
            I = sp.floor((p[:,1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
            J = sp.floor((p[:,0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
            all_people[name]["I"] = I
            all_people[name]["J"] = J

        ## Remove people inside the walls
        in_wall = sp.where( dom.wall_mask[I,J]>0 )[0]
        if (in_wall.shape[0]>0):
            print("** WARNING: people_update_destination **")
            print("** WARNING: "+str(in_wall.shape[0])+" persons are deleted **")
            print("** WARNING: because inside the walls... **")
            print("** WARNING: Are there too many people in the domain ? **")
            for ik,key in enumerate(all_people[name]):
                try:
                    if (all_people[name][key].shape[0]==p.shape[0]):
                        all_people[name][key] = \
                            np.delete(all_people[name][key],in_wall,0)
                except:
                    pass
            destinations = all_people[name]["destinations"]
            p = all_people[name]["xyrv"]
            I = sp.floor((p[:,1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
            J = sp.floor((p[:,0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
            all_people[name]["I"] = I
            all_people[name]["J"] = J

        ## Keep all people paths
        for ip,pid in enumerate(all_people[name]["id"]):
            try:
                all_people[name]["paths"][pid] = \
                    np.hstack((all_people[name]["paths"][pid],p[ip,:2]))
            except:
                all_people[name]["paths"][pid] = p[ip,:2]
        #print("path: ",all_people[name]["paths"])

        ## Change destination or remove people who are at their final destination
        id_rm = None
        for id,dest_name in enumerate(np.unique(destinations)):
            ind = np.where(destinations==dest_name)[0]
            D = dom.destinations[dest_name].distance[I[ind],J[ind]]-p[ind,2]
            # Test if people are close to their destination
            ind2 = sp.where(D<=thld)[0]
            #print("     ",ind2," are close to their destination ",dest_name)
            if (ind2.shape[0]>0):
                if (dom.destinations[dest_name].next_destination is not None):
                    next_domain = dom.destinations[dest_name].next_domain
                    if (next_domain != name):
                        # In this case we have to copy these persons in their new
                        # domain, and then remove them
                        #print("     => change domain and destination...")
                        # update destination
                        all_people[name]["destinations"][ind[ind2]] \
                            = dom.destinations[dest_name].next_destination
                        #print(dom.destinations[dest_name].next_destination)
                        # move
                        all_people[next_domain]["xyrv"] = np.concatenate(
                                 (all_people[next_domain]["xyrv"],
                                  all_people[name]["xyrv"][ind[ind2]]))
                        all_people[next_domain]["destinations"] = np.concatenate(
                                 (all_people[next_domain]["destinations"],
                                  all_people[name]["destinations"][ind[ind2]]))
                        all_people[next_domain]["id"] = np.concatenate(
                                 (all_people[next_domain]["id"],
                                  all_people[name]["id"][ind[ind2]]))
                        all_people[next_domain]["Vd"] = np.concatenate(
                                 (all_people[next_domain]["Vd"],
                                  all_people[name]["Vd"][ind[ind2]]))
                        all_people[next_domain]["U"] = np.concatenate(
                                 (all_people[next_domain]["U"],
                                  all_people[name]["U"][ind[ind2]]))
                        all_people[next_domain]["Uold"] = np.concatenate(
                                 (all_people[next_domain]["Uold"],
                                  all_people[name]["Uold"][ind[ind2]]))
                        # remove
                        try:
                            id_rm = np.concatenate((id_rm,ind[ind2]))
                        except:
                            id_rm = ind[ind2]
                    else:
                        # People do not change domains, they juste change their
                        # destination
                        #print("     => change destination...")
                        all_people[name]["destinations"][ind[ind2]] \
                            = dom.destinations[dest_name].next_destination
                else: ## in this case, we will only remove these persons...
                    #print("     => remove...")
                    try:
                        id_rm = np.concatenate((id_rm,ind[ind2]))
                    except:
                        id_rm = ind[ind2]

        if (id_rm is not None):
            for ik,key in enumerate(all_people[name]):
                try:
                    if (all_people[name][key].shape[0]==p.shape[0]):
                        all_people[name][key] \
                            = np.delete(all_people[name][key],id_rm,0)
                except:
                    pass

    return all_people

def sensor(door, xy0, xy1, t0, t1):
    """
    Compute the number of entries/exits through a door as a pedestrian
    sensor could do

    Parameters
    ----------
    door: numpy array
        door coordinates ``[x0,y0,x1,y1]``
    t0: float
        time
    t1: float
        time
    xy0: numpy array
        people coordinates at time ``t0``
    xy1: numpy array
        people coordinates at time ``t1``

    Returns
    -------
    id: numpy array
        index of persons who go through the door
    p: numpy array
        coordinates of intersection points between the door and people trajectories
    io: numpy array
        the exit direction is the normal direction, 1 = exit, -1 = entry
    times: numpy array
        exit or entry times
    entries: int
        number of entries
    exits: int
        number of exits
    """
    #
    #   trajectories:
    #               xy0
    #                |
    #                 |
    #   door:     d0--p--d1
    #                   |
    #                   xy1
    d0 = np.empty(xy0.shape)
    d0[:,0] = door[0]
    d0[:,1] = door[1]
    d1 = np.empty(xy1.shape)
    d1[:,0] = door[2]
    d1[:,1] = door[3]
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
    norm_vdoor_2 = np.sum(vdoor*vdoor,axis=1)
    dot_vdoor_vp0 = np.sum(vdoor*vp0,axis=1)
    is_p_in_door = (dot_vdoor_vp0>=0)*(dot_vdoor_vp0<=norm_vdoor_2)
    # Test if the intersection point is on the person trajectory
    vpxy0 = np.atleast_2d(p - xy0)
    norm_vtraj_2 = np.sum(vtraj*vtraj,axis=1)
    dot_vtraj_vpxy0 = np.sum(vtraj*vpxy0,axis=1)
    is_p_in_traj = (dot_vtraj_vpxy0>=0)*(dot_vtraj_vpxy0<=norm_vtraj_2)
    # Keep only points on the door and on the trajectory
    is_p_intersect = is_p_in_door*is_p_in_traj
    id = np.where(is_p_intersect==True)[0]
    # Test if the direction is the output normal: (d1-d0)_y , (d1-d0)_x
    vn = np.empty(vdoor.shape)
    vn[:,0] = vdoor[:,1]
    vn[:,1] = -vdoor[:,0]
    dot_vn_vtraj = np.sum(vn*vtraj,axis=1)
    is_normal_dir = (dot_vn_vtraj>0)
    io = (is_normal_dir[id]==True)*1+(is_normal_dir[id]==False)*(-1)
    exits = np.sum(io==1)
    entries = np.sum(io==-1)
    # Compute the distance from the intersection point p to xy0
    norm_vpxy0_2= np.sqrt(np.sum(vpxy0 * vpxy0, axis=1))
    # Compute the distance from the intersection point p to xy1
    vpxy1 = np.atleast_2d(p - xy1)
    norm_vpxy1_2 = np.sqrt(np.sum(vpxy1 * vpxy1, axis=1))
    # Compute the intersection time
    norm_vtraj = np.sqrt(norm_vtraj_2)
    dt = t1-t0
    times = t0 + (is_normal_dir==True)*(norm_vpxy0_2*dt/norm_vtraj) + \
            (is_normal_dir==False)*(norm_vpxy1_2*dt/norm_vtraj)
    return id, p[id,:], io, times[id], entries, exits


def plot_people(ifig, dom, people, contacts, colors, time=-1, axis=None,
        virtual_people=None, plot_people=True, plot_contacts=True,
        plot_velocities=False, plot_desired_velocities=False,
        plot_paths=False, plot_sensors=False,
        sensors=[], savefig=False, filename='fig.png', dpi = 150,
        cmap='winter'):
    """
    This function draws spheres for the individuals, \
    lines for the active contacts and arrows for the \
    (desired or real) velocities.

    Parameters
    ----------
    ifig: int
        figure number
    dom: Domain
        contains everything for managing the domain
    people: dict
        contains everything concerning people: ``x, y, r, v, U,...``
    contacts: numpy array
        all the contacts: ``i, j, dij, eij_x, eij_y``
    colors: numpy array
        scalar field used to define people colors
    time: float
        time in seconds
    virtual_people: dict
        contains everything concerning virtual people: ``x, y, r,...``
    axis: numpy array
        matplotlib axis: ``[xmin, xmax, ymin, ymax]``
    plot_people: boolean
        draws spheres for people if true
    plot_paths: boolean
        draws people paths if true
    plot_sensors: boolean
        draws sensor lines if true
    sensors: numpy array
        sensor line coordinates (see also the sensor function below)
    savefig: boolean
        writes the figure as a png file if true
    filename: string
        png filename used to write the figure
    dpi: integer
        number of pixel per inch for the saved figure
    cmap: string
        matplotlib colormap name
    """
    fig = plt.figure(ifig)
    plt.clf()
    ax1 = fig.add_subplot(111)
    # Domain
    ax1.imshow(dom.image,interpolation='nearest',
               extent=[dom.xmin,dom.xmax,dom.ymin,dom.ymax], origin='lower')
    if (plot_people):
        try:
            # People
            offsets = people["xyrv"][:,:2]
            ec = EllipseCollection(widths=2*people["xyrv"][:,2],
                                   heights=2*people["xyrv"][:,2],
                                   angles=0, units='xy',
                                   cmap=plt.get_cmap(cmap),
                                   offsets=offsets, transOffset=ax1.transData)
            ec.set_array(colors)
            ax1.add_collection(ec)
        except:
            pass
        try:
            # Virtual people
            offsets = virtual_people["xyrv"][:,:2]
            ecv = EllipseCollection(widths=2*virtual_people["xyrv"][:,2],
                                   heights=2*virtual_people["xyrv"][:,2],
                                   angles=0, units='xy',
                                   facecolors='black',
                                   alpha=0.5,
                                   offsets=offsets, transOffset=ax1.transData)
            ax1.add_collection(ecv)
        except:
            pass
    if (plot_contacts):
        try:
            if (virtual_people["xyrv"] is not None):
                xyrv = np.concatenate((people["xyrv"],virtual_people["xyrv"]))
            else:
                xyrv = people["xyrv"]
            # Contacts
            Nc = contacts.shape[0]
            if (Nc>0):
                for ic in np.arange(Nc):
                    i = np.int64(contacts[ic,0]); j = np.int64(contacts[ic,1])
                    if (j!=-1):
                        line = Line2D([ xyrv[i,0], xyrv[j,0] ],
                                      [ xyrv[i,1], xyrv[j,1] ],
                                      lw=1., alpha=0.6, color='k')
                    else:
                        line = Line2D([ xyrv[i,0],
                                        xyrv[i,0]
                                        - (xyrv[i,2]
                                        + contacts[ic,2])*contacts[ic,3] ],
                                      [ xyrv[i,1],
                                        xyrv[i,1]
                                        - (xyrv[i,2]
                                        + contacts[ic,2])*contacts[ic,4] ],
                                      lw=1., alpha=0.6, color='k')
                    ax1.add_line(line)
        except:
            pass
    if (plot_velocities):
        try:
            # Velocities
            Np = people["xyrv"].shape[0]
            for ip in np.arange(Np):
                arrow = Arrow( people["xyrv"][ip,0],
                               people["xyrv"][ip,1],
                               people["U"][ip,0],
                               people["U"][ip,1],
                               width=0.3, color="red")
                ax1.add_patch(arrow)
            if (virtual_people["xyrv"] is not None):
                Np = virtual_people["xyrv"].shape[0]
                for ip in np.arange(Np):
                    arrow = Arrow( virtual_people["xyrv"][ip,0],
                                   virtual_people["xyrv"][ip,1],
                                   virtual_people["U"][ip,0],
                                   virtual_people["U"][ip,1],
                                   width=0.3, color="red",alpha=0.6)
                    ax1.add_patch(arrow)
        except:
            pass
    if (plot_desired_velocities):
        try:
            # Velocities
            Np = people["xyrv"].shape[0]
            for ip in np.arange(Np):
                arrow = Arrow( people["xyrv"][ip,0],
                               people["xyrv"][ip,1],
                               people["Vd"][ip,0],
                               people["Vd"][ip,1],
                               width=0.3, color="blue")
                ax1.add_patch(arrow)
            if (virtual_people["xyrv"] is not None):
                Np = virtual_people["xyrv"].shape[0]
                for ip in np.arange(Np):
                    arrow = Arrow( virtual_people["xyrv"][ip,0],
                                   virtual_people["xyrv"][ip,1],
                                   virtual_people["Vd"][ip,0],
                                   virtual_people["Vd"][ip,1],
                                   width=0.3, color="blue",alpha=0.6)
                    ax1.add_patch(arrow)
        except:
            pass
    if (plot_paths):
        pathlines = []
        for ip,pid in enumerate(people["paths"]):
            pathlines.append(np.reshape(people["paths"][pid],(-1,2)))
        pathlinecoll = LineCollection(pathlines, color="black",
                       linewidths=0.5, linestyle="solid", alpha=0.5,
                       cmap=plt.get_cmap(cmap))
        #pathlinecoll.set_array(colors)
        ax1.add_collection(pathlinecoll)
    if (plot_sensors):
        for ss in sensors:
            line = Line2D( ss["line"][::2], ss["line"][1::2],
                           lw=1., alpha=0.6, color='g')
            ax1.add_line(line)
    if (axis):
        ax1.set_xlim(axis[0],axis[1])
        ax1.set_ylim(axis[2],axis[3])
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.axis('off')
    if (time>=0):
        ax1.set_title('time = {:.2F}'.format(time)+' s')
    fig.canvas.draw()
    if (savefig):
        fig.savefig(filename,dpi=dpi,bbox_inches='tight',pad_inches=0)


def plot_sensors(ifig, sensors, time, flux_timestep=1, savefig=False,
        filename='fig.png', dpi=150, cmap='winter'):
    """

    This function traces, for each sensor, a graph with the passage times
    of the individuals crossing the counting line, as well as the incoming
    and outgoing flows.

    Parameters
    ----------

    ifig: int
        figure number
    sensors: dict
        contains for each sensor, the date, position and direction of an
        individual who cut the counting line
    time: float
        time in seconds
    flux_timestep: float
        timestep for the fluxes: number of persons per flux_timestep seconds
    savefig: boolean
        writes the figure as a png file if true
    filename: string
        png filename used to write the figure
    dpi: integer
        number of pixel per inch for the saved figure
    cmap: string
        matplotlib colormap name
    """
    tmin = 0
    tmax = time

    Ns = len(sensors)
    if (Ns>0):
        fig = plt.figure(ifig)
        # ncol = np.floor(np.sqrt(Ns))
        # nrow = np.ceil(Ns / ncol)
        ncol = 2
        nrow = Ns
        for i,s in enumerate(sensors):
            name = s["name"] # sensor name
            tt = np.array(s["times"])
            dir = np.array(s["dir"])
            Np = tt.shape[0]
            ax1 = fig.add_subplot(nrow , ncol, 2*i+1)
            ax1.plot(np.arange(Np),tt,'b+')
            ax1.set_title(name+": crossing times (s)", fontsize="small")

            tgrid = np.arange(tmin, tmax, step=flux_timestep)
            tgrid = np.append(tgrid, tgrid[-1]+flux_timestep)
            flux_exits = np.zeros(tgrid.shape)
            flux_entries = np.zeros(tgrid.shape)
            exits = np.where(dir==1)[0]
            entries = np.where(dir==-1)[0]
            t_exits = np.ceil((tt[exits]-tmin)/flux_timestep)
            t_entries = np.ceil((tt[entries]-tmin)/flux_timestep)
            #t_exits = np.floor((tt[exits]-tmin)/flux_timestep)
            #t_entries = np.floor((tt[entries]-tmin)/flux_timestep)
            unique_exits, counts_exits = np.unique(t_exits, return_counts=True)
            unique_entries, counts_entries = np.unique(t_entries, return_counts=True)
            flux_exits[unique_exits.astype(int)] = counts_exits
            flux_entries[unique_entries.astype(int)] = counts_entries
            ax2 = fig.add_subplot(nrow , ncol, 2*i+2)
            ax2.plot(tgrid, flux_entries,':og',
                     tgrid, flux_exits, ':or')
            ax2.set_title(name + ": entries (green) and exits (red) per "
                + str(flux_timestep)+"s",fontsize="small")
        fig.set_tight_layout(True)
        fig.canvas.draw()
        if (savefig):
            fig.savefig(filename,dpi=dpi)
