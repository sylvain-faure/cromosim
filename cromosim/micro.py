# Authors:
#     Sylvain Faure <sylvain.faure@math.u-psud.fr>
#     Bertrand Maury <bertrand.maury@math.u-psud.fr>
# License: GPL

import sys
import scipy as sp
import numpy as np
try:
    from scipy.misc import imread
except:
    from imageio import imread
#from scipy.spatial import KDTree
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle, Rectangle, Polygon, Arrow
from matplotlib.lines import Line2D
from matplotlib.collections import EllipseCollection, LineCollection


def compute_contacts(dom, people, dmax):
    """
    This function uses a KDTree method to find the contacts \
    between individuals. Moreover the contacts with the walls \
    are also determined from the wall distance (obtained by the \
    fast-marching method).

    Parameters
    ----------
    dom: Domain
        contains everything for managing the domain
    people: numpy array
        people coordinates and radius : x,y,r
    dmax: float
        threshold value used to consider a contact as \
        active (dij<dmax)

    Returns
    -------
    contacts: numpy array
        all the contacts i,j,dij,eij_x,eij_y such that dij<dmax \
        and i<j (no duplication)
    """
    # lf : the number of points at which the algorithm
    # switches over to brute-force. Has to be positive.
    lf = 100
    if (lf>sys.getrecursionlimit()):
        sys.setrecursionlimit(lf)
    kd = cKDTree(people[:,:2],leafsize = lf)
    ## Find all pairs of points whose distance is at most dmax+2*rmax
    rmax = people[:,2].max()
    neighbors = kd.query_ball_tree(kd,dmax+2*rmax)
    ## Create the contact array : i,j,dij,eij_x,eij_y
    first_elements = np.arange(people.shape[0]) ## i.e. i
    other_elements = list(map(lambda x: x[1:], neighbors)) ## i.e. all the j values for each i
    lengths = list(map(len, other_elements))
    tt = np.stack([first_elements,lengths],axis=1)
    I = np.concatenate(list(map(lambda x: np.full((x[1],), x[0]), tt))).astype(int)
    J = np.concatenate(other_elements).astype(int)
    ind = np.where(I<J)[0]
    I = I[ind] ; J = J[ind]
    DP = people[J,:2]-people[I,:2]
    Norm = np.linalg.norm(DP,axis=1,ord=2)
    Dij = Norm - people[I,2]-people[J,2]
    ind = np.where(Dij<dmax)[0]
    Dij = Dij[ind]
    I = I[ind]
    J = J[ind]
    Norm = Norm[ind]
    DP = DP[ind]
    contacts = np.stack([I,J,Dij,DP[:,0]/Norm,DP[:,1]/Norm],axis=1)
    # Add contacts with the walls
    II = np.floor((people[:,1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
    JJ = np.floor((people[:,0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
    DD = dom.wall_distance[II,JJ] - people[:,2]
    ind = np.where(DD<dmax)[0]
    wall_contacts = np.stack([ind,-1*np.ones(ind.shape),DD[ind],
                              dom.wall_grad_X[II[ind],JJ[ind]],
                              dom.wall_grad_Y[II[ind],JJ[ind]]  ],axis=1)
    contacts = np.vstack([contacts,wall_contacts])
    return np.array(contacts)

def compute_desired_velocity(dom, people):
    """
    This function determines people desired velocities from the desired \
    velocity array computed by Domain thanks to a fast-marching method.

    Parameters
    ----------
    dom: Domain
        contains everything for managing the domain
    people: numpy array
        people coordinates and radius : x,y,r

    Returns
    -------
    I : numpy array
        people index i
    J : numpy array
        people index j
    Vd : numpy array
        people desired velocity
    """
    I = np.floor((people[:,1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
    J = np.floor((people[:,0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
    Vd = np.zeros( (people.shape[0],2) )
    Vd[:,0] = dom.desired_velocity_X[I,J]
    Vd[:,1] = dom.desired_velocity_Y[I,J]
    return I,J,Vd

def compute_forces(F, Fwall, people, contacts, U, Vd, lambda_, delta, k, eta):
    """
    This function computes all the forces (isentropic interaction and \
    friction) and sums them. The correcting pre-factor due to the vision \
    angle is also used into the social force term.

    Parameters
    ----------
    F: float
        social trend of an individual to keep apart from \
        another (homogeneous to a force)
    Fwall: float
        social trend of an individual to keep apart from \
        a wall (homogeneous to a force)
    people: numpy array
        people coordinates and radius : x,y,r
    contacts: numpy array
        all the contacts : i,j,dij,eij_x,eij_y
    U: numpy array
        people velocities
    Vd: numpy array
        people desired velocities
    lambda_: float
        quantifies the directional dependence when the vision angle \
        is considered (between [0,1], if equal to 1 the fully isotropic \
        case is recovered)
    delta: float
        maintains a certain distance between neighbors
    k: float
        used when there is overlapping, k is a stiffness constant \
        of individuals seen as deformable bodies
    eta: float
        friction coefficient

    Returns
    -------
    Forces: numpy array
        sum of all forces for each individual
    """
    Np = people.shape[0]
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
            fij_friction = eta*dij_moins*( -(U[i,0]-U[j,0])*eij_y + \
                (U[i,1]-U[j,1])*eij_x )
            Forces[i,0] -= fij_friction*eij_y
            Forces[i,1] += fij_friction*eij_x
            Forces[j,0] += fij_friction*eij_y
            Forces[j,1] -= fij_friction*eij_x
        else: ## contact person/walls
            fij = -Fwall*np.exp(-dij/delta) + k*dij_moins
            Forces[i,0] -= fij*eij_x
            Forces[i,1] -= fij*eij_y
    return Forces


def projection(dt, people, contacts, Vd, dmin = 0.0, \
               nb_iter_max = 100000, rho=0.1, tol = 0.01, log=False, method="cvxopt"):
    """
    From the desired velocities Vd, this projection step consists of computing \
    the global velocity field defined as the closest velocity to the \
    desired one among all the feasible fields (i.e. fields which do not lead \
    to disks overlapping).

    Parameters
    ----------
    dt: float
        time step
    people: numpy array
        people coordinates and radius : x,y,r
    contacts: numpy array
        all the contacts : i,j,dij,eij_x,eij_y
    Vd: numpy array
        people desired velocities
    dmin: float
        minimum distance guaranteed between individuals
    nb_iter_max: integer
        maximum number of iterations allowed
    rho: float
        parameter of the Uzawa method
    tol: float
        tolerance wished
    log: boolean
        to print the final accuracy, number of iterations,...
    method: string
        optimization algorithm : 'cvxopt' (default) or 'uzawa' (or 'mosek' if installed \
        with a license file).

    Returns
    -------
    B: numpy array
        constraint matrix
    U: numpy array
        new people velocities ensuring that there is no overlap \
        between individuals
    L: numpy array
        Lagrange multipliers (only when method='uzawa', None otherwise)
    P: numpy array
        pressure on each individual (only when method='uzawa', None otherwise)
    info: integer
        number of iterations needed
    """
    Np = people.shape[0]
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
            Z = (contacts[:,2]-dmin)/dt ## ie Dij/dt
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
            U = solution['x']
            if log:
                C = Z - B@U
                if (method == "mosek"):
                    print("    projection (mosek) : nb of contacts = ",Nc,
                          ", contrainte (Z-B@U).min() = ",C.min())
                else:
                    print("    projection (cvxopt) : nb of contacts = ",Nc,
                          ", nb of iterations = ",solution["iterations"],
                          ", status = ",solution["status"],
                          ", contrainte (Z-B@U).min() = ",C.min())
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
            while (( dt*R.max()>tol*2*people[:,2].min()) and (k<nb_iter_max)):
                U[:] = V[:] - B.transpose()@L[:]
                R[:] = B@U[:] - (D[:]-dmin)/dt
                L[:] = np.maximum(L[:] + rho*R[:], 0)
                k += 1
            P = np.zeros(Np) ## Pressure
            P[II[Jpos]] +=  3/(4*np.pi*people[II[Jpos],2]**2)*L[Jpos]
            P[JJ[Jpos]] +=  3/(4*np.pi*people[JJ[Jpos],2]**2)*L[Jpos]
            P[II[Jneg]] +=  3/(4*np.pi*people[II[Jneg],2]**2)*L[Jneg]
            if log:
                print("    projection (uzawa) : nb of contacts = ",Nc,
                      ", nb of iterations = ",k,", min = ",R.min(),
                      ", max = ",R.max(),", tol = ",tol)
            if (k==nb_iter_max):
                print("** WARNING : Method projection **")
                print("** WARNING : you have reached the maximum number of iterations,")
                print("** WARNING : it remains unsatisfied constraints !! ")
                info = -1
            else:
                info = k

        return info, B, U.reshape((Np, 2)), L, P


def move_people(time, dt, people, U, crosslines=[]):
    """
    Updates the people positions according to the new velocities U. \
    If there exists crosslines (i.e. sensors), computes also the id, time, direction and \
    impact points for the individuals who cross the lines.

    Parameters
    ----------
    time: float
        current time
    dt: float
        time step
    people: numpy array
        people coordinates and radius : x,y,r
    U: numpy array
        people velocities
    crosslines: list of numpy array
        list of lines [[x0,y0,x1,y1], [x0,y0,x1,y1], ...]

    Returns
    -------
    people: numpy array
        new people positions after moving
    io_id: list of integers
        indices of people who cross the lines (sensors)
    io_times: list of floats
        times when people cross the lines (sensors)
    io_pts: list of floats
        impact points for people crossing the lines (sensors)
    io_dir: list of floats
        directions for people crossing the lines (sensors), i.e. entry or exit
    """
    if (len(crosslines)>0):
        people_old = people.copy()
    ## Update people coordinates
    people[:,0] += dt*U[:,0]
    people[:,1] += dt*U[:,1]
    io_id = []
    io_times = []
    io_pts = []
    io_dir = []
    for l in crosslines:
        # line [x0,y0,x1,y1],
        id, pts, dir, times, entries, exits = sensor(l, people_old[:,0:2],
                    people[:,0:2], time, time+dt)
        io_id.append(id)
        io_times.append(times)
        io_pts.append(pts)
        io_dir.append(dir)
    if (len(crosslines)>0):
        return people, io_id, io_times, io_dir, io_pts
    else:
        return people


def exit_door(sexit, dom, people, U, arrays=[]):
    """
    Removes individuals who are at a distance less than sexit \
    to the closest door

    Parameters
    ----------
    sexit: float
        threshold value used to remove individuals close to the exit
    dom: Domain
        contains everything for managing the domain
    people: numpy array
        people coordinates and radius : x,y,r
    U: numpy array
        people velocities
    arrays: list of numpy array
        other arrays to resize similarly as people and U arrays

    Returns
    -------
    people: numpy array
        new people positions (outside individuals had been removed)
    U: numpy array
        new poeple velocities (outside individuals had been removed)
    arrays: list of numpy array
        new array resized
    """
    I = np.floor((people[:,1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
    J = np.floor((people[:,0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
    D = dom.door_distance[I,J]
    ind = np.where(D>sexit)
    if (len(arrays)>0):
        return people[ind[0],:], U[ind[0],:], [ a[ind[0]] for a in arrays]
    else:
        return people[ind[0],:], U[ind[0],:]


def exit_out_of_domain(dom, people, arrays=[], box=None):
    """
    Removes individuals who are outside the domain or outside a given box

    Parameters
    ----------
    dom: Domain
        contains everything for managing the domain
    people: numpy array
        people coordinates and radius : x,y,r
    arrays: list of numpy array
        other arrays to resize similarly as people and U
    box: numpy array
        box coordinates [xmin,xmax,ymin,ymax] which replace the \
        domain minimum and maximum coordinates

    Returns
    -------
    people: numpy array
        new people array (outside individuals had been removed)
    arrays: list of numpy array
        new arrays resized similarly as people array
    """
    if box is None:
        ## Remove people who are outside the domain
        S = (people[:,0]-people[:,2]<=dom.xmin+dom.pixel_size) + \
            (people[:,0]-people[:,2]>=dom.xmax-dom.pixel_size) + \
            (people[:,1]-people[:,2]<=dom.ymin+dom.pixel_size) + \
            (people[:,1]-people[:,2]>=dom.ymax-dom.pixel_size)
    else:
        ## Remove people who are outside the given box
        S = (people[:,0]-people[:,2]<=box[0]+dom.pixel_size) + \
            (people[:,0]-people[:,2]>=box[1]-dom.pixel_size) + \
            (people[:,1]-people[:,2]<=box[2]+dom.pixel_size) + \
            (people[:,1]-people[:,2]>=box[3]-dom.pixel_size)
    ind = np.where(S==False)[0]
    people = people[ind,:]
    if (len(arrays)>0):
        for a in arrays:
            a = a[ind]
    ## Remove people who are too close to walls or with a masked door distance
    I = np.floor((people[:,1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
    J = np.floor((people[:,0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
    Dwall = dom.wall_distance[I,J]-people[:,2]
    Ddoor = dom.door_distance[I,J]
    indDwall = np.where(Dwall<=dom.pixel_size)[0]
    indDdoor = np.where(Ddoor.mask==True)[0]
    ind = np.unique(np.concatenate((indDwall,indDdoor)))
    comp_ind = np.setdiff1d(np.arange(people.shape[0]), ind)
    if (len(arrays)>0):
        return people[comp_ind,:], [ a[comp_ind] for a in arrays]
    else:
        return people[comp_ind,:]


def exit_box(box, sexit, people, U, arrays=[]):
    """
    Removes individuals outside the box according to a threshold value.

    Parameters
    ----------
    box: numpy array
        box vertice coordinates [xmin, xmax, ymin, ymax]
    sexit: float
        threshold value used to remove individuals close to the exit
    people : numpy array
        people coordinates and radius : x,y,r
    U: numpy array
        people velocities
    arrays: list of numpy array
        other arrays to resize similarly as people and U

    Returns
    -------
    people: numpy array
        new people coordinates (outside individuals had been removed)
    U: numpy array
        new people velocities (outside individuals had been removed)
    arrays: list of numpy array
        new resized arrays taking into account deleted individuals
    """
    S = (people[:,0]-people[:,2]<=box[0]+sexit) + \
        (people[:,0]-people[:,2]>=box[1]-sexit) + \
        (people[:,1]-people[:,2]<=box[2]+sexit) + \
        (people[:,1]-people[:,2]>=box[3]-sexit)
    ind = np.where(S==False)
    if (len(arrays)>0):
        return people[ind[0],:], U[ind[0],:], [ a[ind[0]] for a in arrays]
    else:
        return people[ind[0],:], U[ind[0],:]


def periodic_bc_vertical(ymin, ymax, people, U, xmin=None, xmax=None, rng=None):
    """
    Does not exactly correspond to periodic boundary conditions (in the
    mathematical sense): the persons having an y-coordinate y less than ymin
    (respectively greater than ymax) are reinjected at y+(ymax-ymin) (respectively
    at y-(ymax-ymin)) with (optionally) random x-coordinates (between xmin and xmax,
    and with velocities equal to 0.

    Parameters
    ----------
    ymin: float
        minimal ordinate of the box
    ymax: float
        maximal ordinate of the box
    people: numpy array
        people coordinates and radius : x,y,r
    U: numpy array
        people velocities
    xmin: float
        minimal abscissa of the box
    xmax: float
        maximal abscissa of the box
    rng: RandomState
        scipy random state object (see scipy.random.RandomState)

    Returns
    -------
    people: numpy array
        new people coordinates (outside individuals had been moved)
    U: numpy array
        new people velocities
    """
    Smin = (people[:,1]<ymin)
    Smax = (people[:,1]>ymax)
    indmin = np.where(Smin==True)[0]
    indmax = np.where(Smax==True)[0]
    people[indmin,1] += (ymax-ymin)
    people[indmax,1] -= (ymax-ymin)
    U[indmin,:]=0
    U[indmax,:]=0
    if xmin is not None:
        if xmax is not None:
            if rng is None:
                rng = np.random.RandomState()
            x_indmin = rng.uniform(xmin, xmax, indmin.shape[0])
            people[indmin,0] = x_indmin
            x_indmax = rng.uniform(xmin, xmax, indmax.shape[0])
            people[indmax,0] = x_indmax
    return people, U


def create_people_in_box(nn, box, rmin, rmax, dom, rng):
    """
    To create nn persons in a given box. The overlaps are not treated for \
    the moment but one checks that the individuals are located in an area where \
    the desired velocity is well defined (outside inaccessible areas or \
    obstacles). If it is not the case, one changes their coordinates in consequence.

    Parameters
    ----------
    nn: integer
        number of individuals to create
    box: list
        coordinates of the box [xmin, xmax, ymin, ymax]
    rmin: float
        people minimal radius
    rmax: float
        people maximal radius
    dom: Domain
        contains everything for managing the domain
    rng: RandomState
        scipy random state object (see scipy.random.RandomState)

    Returns
    -------
    p: numpy array
        new people coordinates x y r
    """
    print("------ create_people_in_box --> Create "+str(nn)+ \
          " individuals in the box "+str(box)+", overlaps can occur...")
    p = np.zeros((nn,3))  # x y r
    p[:,2] = rng.uniform(rmin, rmax, nn)
    p_rmax = p[:,2].max()
    xmin, xmax, ymin, ymax = box
    p[:,0] = rng.uniform(xmin+p_rmax, xmax-p_rmax, nn)
    p[:,1] = rng.uniform(ymin+p_rmax, ymax-p_rmax, nn)
    ## To check if people are well localized in the domain
    ## i.e. where a desired velocity is defined
    while True:
        I, J, Vd = compute_desired_velocity(dom, p)
        normVd = Vd[:,0]**2+Vd[:,1]**2
        ind = np.where(normVd==0)[0]
        if (ind.shape[0]>0):
            p[ind,0] = rng.uniform(xmin+p_rmax, xmax-p_rmax, ind.shape[0])
            p[ind,1] = rng.uniform(ymin+p_rmax, ymax-p_rmax, ind.shape[0])
        else:
            break
    return p


def check_people_in_box(dom, box, p, rng):
    """
    To check that people coordinates are in the given box (test 1) and in an \
    usable space i.e. in an area accessible and not concerned by obstacles \
    (test 2). On the other hand, one moves the individuals which do not satisfy \
    these two tests.

    Parameters
    ----------
    dom: Domain
        contains everything for managing the domain
    box: list
        coordinates of the box [xmin, xmax, ymin, ymax]
    p: numpy array
        people coordinates x y r
    rng: RandomState
        scipy random state object (see scipy.random.RandomState)

    Returns
    -------
    p: numpy array
        new people coordinates x y r

    """
    print("------ check_people_in_box --> To verify that "+str(p.shape[0])+ \
          " individuals are in the domain, in the box and with a defined"+ \
          " desired velocity")
    p_rmax = p[:,2].max()
    xmin, xmax, ymin, ymax = box
    info = False
    while True:
        ## test 1
        I = np.floor((p[:,1]-dom.ymin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
        J = np.floor((p[:,0]-dom.xmin-0.5*dom.pixel_size)/dom.pixel_size).astype(int)
        test1 = (I>=0)*(I<dom.height)*(J>=0)*(J<dom.width)
        ind1 = np.where(test1==0)[0]
        if (ind1.shape[0]>0):
            print("------ check_people_in_box --> "+str(ind1.shape[0])+ \
                  " individuals outside the domain")
            info = True
            p[ind1,0] = rng.uniform(xmin+p_rmax, xmax-p_rmax, ind1.shape[0])
            p[ind1,1] = rng.uniform(ymin+p_rmax, ymax-p_rmax, ind1.shape[0])
        else:
            ## test 2
            I, J, Vd = compute_desired_velocity(dom, p)
            normVd = Vd[:,0]**2+Vd[:,1]**2
            test2 = (p[:,0]>xmin+p_rmax)*(p[:,0]<xmax-p_rmax) \
                  *(p[:,1]>ymin+p_rmax)*(p[:,1]<ymax-p_rmax) \
                  *(normVd>0)
            ind2 = np.where(test2==0)[0]
            if (ind2.shape[0]>0):
                print("------ check_people_in_box --> "+str(ind2.shape[0])+ \
                      " individuals with an undefined desired velocity ")
                info = True
                p[ind2,0] = rng.uniform(xmin+p_rmax, xmax-p_rmax, ind2.shape[0])
                p[ind2,1] = rng.uniform(ymin+p_rmax, ymax-p_rmax, ind2.shape[0])
            else:
                print("------ check_people_in_box --> OK !")
                break
    return info, p

def remove_overlaps_in_box(dom, box, p, dt, rng, dmin, itermax=10):
    """
    To remove the overlaps between individuals (spheres) in the given box. \
    Several Uzawa projections are used to give a better robustness to this \
    process when the number of overlaps is very high e.g. during the \
    initialization when the individuals have random positions.

    Parameters
    ----------
    dom: Domain
        contains everything for managing the domain
    box: list
        coordinates of the box [xmin, xmax, ymin, ymax]
    p: numpy array
        people coordinates x y r
    dt: float
        time step
    rng: RandomState
        scipy random state object (see scipy.random.RandomState)
    dmin: float
        minimal distance allowed between individuals
    itermax: integer
        maximal number of Uzawa projections (10 by default)

    Returns
    -------
    p: numpy array
        new people coordinates x y r

    """
    print("------ remove_overlaps_in_box --> Remove overlaps...")
    dmax_init = p[:,2].max()
    it1 = 0
    t = 0
    while (True):
        print("------ remove_overlaps_in_box --> Remove overlaps in box "+ \
              str(box)+": iteration "+str(it1)+" / "+str(itermax))
        c = compute_contacts(dom, p, dmax_init)
        I, J, Vd = compute_desired_velocity(dom, p)
        info1, B, U, L, P = projection(dt, p, c, Vd, dmin = dmin, \
                                             nb_iter_max = 10000, log = False)
        p = move_people(t, dt, p, U)
        it1 += 1
        info2, people = check_people_in_box(dom, box, p, rng)
        if ((info1>=0) and (info2==False)):
            print("------ remove_overlaps_in_box --> Successful... No overlaps")
            break
        if (it1>=itermax):
            print("------ remove_overlaps_in_box --> ** WARNING : Unsuccessful initialization **")
            print("------ remove_overlaps_in_box --> ** WARNING : Let us continue...")
            break
    return p


def people_initialization(N, init_people_box, dom, dt, rmin, rmax, dmin=0,
                          seed=0, itermax=10):
    """
    To initialize people array (xyr) with coordinates in several boxes and \
    without overlaps between them.

    Parameters
    ----------
    N: list
        list of the numbers of persons to add in each box
    init_people_box: list
        list of the boxes [[xmin, xmax, ymin, ymax],...]
    dom: Domain
        contains everything for managing the domain
    dt: float
        time step
    rmin: float
        people minimal radius
    rmax: float
        people maximal radius
    dmin: float
        minimal distance allowed between individuals (0 by default)
    seed: integer
        seed to initialize the RandomState (0 by default means different seed \
        at each run)
    itermax: integer
        maximal number of Uzawa projections (10 by default)

    Returns
    -------
    people: numpy array
        new people coordinates x y r
    people_init_box_id: numpy array
        box number for each individual
    rng: RandomState
        scipy random state object (see scipy.random.RandomState)
    """
    print("\n =================> INITIALIZATION : PEOPLE POSITIONS")
    ## To initialize people positions
    rng = np.random.RandomState()
    if (seed>0):
        rng = np.random.RandomState(seed)
    print("=================> INITIALIZATION : SEED = ",rng.get_state()[1][0])
    ## People properties : radius and initial random coordinates
    ## overlaps can occur...
    people = np.zeros((N.sum(),3))  # x y r
    people_init_box_id = np.zeros((N.sum(),),dtype=int)
    pos = 0
    ## Loop over all the boxes where we are supposed to put individuals at
    ## the initializations
    for ip,nn in enumerate(N):
        print("=================> INITIALIZATION : "+str(nn)+" IN BOX "+ \
              str(init_people_box[ip]))
        pp = create_people_in_box(nn, init_people_box[ip], rmin, rmax, dom, rng)
        pp = remove_overlaps_in_box(dom, init_people_box[ip], pp, dt, rng,
                                    dmin, itermax=itermax)
        people[pos:pos+nn,0] = pp[:,0]
        people[pos:pos+nn,1] = pp[:,1]
        people[pos:pos+nn,2] = pp[:,2]
        people_init_box_id[pos:pos+nn] = ip
        pos += nn
    print("=================> INITIALIZATION : LAST STEP (remove the overlaps"+ \
          " between individuals in different boxes...)")
    people = remove_overlaps_in_box(dom, [dom.xmin, dom.xmax, dom.ymin, dom.ymax],
                                    people, dt, rng, dmin, itermax=itermax)
    print("=================> INITIALIZATION : END ! \n")
    return people, people_init_box_id, rng


def sensor(door, xy0, xy1, t0, t1):
    """
    Compute the number of entries/exits through a door as a pedestrian
    sensor could do

    Parameters
    ----------
    door: numpy array
        door coordinates [x0,y0,x1,y1]
    t0: float
        time
    t1: float
        time
    xy0: numpy array
        people coordinates at time t0
    xy1: numpy array
        people coordinates at time t1

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
    #   trajectories :
    #               xy0
    #                |
    #                 |
    #   door :     d0--p--d1
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
    # Test if the direction is the output normal : (d1-d0)_y , (d1-d0)_x
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

def add_people_in_box(Np, dom, xmin, xmax, ymin, ymax, rmin, rmax, rng):
    """
    Adds new persons in the box [xmin,xmax]x[ymin,ymax]
    Be careful : overlaps can occur...

    Parameters
    ----------
    Np: int
        Number of persons
    dom: Domain
        contains everything for managing the domain
    xmin: float
        minimal abscissa of the box
    xmax : float
        maximal abscissa of the box
    ymin: float
        minimal ordinate of the box
    ymax: float
        maximal ordinate of the box
    rmin: float
        minimum radius for the individuals
    rmax: float
        maximum radius for the individuals
    rng: scipy.random.RandomState
        container for the Mersenne Twister pseudo-random number generator
    Returns
    -------
    people : numpy array
        people coordinates and radius
    """
    px = dom.pixel_size
    people = np.zeros((Np,3))  # x y r
    people[:,0] = rng.uniform(xmin, xmax, Np)
    people[:,1] = rng.uniform(ymin, ymax, Np)
    people[:,2] = rng.uniform(rmin, rmax, Np)
    I = np.floor((people[:,1]-dom.ymin-0.5*px)/px).astype(int)
    J = np.floor((people[:,0]-dom.xmin-0.5*px)/px).astype(int)

    D = dom.wall_distance[I,J]-people[:,2]
    ind = np.where(D<=px)
    while (ind[0].shape[0]>0):
        people[ind[0],0] = rng.uniform(xmin, xmax, ind[0].shape[0])
        people[ind[0],1] = rng.uniform(ymin, ymax, ind[0].shape[0])
        I = np.floor((people[:,1]-dom.ymin-0.5*px)/px).astype(int)
        J = np.floor((people[:,0]-dom.xmin-0.5*px)/px).astype(int)
        D = dom.wall_distance[I,J]-people[:,2]
        ind = np.where(D<=px)
    return people


def plot_people(ifig, dom, people, contacts, U, colors, time=-1, axis=None,
                plot_people=True, plot_contacts=True, plot_velocities=True,
                plot_paths=False,paths=None,plot_sensors=False,sensors=[],
                savefig=False, filename='fig.png', cmap='winter'):
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
    people: numpy array
        people coordinates and radius : x,y,r
    contacts: numpy array
        all the contacts : i,j,dij,eij_x,eij_y
    U: numpy array
        people velocities
    colors: numpy array
        scalar field used to define people colors
    time: float
        time in seconds
    axis: numpy array
        matplotlib axis : [xmin, xmax, ymin, ymax]
    plot_people: boolean
        draws spheres for people if true
    plot_paths: boolean
        draws people paths if true
    paths: numpy array
        coordinates of the people paths
    plot_sensors: boolean
        draws sensor lines if true
    sensors: numpy array
        sensor line coordinates (see also the sensor function below)
    savefig: boolean
        writes the figure as a png file if true
    filename: string
        png filename used to write the figure
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
        # People
        #offsets = list(zip(people[:,:2])) ## for older versions of matplotlib...
        offsets = people[:,:2]
        ec = EllipseCollection(widths=2*people[:,2], heights=2*people[:,2],
                               angles=0, units='xy', cmap=plt.get_cmap(cmap),
                               offsets=offsets, transOffset=ax1.transData)
        ec.set_array(colors)
        ax1.add_collection(ec)
    if (plot_contacts):
        # Contacts
        Nc = contacts.shape[0]
        if (Nc>0):
            for ic in np.arange(Nc):
                i = np.int64(contacts[ic,0]); j = np.int64(contacts[ic,1])
                if (j!=-1):
                    line = Line2D([ people[i,0],people[j,0] ],
                              [ people[i,1],people[j,1] ],
                              lw=1., alpha=0.6, color='k')
                else:
                    line = Line2D([ people[i,0],
                        people[i,0]-(people[i,2]+contacts[ic,2])*contacts[ic,3] ],
                        [ people[i,1],
                        people[i,1]-(people[i,2]+contacts[ic,2])*contacts[ic,4] ],
                              lw=1., alpha=0.6, color='k')
                ax1.add_line(line)
    if (plot_velocities):
        # Velocities
        Np = people.shape[0]
        for ip in np.arange(Np):
            arrow = Arrow(people[ip,0], people[ip,1], U[ip,0], U[ip,1], width=0.3)
            ax1.add_patch(arrow)
    if ((plot_paths) and (paths is not None)):
        mpaths = np.ma.masked_values(paths, 1e99)
        pathlines = []
        for id in np.arange(paths.shape[0]):
            pathlines.append(np.swapaxes(mpaths[id,:,:],0,1))
        pathlinecoll = LineCollection(pathlines, linewidths=0.5, linestyle="solid",
                       cmap=plt.get_cmap(cmap))
        pathlinecoll.set_array(colors)
        ax1.add_collection(pathlinecoll)
    if (plot_sensors):
        # Sensors
        for ss in sensors:
            line = Line2D([ss[0],ss[2]], [ss[1],ss[3]],lw=1., alpha=0.6, color='g')
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
        fig.savefig(filename,dpi=600,bbox_inches='tight',pad_inches=0)


def plot_sensor_data(ifig, sensor_data, time, initial_door_dist=None, axis = None, \
                     flux_timestep=1, \
                     savefig=False, filename='fig.png', cmap='winter'):
    """
    When a sensor line is defined this function allows to draw the \
    repartition of the people exit times.

    Parameters
    ----------

    ifig: int
        figure number
    sensor_data : numpy array
        [time, direction, intersection_point[2]] for each individual
    time: float
        time in seconds
    initial_door_dist: numpy array
        people initial distance to the door
    axis: numpy array
        matplotlib axis : [xmin, xmax, ymin, ymax]
    flux_timestep: float
        timestep for the fluxes : number of persons per flux_timestep seconds
    savefig: boolean
        writes the figure as a png file if true
    filename: string
        png filename used to write the figure
    cmap: string
        matplotlib colormap name
    """
    Np = sensor_data.shape[0]
    tmin = 0
    tmax = time

    fig = plt.figure(ifig)
    plt.clf()

    ax1 = fig.add_subplot(211)
    if (initial_door_dist is None):
        ax1.plot(np.arange(Np),sensor_data[:,0],'b+')
        ax1.set_title('Crossing time (s) vs people id')
    else:
        ax1.plot(initial_door_dist,sensor_data[:,0],'b+')
        ax1.set_title('Crossing time (s) vs initial door distance (m)')
    if (axis):
        ax1.set_xlim(axis[0],axis[1])
        ax1.set_ylim(axis[2],axis[3])
    #ax1.set_xticks([])
    #ax1.set_yticks([])
    #ax1.axis('off')

    tgrid = np.arange(tmin, tmax, step=flux_timestep)
    tgrid = np.append(tgrid, tgrid[-1]+flux_timestep)
    flux_exits = np.zeros(tgrid.shape)
    flux_entries = np.zeros(tgrid.shape)
    exits = np.where(sensor_data[:,1]==1)[0]
    entries = np.where(sensor_data[:,1]==-1)[0]
    t_exits = np.ceil((sensor_data[exits,0]-tmin)/flux_timestep)
    t_entries = np.ceil((sensor_data[entries,0]-tmin)/flux_timestep)
    #t_exits = np.floor((sensor_data[exits,0]-tmin)/flux_timestep)
    #t_entries = np.floor((sensor_data[entries,0]-tmin)/flux_timestep)
    unique_exits, counts_exits = np.unique(t_exits, return_counts=True)
    unique_entries, counts_entries = np.unique(t_entries, return_counts=True)
    flux_exits[unique_exits.astype(int)] = counts_exits
    flux_entries[unique_entries.astype(int)] = counts_entries

    ax2 = fig.add_subplot(212)
    ax2.plot(tgrid, flux_entries,':og',
             tgrid, flux_exits, ':or')
    ax2.set_title("Entries (green) and exits (red) per "+str(flux_timestep)+" s")
    if (axis):
        ax2.set_xlim(axis[0],axis[1])
        ax2.set_ylim(axis[2],axis[3])
    #ax2.set_xticks([])
    #ax2.set_yticks([])
    #ax2.axis('off')
    # Optionally : adds some histograms
    # if (exits.shape[0]>0):
    #     ax3 = fig.add_subplot(413)
    #     t_exits_sorted = np.sort(sensor_data[exits,0])
    #     #print("t_exits_sorted = ",t_exits_sorted)
    #     tmp = np.concatenate(([0],t_exits_sorted))
    #     bins = 0.5*(tmp[:-1]+tmp[1:])
    #     widths = tmp[1:]-tmp[:-1]
    #     heights = 1/widths
    #     ax3.bar(bins, heights, width=widths,color='r',align='center')
    #
    # if (entries.shape[0]>0):
    #     ax4 = fig.add_subplot(414)
    #     t_entries_sorted = np.sort(sensor_data[entries,0])
    #     tmp = np.concatenate(([0],t_entries_sorted))
    #     bins = 0.5*(tmp[:-1]+tmp[1:])
    #     widths = tmp[1:]-tmp[:-1]
    #     heights = 1/widths
    #     ax4.bar(bins, heights, width=widths,color='r',align='center')
    fig.set_tight_layout(True)
    fig.canvas.draw()
    if (savefig):
        fig.savefig(filename,dpi=300)
