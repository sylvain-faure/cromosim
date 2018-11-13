# Authors:
#     Sylvain Faure <sylvain.faure@math.u-psud.fr>
#     Bertrand Maury <bertrand.maury@math.u-psud.fr>
# License: GPL

import sys
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection

def update_positions_ftl_order_1(L, X, t, dt, Phi, periodicity=True,
                                 V_leader=None):
    """
    To update the positions according to the follow the leader model (order 1)

    Parameters
    ----------
    L: float
        width of the domain
    X: numpy array
        coordinates of the individuals
    t: float
        time
    dt: float
        time step
    Phi: function
        speed as a function of the distance
    periodicity: boolean
        if true the domain is periodic
    V_leader: function
        leader velocity as a function of the time

    Returns
    -------
    Xnew: numpy array
        new coordinates of the individuals
    V: numpy array
        new velocities of the individuals

    """
    if (periodicity == True):
        # No leader
        S = sp.argsort(X)
        W = X[S[1:]] - X[S[:-1]]
        W = sp.append(W,(X[S[0]]+L) - X[S[-1]])
        V = Phi(W)
        Xnew = X[S] + dt*V
        Xnew1 = Xnew*(Xnew<L) + (Xnew-L)*(Xnew>=L)
        Xnew[S] = Xnew1
        V1 = V ; V[S] = V1
    else:
        # The last position in the array is the leader...
        W = X[1:]-X[:-1]
        VV = Phi(W)
        VL = V_leader(t)
        V = sp.append(VV,VL)
        Xnew = X+dt*V
    return Xnew, V


def update_positions_ftl_order_2(L, X, U, t, dt, tau, Phi, periodicity=True,
                                 V_leader=None):
    """
    To update the positions according to the follow the leader model (order 2)

    Parameters
    ----------
    L: float
        width of the domain
    X: numpy array
        coordinates of the individuals
    U: numpy array
        velocities of the individuals
    t: float
        time
    dt: float
        time step
    tau: float
        the actual velocity of an individual relaxes toward the velocity that \
        is associated to the current distance from the person in front of them \
        with a characteristic time tau
    Phi: function
        speed as a function of the distance
    periodicity: boolean
        if true the domain is periodic
    V_leader: function
        leader velocity as a function of the time

    Returns
    -------
    Xnew: numpy array
        new coordinates of the individuals
    Unew: numpy array
        new velocities of the individuals

    """
    if (periodicity == True):
        # No leader
        S = sp.argsort(X)
        W = X[S[1:]]-X[S[:-1]]
        W = sp.append(W,(X[S[0]]+L)-X[S[-1]])
        V = Phi(W)
        Unew = U[S]+dt*(V-U[S])/tau
        Xnew = X[S]+dt*Unew
        Xnew1 = Xnew*(Xnew<L)+(Xnew-L)*(Xnew>=L)
        Xnew[S] = Xnew1
        V1 = V ; V[S] = V1
        Unew1 = Unew ; Unew[S] = Unew
    else:
        # The last position in the array is the leader...
        Xnew = X+dt*U
        W = X[1:]-X[:-1]
        Unew = sp.zeros(U.shape)
        Unew[:-1] = U[:-1] + dt*(Phi(W)-U[:-1])/tau
        Unew[-1] = V_leader(t)
    return Xnew, Unew


def plot_people(L, X, t, data, tgrid, periodicity=True, V_leader=None,
                speed_fct=None, savefig=False, filename='fig.png',
                ifig=100, spheresize=300):
    """
    To plot people and their individual paths

    Parameters
    ----------
    L: float
        width of the domain
    X: numpy array
        positions of the individuals
    data: numpy array
        coordinates of indviduals at each time
    tgrid: numpy array
        temporal discretization
    periodicity: boolean
        if true the domain is periodic
    V_leader: function
        leader velocity as a function of the time
    speed_fct: function
        speed as a function of the distance
    savefig: boolean
        writes the figure as a png file if true
    filename: string
        png filename used to write the figure
    ifig: int
        figure number
    spheresize: int
        size of the spheres used to represent the individuals
    """
    N = X.shape[0]
    XY = sp.zeros((N,2))
    XY[:,0] = X
    itmax = sp.where(tgrid<=t)[0].max()

    fig = plt.figure(ifig,figsize=(16,9))
    fig.clf()

    if periodicity:
        ax1 = fig.add_subplot(221)
    else:
        ax1 = fig.add_subplot(231)
    ax1.set_xlim(0,L)
    ax1.set_ylim(-0.5,0.5)
    #ax1.set_xticks([])
    ax1.set_yticks([])
    #ax1.axis('off')
    ax1.set_xlabel('x')
    offsets = list(zip(XY))
    colors = sp.arange(N)/N
    colors = sp.ones(N)#sp.arange(N)/N
    scale = 20
    sc = ax1.scatter(XY[:,0],XY[:,1], c=colors, s=spheresize, edgecolor='None',
        marker='.',cmap=plt.get_cmap('Greys'),norm=plt.Normalize(-0.2, 1))
    line = Line2D([ 0, 0 ],[ -0.05, 0.05 ], lw=3., alpha=1, c='k')
    ax1.add_line(line)
    line = Line2D([ L, L ],[ -0.05, 0.05 ], lw=3., alpha=1, c='k')
    ax1.add_line(line)

    if (periodicity==False):
        ax3 = fig.add_subplot(233)
        ax3.plot(tgrid,V_leader(tgrid))
        ax3.set_xlabel('time')
        ax3.set_ylabel('Leader velocity')
        ax3.plot(t,V_leader(t),'ko')

    wgrid = sp.linspace(0,L,500)
    if periodicity:
        ax2 = fig.add_subplot(222)
    else:
        ax2 = fig.add_subplot(232)
    ax2.plot(wgrid,speed_fct(wgrid))
    S = sp.argsort(X)
    W = X[S[1:]] - X[S[:-1]]
    if periodicity:
        W = sp.append(W,(X[S[0]]+L) - X[S[-1]])
    ax2.plot(W,speed_fct(W),'ko')
    ax2.set_xlabel('distance')
    ax2.set_ylabel('speed')

    if periodicity:
        ax4 = fig.add_subplot(223)
    else:
        ax4 = fig.add_subplot(234)
    for ip in sp.arange(N):
        ax4.plot(tgrid[:itmax], data[ip,0,:itmax])
    ax4.set_ylim(0,max(L,data[:,0,:itmax].max()))
    ax4.set_xlim(tgrid.min(),tgrid.max())
    ax4.set_xlabel('time')
    ax4.set_ylabel('x')

    if periodicity:
        ax5 = fig.add_subplot(224)
    else:
        ax5 = fig.add_subplot(235)
    Y = sp.linspace(0, 0.25*L*N, num=N)
    for ip in sp.arange(N-1):
        ax5.plot(tgrid[:itmax],Y[ip]+(data[ip+1,0,:itmax]-data[ip,0,:itmax]))
    if periodicity:
        d0 = data[0,0,:itmax]
        dN = data[N-1,0,:itmax]
        ax5.plot(tgrid[:itmax],Y[-1]+(d0<dN)*(d0+L-dN)+(d0>=dN)*(d0-dN))
    ax5.set_yticks([])
    #ax5.set_ylim(0,1.1*Y.max())
    ax5.set_xlim(tgrid.min(),tgrid.max())
    ax5.set_xlabel('time')
    ax5.set_ylabel('distances')

    fig.set_tight_layout(True)

    if (savefig):
        fig.savefig(filename,dpi=150)
