# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
#
#     cromosim/examples/follow_the_leader/follow_the_leader.py
#     python follow_the_leader.py --json input_ftl_order1.json
#     python follow_the_leader.py --json input_ftl_order1_periodic.json
#     python follow_the_leader.py --json input_ftl_order2.json
#     python follow_the_leader.py --json input_ftl_order2_periodic.json
#
# License: GPL

import sys
import os
from optparse import OptionParser
import json
import numpy as np
import matplotlib.pyplot as plt

from cromosim.ftl import update_positions_ftl_order_1, update_positions_ftl_order_2
from cromosim.ftl import plot_people
plt.ion()

"""
 Follow The Leader model
 Models :
 - "ftl_order1"
 - "ftl_order2" with inertia
 Domain :
 - periodic
 - non-periodic
"""

"""
 Read the input JSON file
"""

parser = OptionParser(usage="usage: %prog [options] filename", version="%prog 1.0")
parser.add_option('--json', dest="jsonfilename",
                  default="input.json",
                  type="string",
                  action="store",
                  help="Input json filename")
opt, remainder = parser.parse_args()
with open(opt.jsonfilename) as json_file:
    input = json.load(json_file)

"""
 Parameters obtained from the json file
"""

# Prefix for the result path
model = input["model"]
# tau :
if (model == "ftl_order_2"):
    tau = input["tau"]
prefix = input["prefix"]
if not os.path.exists(prefix):
    os.makedirs(prefix)
# Number of persons : 0, 1, ..., N-1
N = input["N"]
# Length of the way [0,L]
L = input["L"]
# Periodicity : True or False
# If True : no Leader
# If False : the last person (number N-1) is the leader
periodicity = input["periodicity"]
# Final time
T = input["T"]
# Timestep
dt = input["dt"]
# The positions of the persons are drawn every "drawper" iterations
drawper = input["drawper"]
# Prescribed velocity for the Leader (only used in the non-periodic case)
if (periodicity is False):
    V_leader = lambda t: eval(input["V_leader"])
# Speed function
Phi = lambda w: eval(input["speed_function"])
# To initialize the positions of the persons (regularly):
xmin_t0 = input["xmin_t0"]
xmax_t0 = input["xmax_t0"]

# To build tgrid, Nt is the number of time iterations
Nt = int(np.floor(T/dt))+1
Nt += (Nt*dt < T)
tgrid = dt*np.arange(Nt)
tgrid[-1] = min(T, tgrid[-1])

# data : array where the positions and the velocities will be stored
data = np.zeros((N, 2, Nt))

# Iteration counter
counter = 0

# Initialization
time = 0.0
Xold = np.linspace(xmin_t0, xmax_t0, N)
Vold = np.zeros(Xold.shape)
data[:, 0, counter] = Xold
shift = np.zeros(Xold.shape)


# Time loop
while (time < T-0.5*dt):
    dt = min(dt, T-time)
    if (model == "ftl_order_1"):
        if periodicity:
            X, V = update_positions_ftl_order_1(L, Xold, time, dt, Phi)
        else:
            X, V = update_positions_ftl_order_1(L, Xold, time, dt, Phi,
                                                V_leader=V_leader, periodicity=periodicity)
    elif (model == "ftl_order_2"):
        if periodicity:
            X, V = update_positions_ftl_order_2(L, Xold, Vold, time, dt,
                                                tau, Phi)
        else:
            X, V = update_positions_ftl_order_2(L, Xold, Vold, time, dt,
                                                tau, Phi, V_leader=V_leader,
                                                periodicity=periodicity)
    else:
        print("Bad model name... EXIT")
        sys.exit()
    ind = np.where((X-Xold < 0))[0]
    shift[ind] += L
    data[:, 0, counter] = X + shift
    data[:, 1, counter] = V
    time += dt
    Xold = X
    Vold = V
    counter += 1
    if (counter % drawper == 0):
        if periodicity:
            plot_people(L, X, time, data, tgrid, speed_fct=Phi,
                        ifig=10, savefig=True,
                        filename=prefix+"fig_"+str(counter).zfill(6)+".png")
        else:
            plot_people(L, X, time, data, tgrid, V_leader=V_leader,
                        periodicity=periodicity, speed_fct=Phi, ifig=10,
                        savefig=True,
                        filename=prefix+"fig_"+str(counter).zfill(6)+".png")
        plt.pause(0.05)
        print("==> Time : ", time, " s, dt = ", dt, " counter = ", counter)

plt.ioff()
plt.show()
