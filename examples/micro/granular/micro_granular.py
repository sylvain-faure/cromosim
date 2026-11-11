# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
#
#     cromosim/examples/micro/granular/micro_granular.py
#     python micro_granular.py --json input.json
#
# License: GPL

import sys
import os
import json
import time
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Rectangle, Polygon
from matplotlib.lines import Line2D

import logging
logging.basicConfig(
    level=logging.INFO,
    # format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
    format="%(levelname)-8s | %(name)s | %(message)s",
)
LOGGER = logging.getLogger(__name__)
# LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"

from cromosim.domain import Destination, Domain
from cromosim.micro import (
    people_initialization,
    add_people_flows,
    append_people,
    compute_contacts,
    find_duplicate_people,
    move_people,
    people_update_destination,
    projection,
)
from cromosim.visu import (
    plot_sensors
)

# plt.ion()

"""
    python3 micro_granular.py --json input.json
"""

parser = OptionParser(usage="usage: %prog [options] filename",
                      version="%prog 1.0")
parser.add_option('--json', dest="jsonfilename", default="input.json",
                  type="string", action="store", help="Input json filename")
opt, remainder = parser.parse_args()

LOGGER.info(f"Loading JSON configuration: {opt.jsonfilename}")

with open(opt.jsonfilename) as json_file:
    try:
        config = json.load(json_file)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {opt.jsonfilename}: {exc}") from exc

LOGGER.info(
    "Simulation: Tf=%s, dt=%s, projection_method=%s",
    config["Tf"],
    config["dt"],
    config["projection_method"],
)
LOGGER.info(
    "Contacts: dmax=%s, dmin_people=%s, dmin_walls=%s",
    config["dmax"],
    config["dmin_people"],
    config["dmin_walls"],
)
LOGGER.info(
    "Visualization: enabled=%s, backend=%s, people=%s " 
    "contacts=%s velocities=%s desired_velocities=%s "
    "paths=%s sensors=%s plot_period=%s",
    config["with_graphes"],
    config["plot_module"],
    config["plot_people"],
    config["plot_contacts"],
    config["plot_velocities"],
    config["plot_desired_velocities"],
    config["plot_paths"],
    config["plot_sensors"],
    config["plot_period"],
)
LOGGER.info(
    "Miscellaneous: prefix=%s, seed=%s, dmin_walls=%s",
    config["prefix"],
    config["seed"],
    config["dmin_walls"],
)

os.makedirs(config["prefix"], exist_ok=True)

"""
    Build the Domain objects
"""

domains = {}

domain_configs = config["domains"]

for domain_index, domain_config in enumerate(domain_configs):
    name = domain_config["name"]
    LOGGER.info(f"Start building domain {domain_index}: {name}")
    if domain_config["background"] == "":
        domain = Domain(
            name=name,
            pixel_size=domain_config["px"],
            width=domain_config["width"],
            height=domain_config["height"],
            wall_colors=domain_config["wall_colors"],
        )
    else:
        domain = Domain(
            name=name,
            background=domain_config["background"],
            pixel_size=domain_config["px"],
            wall_colors=domain_config["wall_colors"],
        )
        
    LOGGER.info(
        f"Building domain {domain_index}: {name} => "
        f"Add geometric objects described in the JSON "
        f"configuration")
    
    for shape_config in domain_config["shape_lines"]:
        shape = Line2D(
            shape_config["xx"], shape_config["yy"],
            linewidth=shape_config["linewidth"],
        )
        domain.add_shape(
            shape,
            outline_color=shape_config["outline_color"],
            fill_color=shape_config["fill_color"],
        )

    for shape_config in domain_config["shape_circles"]:
        shape = Circle(
            (shape_config["center_x"], shape_config["center_y"]),
            shape_config["radius"],
        )
        domain.add_shape(
            shape,
            outline_color=shape_config["outline_color"],
            fill_color=shape_config["fill_color"],
        )

    for shape_config in domain_config["shape_ellipses"]:
        shape = Ellipse(
            (shape_config["center_x"], shape_config["center_y"]),
            shape_config["width"], shape_config["height"],
            shape_config["angle_in_degrees_anti-clockwise"],
        )
        domain.add_shape(
            shape,
            outline_color=shape_config["outline_color"],
            fill_color=shape_config["fill_color"],
        )

    for shape_config in domain_config["shape_rectangles"]:
        shape = Rectangle(
            (shape_config["bottom_left_x"], shape_config["bottom_left_y"]),
            shape_config["width"], shape_config["height"],
            shape_config["angle_in_degrees_anti-clockwise"],
        )
        domain.add_shape(
            shape,
            outline_color=shape_config["outline_color"],
            fill_color=shape_config["fill_color"],
        )

    for shape_config in domain_config["shape_polygons"]:
        shape = Polygon(shape_config["xy"])
        domain.add_shape(
            shape,
            outline_color=shape_config["outline_color"],
            fill_color=shape_config["fill_color"],
        )
        
    LOGGER.info(f"Building domain {domain_index}: {name} => build_domain()")    
    domain.build_domain()
    
    LOGGER.info(
        f"Building domain {domain_index}: {name} => "
        f"Add the destinations described in the JSON "
        f"configuration")

    for destination_index, destination_config in enumerate(
            domain_config["destinations"]
        ):
        desired_velocity_from_color = [
            np.concatenate((cdv["color"], cdv["desired_velocity"]))
            for cdv in destination_config["desired_velocity_from_color"]
        ]
        destination = Destination(
            name=destination_config["name"], 
            colors=destination_config["colors"],
            excluded_colors=destination_config["excluded_colors"],
            desired_velocity_from_color=desired_velocity_from_color,
            velocity_scale=destination_config["velocity_scale"],
            next_destination=destination_config["next_destination"],
            next_domain=destination_config["next_domain"],
            next_destination_weights=destination_config["next_destination_weights"],
            next_transit_box=destination_config["next_transit_box"],
            next_transit_box_xshift=destination_config["next_transit_box_xshift"],
            next_transit_box_yshift=destination_config["next_transit_box_yshift"],
        )
        LOGGER.info(
                f"Building domain {domain_index}: {name} => "
                f"Add destination {destination_index}: {destination_config['name']}")
        domain.add_destination(destination)
        if (config["with_graphes"]):
            domain.plot_desired_velocity(
                destination_config["name"], 
                id=100*domain_index+10+destination_index, 
                step=20,
                title=f"Destination '{destination_config['name']}' "
                      f"in '{name}' domain")
     
    if (config["with_graphes"]):
        domain.plot(
            id=100*domain_index, 
            title=f"Geometry of the '{name}' domain")
        domain.plot_wall_dist(
            id=100*domain_index+1, 
            step=20,
            title=f"Distance from walls for the '{name}' domain")
       
    domains[name] = domain


"""
    To create the sensors to measure the pedestrian flows
"""

all_sensors = {domain_name: [] for domain_name in domains}

for sensor_config in config["sensors"]:
    sensor = dict(sensor_config)
    sensor["id"] = []
    sensor["times"] = []
    sensor["xy"] = []
    sensor["dir"] = []
    LOGGER.info(
        f"Initializing the '{sensor['name']}' sensor "
        f"located in the '{sensor['domain']}' domain")
    all_sensors[sensor["domain"]].append(sensor)


"""
    Initialize pedestrian populations in all configured domains
"""

# Current time
t = 0.0
counter = 0
dt = config["dt"]
Tf = config["Tf"]

all_people = {}

for people_index, people_config in enumerate(config["people_init"]):
    domain_name = people_config["domain"]
    domain = domains[domain_name]
    LOGGER.info(
        f"Initializing people configuration in the '{domain_name}' domain"
    )

    people = people_initialization(
        domain,
        people_config["groups"],
        dt,
        dmin_people=config["dmin_people"],
        dmin_walls=config["dmin_walls"],
        seed=config["seed"],
        itermax=10,
        projection_method="uzawa",
        flow_duration=60,
        log=True,
    )
    if people is None:
        raise RuntimeError(
            f"No population was created for domain {domain_name}."
        )
    else:
        all_people[domain_name] = people


"""
    Visualization: either using Matplotlib or pyvista 
"""

all_visu = {} # one per domain...

for domain_index, (name, domain) in enumerate(domains.items()):

    LOGGER.info(
        f"Creating the view vindow for the '{domain_name}' domain "
        f"(render={config['plot_module']})"
    )
    
    if not config["with_graphes"]:
    
        all_visu[name] = None
    
    else:

        people = all_people[name]
        colors = people["xyrv"][:, 2]
        if colors.size>0:
            cmin, cmax = colors.min(), colors.max()
        else:
            cmin, cmax = 0, 1
                
        II, JJ, Vd = domain.people_desired_velocity(
            people["xyrv"], 
            people["destinations"]
        )
        people["Vd"] = Vd
                
        contacts = compute_contacts(domain, people["xyrv"], config["dmax"])      

        if config["plot_module"] == "matplotlib":
            
            from cromosim.visu import matplotlibVisualization
            
            visu = matplotlibVisualization(
                domain,
                figure_id=100 * (1+domain_index) + 20,
                cmap="viridis",
                velocity_reference=config["plot_velocity_reference"],
                velocity_arrow_length=config["plot_velocity_arrow_length"]
            )
            visu.update(
                people=people,
                contacts=contacts,
                sensors=all_sensors[name],
                time=t,
                plot_people=config["plot_people"],
                plot_contacts=config["plot_contacts"],
                plot_paths=config["plot_paths"],
                plot_velocities=config["plot_velocities"],
                plot_desired_velocities=config["plot_desired_velocities"],
                plot_sensors=config["plot_sensors"],
                color_scalar_values=colors,
                color_scalar_name="r",
                color_minmax=(cmin, cmax),
                color_cmap="viridis",
                show_scalar_bar=config["plot_people"],
                savefig=True,
                filename=os.path.join(
                    config["prefix"],
                    f"{name}_fig_{counter:06d}.png",
                ),
            )
            plt.pause(0.001)
            visu.show()
            all_visu[name] = visu

        elif config["plot_module"] == "pyvista":
        
            from cromosim.visu import pyvistaVisualization
        
            target_width = 1920
            target_height = round(
                target_width * domain.height / domain.width
            )    
            # scale = max(1, round(1920 / domain.width))
            movie_path = os.path.join(config["prefix"], f"{name}.mp4")
            visu = pyvistaVisualization(
                domain,
                image_path=domain._Domain__image_filename,
                pixel_size=domain.pixel_size,
                movie_path=movie_path,
                fps=4, # max(1,int(round(1/(config["plot_period"]*dt)))),
                # window_size=(scale * domain.width, scale * domain.height),
                window_size=(target_width, target_height),
                movie_macro_block_size=16,
                velocity_reference=config["plot_velocity_reference"],
                velocity_arrow_length=config["plot_velocity_arrow_length"]
            )
            visu.update(
                people=people,
                contacts=contacts,
                sensors=all_sensors[name],
                time=t,
                plot_people=config["plot_people"],
                plot_contacts=config["plot_contacts"],
                plot_velocities=config["plot_velocities"],
                plot_desired_velocities=config["plot_desired_velocities"],
                plot_paths=config["plot_paths"],
                plot_sensors=config["plot_sensors"],
                color_scalar_values=colors,
                color_scalar_name="r",
                color_minmax=(cmin, cmax),
                color_cmap="viridis",
                show_scalar_bar=config["plot_people"],
                write_frame=True,
            )
            visu.show()
            all_visu[name] = visu
                

# It may be helpful to pause here to check the calculation 
# domains and the positions of the individuals before 
# continuing the simulation...
# sys.exit()

"""
    Main loop
"""

cc = 0
draw = True

while (t < Tf):
    
    LOGGER.info( f"-------------------------------------\n" )
    LOGGER.info( f"===> Time = {t}" )

    # If there are inflows of people arriving in the area, we calculate 
    # the number of people to be added.
    
    for people_index, people_config in enumerate(config["people_init"]):
        domain_name = people_config["domain"]
        domain = domains[domain_name]
        people_to_add = add_people_flows(
            domain,
            people_config["groups"],
            dt,
            dmin_people=config["dmin_people"],
            dmin_walls=config["dmin_walls"],
            rng=all_people[domain_name]["rng"],
            itermax=10,
            projection_method="uzawa",
            flow_duration=config["dt"],
            log=True,
        )
        if people_to_add is not None: 
            LOGGER.info(
                f"People flow in the '{domain_name}' domain: "
                f"number of people added = {people_to_add['xyrv'].shape[0]}"
            )
            append_people(all_people[domain_name], people_to_add)
                            
 
     # Compute people desired velocity
     
    for domain_index, (name, domain) in enumerate(domains.items()):
        LOGGER.info(
            f"Compute people desired velocity in the '{name}' domain"
            )
        people = all_people[name]
        II, JJ, Vd = domain.people_desired_velocity(
            people["xyrv"],
            people["destinations"])
        people["Vd"][:,0] = people["beta"]*Vd[:,0]
        people["Vd"][:,1] = people["beta"]*Vd[:,1]
        people["I"] = II
        people["J"] = JJ

    # Look at if there are people in the transit boxes
    
    LOGGER.info(
                f"Find the people who need to be duplicated when "
                f"they are in an area where two domains overlap "
                f"(such as when changing floors)"
            )
    virtual_people = find_duplicate_people(all_people, domains)
 
    # Projection
    
    for domain_index, (name, domain) in enumerate(domains.items()):
        LOGGER.info(
                    f"Start projection step for the '{name}' domain "
                    f"to prevent individuals from overlapping..."
                )
        people = all_people[name]

        if virtual_people[name]["xyrv"] is not None:
            LOGGER.info(
                        f"Projection step for the '{name}' domain: "
                        f"{virtual_people[name]['xyrv'].shape[0]} "
                        f"virtual guys found"
                    )
            xyrv = np.concatenate(
                (people["xyrv"], virtual_people[name]["xyrv"]))
            Vd = np.concatenate(
                (people["Vd"], virtual_people[name]["Vd"]))
            Mass = np.concatenate(
                (people["mass"], virtual_people[name]["mass"]))
        else:
            xyrv = people["xyrv"]
            Vd = people["Vd"]
            Mass = people["mass"]

        if (xyrv.shape[0] > 0):

            if (np.unique(xyrv, axis=0).shape[0] != xyrv.shape[0]):
                raise ValueError( 
                                 "There are two identical lines in "
                                 "the array xyrv used to determine the "
                                 "contacts between individuals and this "
                                 "is not normal!"
                            )
                         
            contacts = compute_contacts(domain, xyrv, config["dmax"])
            LOGGER.info(
                f"Projection step for the '{name}' domain: "
                f"{contacts.shape[0]} contacts computed"
            )

            start = time.perf_counter()
            info, B, U, L, P, frustrations = projection(
                dt, xyrv,
                contacts, Vd, Mass,
                dmin_people=config["dmin_people"],
                dmin_walls=config["dmin_walls"],
                method=config["projection_method"],
                log=True)
            cputime_projection = time.perf_counter() - start
            LOGGER.info(
                        f"Projection step for the '{name}' domain: "
                        f"info={info} (iterations for the projection algorithm) "
                        f"cputime={np.round(cputime_projection,4)}"
                    )
              
            LOGGER.info(
                    f"Update the positions of the individuals using "
                    f"the new permissible velocities obtained"
                )
                                    
            nn = people["xyrv"].shape[0]
            all_people[name]["U"] = U[:nn, :]
            
            virtual_people[name]["U"] = U[nn:, :]
               
            # We retain the last frustration values and compute beta 
                 
            all_people[name]["frustration_last_pos"] = (all_people[name]["frustration_last_pos"] + 1) % \
                config["frustration_memory_size"]
            
            all_people[name]["frustration_memory"      ][ np.arange(nn), all_people[name]["frustration_last_pos"] ] = frustrations[:nn]
            all_people[name]["frustration_memory_times"][ np.arange(nn), all_people[name]["frustration_last_pos"] ] = t
            
            # The correcting factor, beta, is activated only 
            # for individual who have neighbors in their cone
            # of vision and at a distance less than a disc diameter.
            
            all_people[name]["beta"] = np.ones((nn,))
            
            active = np.zeros((nn,))
            for ic, cont in enumerate(contacts):
                i = int(cont[0])
                j = int(cont[1])
                dij = cont[2]
                eij_x = cont[3]
                eij_y = cont[4]
                if (j > -1):  # contact person/person
                    # theta_ij is the angle at which i (assumed to look 
                    # ahead along their desired velocity) sees j
                    # Angular dependence
                    norm_Vdi = np.sqrt(Vd[i, 0]**2+Vd[i, 1]**2)
                    if (norm_Vdi > 0) and (i<nn): # not virtual people...
                        theta_ij = np.arccos( max(-1.0, min(1.0, ( (Vd[i,0]*eij_x+Vd[i,1]*eij_y)/norm_Vdi ))))
                        if theta_ij<=all_people[name]["viewing_angle"][i]/180.0*np.pi:
                            active[i] = 1
                                                    
                    norm_Vdj = np.sqrt(Vd[j, 0]**2+Vd[j, 1]**2)
                    if (norm_Vdj > 0) and (j<nn): # not virtual people...
                        theta_ji = np.arccos( max(-1.0, min(1.0, ( -(Vd[j,0]*eij_x+Vd[j,1]*eij_y)/norm_Vdj ))))
                        if theta_ji<=all_people[name]["viewing_angle"][j]/180.0*np.pi:
                            active[j] = 1

            ind_frust = np.where( (all_people[name]["with_frustration"]==1) * \
                (active==1) )[0]
                                                                
            all_people[name]["beta"][ind_frust] = 1.0 - \
                np.mean(all_people[name]["frustration_memory"], axis=1)[ind_frust]
                                    
            # Update people positions
            
            all_people[name], all_sensors[name] = move_people(
                t, dt,
                all_people[name],
                all_sensors[name])

    # Update people destinations
    
    # the distance at which one changes destinations:
    dist_to_change_dest = 2*dt
    all_people = people_update_destination(
        all_people, domains, 
        dist_to_change_dest)

  
    t += dt
    cc += 1
    counter += 1
  
    for domain_index, (name, domain) in enumerate(domains.items()):
        LOGGER.info(
            f"--- Number of persons in the '{name}' "
            f"domain: {all_people[name]['xyrv'].shape[0]}"
        )    

    if (cc >= config["plot_period"]):
        draw = True
        cc = 0
    else:
        draw = False
    
    # Ntot = sum([people["xyrv"].shape[0] for (name,people) in all_people.items()])
    # if (Ntot==0):
    #     LOGGER.info(
    #         f"No one left in any domain! END..."
    #     )
    #     break    
    
    """
    Visualization: either using Matplotlib or pyvista 
    """

    if draw and config["with_graphes"]:

        for domain_index, (name, domain) in enumerate(domains.items()):

            people = all_people[name]
            
            colors = people["xyrv"][:, 2]
            if colors.size>0:
                cmin, cmax = colors.min(), colors.max()
            else:
                cmin, cmax = 0, 1
            
            # Since people have moved, we recalculate Vd and 
            # the contacts to ensure a correct display
            II, JJ, Vd = domain.people_desired_velocity(
                people["xyrv"], 
                people["destinations"]
            )
            people["Vd"][:,0] = people["beta"]*Vd[:,0]
            people["Vd"][:,1] = people["beta"]*Vd[:,1]        
            
            contacts = compute_contacts(domain, people["xyrv"], config["dmax"]) 
            
            if config["plot_module"] == "matplotlib":
            
                all_visu[name].update(
                    people=people,
                    contacts=contacts,
                    virtual_people=virtual_people[name],
                    sensors=all_sensors[name],
                    time=t,
                    plot_people=config["plot_people"],
                    plot_contacts=config["plot_contacts"],
                    plot_paths=config["plot_paths"],
                    plot_velocities=config["plot_velocities"],
                    plot_desired_velocities=config["plot_desired_velocities"],
                    plot_sensors=config["plot_sensors"],
                    color_scalar_values=colors,
                    color_scalar_name="r",
                    color_minmax=(cmin, cmax),
                    color_cmap="viridis",
                    show_scalar_bar=config["plot_people"],
                    savefig=True,
                    filename=os.path.join(
                        config["prefix"],
                        f"{name}_fig_{counter:06d}.png",
                    ),
                )
                plt.pause(0.001)
                all_visu[name].show()

            elif config["plot_module"] == "pyvista":
                       
                all_visu[name].update(
                    people=people,
                    contacts=contacts,
                    virtual_people=virtual_people[name],
                    sensors=all_sensors[name],
                    time=t,
                    plot_people=config["plot_people"],
                    plot_contacts=config["plot_contacts"],
                    plot_velocities=config["plot_velocities"],
                    plot_desired_velocities=config["plot_desired_velocities"],
                    plot_paths=config["plot_paths"],
                    plot_sensors=config["plot_sensors"],
                    color_scalar_values=colors,
                    color_scalar_name="r",
                    color_minmax=(cmin, cmax),
                    color_cmap="viridis",
                    show_scalar_bar=config["plot_people"],
                    write_frame=True,
                )
                all_visu[name].process_events()
                all_visu[name].plotter.update()
                all_visu[name].show()      


# End...

LOGGER.info(
    f"============ It's over! ============\n"
)    


if config["with_graphes"]:
    
    for domain_index, (name, domain) in enumerate(domains.items()):
        
        if config["plot_module"] == "pyvista":
            LOGGER.info(
                f"Close the pyvista visualization window for the '{name}' domain"
            )    
            all_visu[name].close()
        
        if config["plot_sensors"]:
            sensors = all_sensors[domain_name]
            LOGGER.info(
                f"Plot sensor data for the '{name}' domain "
                f"({len(sensors)} sensors)"
            )    
            plot_sensors(
                100*domain_index+40, 
                sensors, 
                t, 
                savefig=True,
                filename=os.path.join(
                    config["prefix"],
                    f"sensors_{domain_name}.png"
                    )
                )
            plt.pause(0.001)

plt.show()