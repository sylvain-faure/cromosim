# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
# License: GPL

import numpy as np
import scipy as sp
import sys
import matplotlib.pyplot as plt

from matplotlib.patches import Ellipse, Circle, Rectangle, Polygon
from matplotlib.lines import Line2D
from PIL import Image
from PIL import ImageDraw
import skfmm


class Destination():
    """
    A Destination object contains all the information necessary to direct
    people to a goal, for example a door, a staircase or another floor.

    Attributes
    ----------

    name: string
       name of the destination
    colors: list
       List of colors ``[ [r,g,b],... ]`` drawing the destination.
       For example, a door can be represented by a red line.
    excluded_colors: list
        List of colors ``[ [r,g,b],... ]`` representing obstacles that cannot be crossed for
        someone wishing to go to this destination: the walls of course, but
        possibly another objects only visible to the people concerned by
        this destination.
    desired_velocity_from_color: list
        Allows you to associate a desired velocity (vx, vy) with a
        color (r, g, b):
            ``[ [r,g,b, vx,vy],... ]``.
    fmm_speed: numpy array
        To automatically calculate the desired velocity leading to this
        destination, a Fast-Marching method is used to calculate the travel
        time to this destination. The opposite of the gradient of this time
        gives the desired velocity. This method solves the Eikonal equation
        ``|grad D| = 1/fmm_speed``. By changing the ``fmm_speed`` (1 everywhere by
        default) you can make certain areas slower and thus modify the
        desired velocity to divert people
    velocity_scale: float
        Multiplying coefficient used in front of the desired velocity vector
        (which is renormalized). For example on a staircase one may wish to
        reduce the speed of people.
    next_destination: string
        Name of the next destination. Useful for linking destinations one
        after the other
    next_domain: string
        Name of the next domain where is the next_destination
    next_transit_box: list
        The people in the current domain present in this box will be duplicated
        in the next_domain. A box is defined by four points:
            ``[x0, y0, x1, y1, x2, y2, x3, y3]``
    distance: numpy array
        Distance (if ``fmm_speed`` == 1) or travel time (if ``fmm_speed`` != 1)
        to the destination
    desired_velocity_X: numpy array
        First component of the desired velocity
    desired_velocity_Y: numpy array
        Second component of the desired velocity

    Examples
    --------

     * Examples with one destination:
        ``cromosim/examples/domain/domain_room.py``
        ``cromosim/examples/domain/domain_stadium.py``
     * Example with several destinations:
        ``cromosim/examples/domain/domain_shibuya_crossing.py``
     * Example with a destination described in a json file
        ``cromosim/examples/domain/domain_from_json.py``
    """
    def __init__(self,
                 name,
                 colors,
                 excluded_colors=[],
                 desired_velocity_from_color=[],
                 fmm_speed=None,
                 velocity_scale=1,
                 next_destination=None,
                 next_domain=None,
                 next_transit_box=None):
        """
        Constructor of a Destination object

        Parameters
        ----------

        name: string
           name of the destination
        colors: list ``[ [r,g,b],... ]``
           List of colors drawing the destination. For example, a door can be
           represented by a red line.
        excluded_colors: = list ``[ [r,g,b],... ]``
            List of colors representing obstacles that cannot be crossed for
            someone wishing to go to this destination: the walls of course, but
            possibly another objects only visible to the people concerned by
            this destination.
        desired_velocity_from_color: list ``[ [r,g,b, vx,vy],... ]``
            Allows you to associate a desired velocity (vx, vy) with a
            color (r, g, b)
        fmm_speed:= numpy array
            To automatically calculate the desired velocity leading to this
            destination, a Fast-Marching method is used to calculate the travel
            time to this destination. The opposite of the gradient of this time
            gives the desired velocity. This method solves the Eikonal equation
            |grad D| = 1/fmm_speed. By changing the fmm_speed (1 everywhere by
            default) you can make certain areas slower and thus modify the
            desired velocity to divert people
        velocity_scale: float
            Multiplying coefficient used in front of the desired velocity vector
            (which is renormalized). For example on a staircase one may wish to
            reduce the speed of people.
        next_destination: string
            Name of the next destination. Useful for linking destinations one
            after the other
        next_domain:
            Name of the next domain where is the next_destination
        next_transit_box:
            The people in the current domain present in this box will be duplicated
            in the ``next_domain``.

        """
        self.name = name
        self.colors = colors
        self.excluded_colors = excluded_colors
        self.desired_velocity_from_color = desired_velocity_from_color
        self.fmm_speed = fmm_speed
        self.velocity_scale = velocity_scale
        self.next_destination = next_destination
        self.next_transit_box = next_transit_box  # [xmin,xmax,ymin,ymax]
        self.next_domain = next_domain
        self.distance = None
        self.desired_velocity_X = None
        self.desired_velocity_Y = None

    def __str__(self):
        """
        Print this Destination object
        """
        return "--> Destination: " \
            + "\n    name: "+str(self.name) \
            + "\n    colors: "+str(self.colors) \
            + "\n    excluded_colors: "+str(self.excluded_colors) \
            + "\n    desired_velocity_from_color: "+str(self.desired_velocity_from_color) \
            + "\n    next_destination: "+str(self.next_destination) \
            + "\n    next_domain: "+str(self.next_domain) \
            + "\n    velocity_scale: "+str(self.velocity_scale)


class Domain():
    """
    To define the computational domain:
     * a background: empty (white) or a PNG image which only
       contains the colors white, red (for the doors) and black
       (for the walls)
     * supplementary doors represented by matplotlib shapes:
        ``line2D``
     * supplementary walls represented by matplotlib shapes:
        ``line2D``, ``circle``, ``ellipse``, ``rectangle`` or ``polygon``

    To compute the obstacle distances and the desired velocities

    Attributes
    ----------

    pixel_size: float
        size of a pixel in meters
    width: int
        width of the background image (number of pixels)
    height: int
        height of the background image (number of pixels)
    xmin: float
        x coordinate of the origin (bottom left corner)
    xmax: float
        ``xmax = xmin + width*pixel_size``
    ymin: float
        y coordinate of the origin (bottom left corner)
    ymax: float
        ``ymax = ymin + height*pixel_size``
    X: numpy array
        x coordinates (meshgrid)
    Y: numpy array
        y coordinates (meshgrid)
    image: numpy array
        pixel array (r,g,b,a)
        The Pillow image is converted to a numpy arrays, then
        using ``flipud``
        the origin of the image is put it down left instead the
        top left
    image_red: numpy array
        red values of the image (r,g,b,a)
    image_green: numpy array
        green values of the image (r,g,b,a)
    image_blue: numpy array
        blue values of the image (r,g,b,a)
    wall_mask: numpy array
        boolean array: true for wall pixels
    wall_mask_id: numpy array
        wall pixel indices
    wall_distance: numpy array
        distance (m) to the wall
    wall_grad_X: numpy array
        gradient of the distance to the wall (first component)
    wall_grad_Y: numpy array
        gradient of the distance to the wall (second component)
    door_distance: numpy array
        distance (m) to the door
    desired_velocity_X: numpy array
        opposite of the gradient of the distance to the door: desired velocity
         (first component)
    desired_velocity_Y: numpy array
        opposite of the gradient of the distance to the door: desired velocity
        (second component)

    Examples
    --------

    * A simple room
        ``cromosim/examples/domain/domain_room.py``
    * A circular domain
        ``cromosim/examples/domain/domain_stadium.py``
    * A domain with several destinations
        ``cromosim/examples/domain/domain_shibuya_crossing.py``
    * A domain built from json file (where is its description)
        ``cromosim/examples/domain/domain_from_json.py``
    """

    def __init__(self, name='Domain', background='White', pixel_size=1.0,
                 xmin=0.0, width=100,
                 ymin=0.0, height=100,
                 wall_colors=[[0, 0, 0]],
                 npzfile=None):
        """
        Constructor of a Domain object

        Parameters
        ----------

        name: string
            domain name (default: 'Domain')
        background: string
            name of the background image (default: 'White', no image)
        pixel_size: float
            size of a pixel in meters (default: 1.0)
        xmin: float
            x coordinate of the origin, bottom left corner (default: 0.0)
        ymin: float
            y coordinate of the origin, bottom left corner (default: 0.0)
        width: int
            width of the background image (default: 100 pixels)
        height: int
            height of the background image (default: 100 pixels)
        npzfile: string
            to build domain from a npz file which contains all variables
        """
        if (npzfile is None):
            self.__shapes = []
            self.__outline_color_shapes = []
            self.__fill_color_shapes = []
            self.__image_filename = ''
            self.name = name
            self.__background = background
            self.destinations = None
            self.pixel_size = pixel_size
            self.xmin, self.ymin = [xmin, ymin]
            if (self.__background != 'White'):
                self.image = Image.open(self.__background)
                self.width = self.image.size[0]
                self.height = self.image.size[1]
            else:
                self.width, self.height = [width, height]
                self.image = Image.new("RGB", (self.width, self.height), "white")

            self.xmax = self.xmin + self.width*pixel_size
            self.ymax = self.ymin + self.height*pixel_size
            self.X, self.Y = np.meshgrid(np.arange(self.width), np.arange(self.height))
            self.X = 0.5*self.pixel_size + self.xmin + self.X*self.pixel_size
            self.Y = 0.5*self.pixel_size + self.ymin + self.Y*self.pixel_size
            self.wall_colors = wall_colors  # walls are black by default
            self.wall_mask = None
            self.wall_id = None
            self.wall_distance = None
            self.wall_grad_X = None
            self.wall_grad_Y = None
            self.draw = ImageDraw.Draw(self.image)

        else:
            data = np.load(npzfile, allow_pickle=True)
            self.__shapes = data["__shapes"].tolist()
            self.__outline_color_shapes = data["__outline_color_shapes"].tolist()
            self.__fill_color_shapes = data["__fill_color_shapes"].tolist()
            self.__image_filename = data["__image_filename"]
            # print("read name: ",data["name"])
            self.name = str(data["name"])
            # print("read __background: ",data["__background"])
            self.__background = str(data["__background"])
            # print("read destinations: ",data["destinations"])
            self.destinations = dict(data["destinations"].tolist())
            # print("read pixel_size: ",data["pixel_size"])
            self.pixel_size = data["pixel_size"]
            # print("read xmin: ",data["xmin"])
            self.xmin = data["xmin"]
            # print("read ymin: ",data["ymin"])
            self.ymin = data["ymin"]
            # print("read xmax: ",data["xmax"])
            self.xmax = data["xmax"]
            # print("read ymax: ",data["ymax"])
            self.ymax = data["ymax"]
            # print("read width: ",data["width"])
            self.width = data["width"]
            # print("read height: ",data["height"])
            self.height = data["height"]
            # print("read X: ",data["X"])
            self.X = data["X"]
            # print("read Y: ",data["Y"])
            self.Y = data["Y"]
            # print("read wall_colors: ",data["wall_colors"])
            self.wall_colors = data["wall_colors"]
            # print("read wall_mask: ",data["wall_mask"])
            self.wall_mask = data["wall_mask"]
            # print("read wall_id: ",data["wall_id"])
            self.wall_id = data["wall_id"]
            # print("read wall_distance: ",data["wall_distance"])
            self.wall_distance = data["wall_distance"]
            # print("read wall_grad_X: ",data["wall_grad_X"])
            self.wall_grad_X = data["wall_grad_X"]
            # print("read wall_grad_Y: ",data["wall_grad_Y"])
            self.wall_grad_Y = data["wall_grad_Y"]
            self.image = data["image"]

    def save(self, outfile):
        """To save the content of the domain in a file

        Parameters
        ----------

        outfile: string
            output filename
        """
        np.savez(
            outfile,
            __shapes=self.__shapes,
            __outline_color_shapes=self.__outline_color_shapes,
            __fill_color_shapes=self.__fill_color_shapes,
            __image_filename=self.__image_filename,
            name=self.name,
            __background=self.__background,
            pixel_size=self.pixel_size,
            xmin=self.xmin,
            ymin=self.ymin,
            xmax=self.xmax,
            ymax=self.ymax,
            width=self.width,
            height=self.height,
            X=self.X,
            Y=self.Y,
            wall_colors=self.wall_colors,
            wall_mask=self.wall_mask,
            wall_id=self.wall_id,
            wall_distance=self.wall_distance,
            wall_grad_X=self.wall_grad_X,
            wall_grad_Y=self.wall_grad_Y,
            destinations=self.destinations,
            image=self.image)

    def add_shape(self, shape, outline_color=[0, 0, 0], fill_color=[255, 255, 255]):
        """To add a matplotlib shape:
        ``line2D``, ``circle``, ``ellipse``, ``rectangle`` or ``polygon``

        Parameters
        ----------

        shape: matplotlib shape
            line2D, circle, ellipse, rectangle or polygon
        outline_color: list
            rgb color
        fill_color: list
            rgb color
        """
        self.__shapes.append(shape)
        self.__outline_color_shapes.append(outline_color)
        self.__fill_color_shapes.append(fill_color)
        if (isinstance(shape, Circle) or isinstance(shape, Ellipse) or
                isinstance(shape, Rectangle) or isinstance(shape, Polygon)):
            xy = shape.get_verts()/self.pixel_size
            xy[:, 1] = self.height - xy[:, 1]
            self.draw.polygon(
                np.around(xy.flatten()).tolist(),
                outline="rgb("+str(outline_color[0])+"," +
                str(outline_color[1])+"," +
                str(outline_color[2])+")",
                fill="rgb("+str(fill_color[0])+"," +
                str(fill_color[1])+"," +
                str(fill_color[2])+")")
            linewidth = shape.get_linewidth()
            self.draw.line(
                np.around(xy.flatten()).tolist(),
                width=int(linewidth),
                fill="rgb("+str(outline_color[0])+"," +
                str(outline_color[1])+"," +
                str(outline_color[2])+")")
        elif isinstance(shape, Line2D):
            linewidth = shape.get_linewidth()
            xy = shape.get_xydata()/self.pixel_size
            xy[:, 1] = self.height - xy[:, 1]
            self.draw.line(
                np.around(xy.flatten()).tolist(),
                width=int(linewidth),
                fill="rgb("+str(outline_color[0])+"," +
                str(outline_color[1])+"," +
                str(outline_color[2])+")")

    def build_domain(self):
        """To build the domain: reads the background image (if supplied) \
        and initializes all the color arrrays
        """
        self.__image_filename = self.name+'_domain.png'
        self.image.save(self.__image_filename)
        # Easy way to convert a Pillow image to numpy arrays...
        # The origin of img is at the top (left) and flipud allows to put it down...
        self.image = np.flipud(np.array(self.image))
        self.image_red = self.image[:, :, 0]
        self.image_green = self.image[:, :, 1]
        self.image_blue = self.image[:, :, 2]
        self.wall_mask = np.zeros_like(self.image_red)
        for c in self.wall_colors:
            self.wall_mask += (self.image_red == c[0]) \
                * (self.image_green == c[1]) \
                * (self.image_blue == c[2])
        self.wall_id = np.where(self.wall_mask > 0)
        # Compute wall distances: wall = "wall_colors" pixels
        if (self.wall_id[0].size > 0):
            self.compute_wall_distance()
        else:
            print("WARNING: Failed to compute wall distance!")
            print("WARNING: Wall colors are ", self.wall_colors)
            print("WARNING: Check that there are pixels with these colors!")
            sys.exit()

    def compute_wall_distance(self):
        """To compute the geodesic distance to the walls in using
        a fast-marching method
        """
        phi = np.ones(self.image_red.shape)
        if (len(self.wall_id[0]) > 0):
            phi[self.wall_id] = 0
            self.wall_distance = skfmm.distance(phi, dx=self.pixel_size)
            grad = np.gradient(self.wall_distance, edge_order=2)
            grad_X = grad[1]/self.pixel_size
            grad_Y = grad[0]/self.pixel_size
            norm = np.sqrt(grad_X**2+grad_Y**2)
            norm = (norm > 0) * norm + (norm == 0) * 0.001
            self.wall_grad_X = grad_X/norm
            self.wall_grad_Y = grad_Y/norm
        else:
            self.wall_distance = 1.0e99*np.ones(self.image_red.shape)

    def add_destination(self, dest):
        """To compute the desired velocities to this destination
        and then to add this Destination object to this domain.

        Parameters
        ----------

        dest: Destination
            contains the Destination object which must be added to this domain
        """

        # Compute mask for dest.excluded_colors (which can be for example the
        # wall colors...)
        excluded_color_mask = np.zeros_like(self.image_red)
        for c in dest.excluded_colors:
            excluded_color_mask += (self.image_red == c[0]) \
                * (self.image_green == c[1]) \
                * (self.image_blue == c[2])
        excluded_color_id = np.where(excluded_color_mask > 0)

        # Compute mask for dest.colors (which correspond which correspond to
        # the destination of the people)
        dest_mask = np.zeros_like(self.image_red)
        for ic, rgb in enumerate(dest.colors):
            dest_mask = (self.image_red == rgb[0]) \
                * (self.image_green == rgb[1]) \
                * (self.image_blue == rgb[2])
        mask_id = np.where((dest_mask >= 1))

        # Define rhs of the Eikonal equation (i.e. the speed)
        if (dest.fmm_speed is not None):
            if (dest.fmm_speed.shape != self.image_red.shape):
                print("Bad speed shape ! Failed to compute the destination distance...")
                sys.exit()
        else:
            dest.fmm_speed = np.ones_like(self.image_red)

        if (mask_id[0].size > 0):
            dest_mask[mask_id] = True
            phi = np.ones(self.image_red.shape)
            phi[mask_id] = 0

            phi = np.ma.MaskedArray(phi, mask=excluded_color_mask)
            dest.distance = skfmm.travel_time(phi, dest.fmm_speed, dx=self.pixel_size)
            # dest.distance = skfmm.distance(phi, dx=self.pixel_size)
            if (excluded_color_id[0].size > 0):
                tmp_dist = dest.distance.filled(9999)
            else:
                tmp_dist = dest.distance
            grad = np.gradient(tmp_dist, edge_order=2)
        else:
            dest.distance = -np.ones_like(self.image_red)
            grad = np.gradient(np.zeros_like(self.image_red), edge_order=2)

        for rgbgrad in dest.desired_velocity_from_color:
            test = (self.image_red == int(rgbgrad[0])) \
                * (self.image_green == int(rgbgrad[1])) \
                * (self.image_blue == int(rgbgrad[2]))
            indices = np.where(test)
            grad[1][indices] = -rgbgrad[3]
            grad[0][indices] = -rgbgrad[4]

        grad_X = -grad[1]/self.pixel_size
        grad_Y = -grad[0]/self.pixel_size
        norm = np.sqrt(grad_X**2+grad_Y**2)
        norm = (norm > 0)*norm + (norm == 0)*0.001
        dest.desired_velocity_X = grad_X/norm
        dest.desired_velocity_Y = grad_Y/norm
        try:
            self.destinations[dest.name] = dest
        except:
            self.destinations = {dest.name: dest}

    def people_desired_velocity(self, xyr, people_dest, II=None, JJ=None):
        """This function determines people desired velocities from the desired \
        velocity array computed by Domain thanks to a fast-marching method.

        Parameters
        ----------
        xyr: numpy array
            people coordinates and radius: x,y,r
        people_dest: list of string
            destination for each individual
        II: numpy array (None by default)
            people index i
        JJ: numpy array (None by default)
            people index j

        Returns
        -------
        II: numpy array
            people index i
        JJ: numpy array
            people index j
        Vd: numpy array
            people desired velocity
        """
        if ((II is None) or (JJ is None)):
            II = np.floor((xyr[:, 1]-self.ymin-0.5*self.pixel_size)/self.pixel_size).astype(int)
            JJ = np.floor((xyr[:, 0]-self.xmin-0.5*self.pixel_size)/self.pixel_size).astype(int)
        Vd = np.zeros((xyr.shape[0], 2))
        for dest_name in np.unique(people_dest):
            ind = np.where(np.array(people_dest) == dest_name)[0]
            scale = self.destinations[dest_name].velocity_scale
            Vd[ind, 0] = xyr[ind, 3]*scale*self.destinations[dest_name].desired_velocity_X[II[ind], JJ[ind]]
            Vd[ind, 1] = xyr[ind, 3]*scale*self.destinations[dest_name].desired_velocity_Y[II[ind], JJ[ind]]
        return II, JJ, Vd

    def people_target_distance(self, xyr, people_dest, II=None, JJ=None):
        """This function determines distances to the current target for all people

        Parameters
        ----------
        xyr: numpy array
            people coordinates and radius: ``x,y,r``
        people_dest: list of string
            destination for each individual
        II: numpy array (None by default)
            people index ``i``
        JJ: numpy array (None by default)
            people index ``j``
        Returns
        -------
        II: numpy array
            people index i
        JJ: numpy array
            people index j
        D: numpy array
            distances to the current target
        """

        if ((II is None) or (JJ is None)):
            II = np.floor((xyr[:, 1]-self.ymin-0.5*self.pixel_size)/self.pixel_size).astype(int)
            JJ = np.floor((xyr[:, 0]-self.xmin-0.5*self.pixel_size)/self.pixel_size).astype(int)
        D = np.zeros(xyr.shape[0])
        for dest_name in np.unique(people_dest):
            ind = np.where(np.array(people_dest) == dest_name)[0]
            D[ind] = self.destinations[dest_name].distance[II[ind], JJ[ind]]-xyr[ind, 2]
        return II, JJ, D

    def people_wall_distance(self, xyr, II=None, JJ=None):
        """This function determines distances to the nearest wall for all people

        Parameters
        ----------
        xyr: numpy array
            people coordinates and radius: ``x,y,r``
        II: numpy array (None by default)
            people index ``i``
        JJ: numpy array (None by default)
            people index ``j``

        Returns
        -------
        II: numpy array
            people index ``i``
        JJ: numpy array
            people index ``j``
        D: numpy array
            distances to the nearest wall
        """
        if ((II is None) or (JJ is None)):
            II = np.floor((xyr[:, 1]-self.ymin-0.5*self.pixel_size)/self.pixel_size).astype(int)
            JJ = np.floor((xyr[:, 0]-self.xmin-0.5*self.pixel_size)/self.pixel_size).astype(int)
        D = self.wall_distance[II, JJ]-xyr[:, 2]
        return II, JJ, D

    def plot(self, id=1, title="", savefig=False, filename='fig.png', dpi=150):
        """To plot the computational domain

        Parameters
        ----------

        id: integer
            Figure id (number)
        title: string
            Figure title
        savefig: boolean
            writes the figure as a png file if true
        filename: string
            png filename used to write the figure
        dpi: integer
            number of pixel per inch for the saved figure
        """
        fig = plt.figure(id)
        ax1 = fig.add_subplot(111)
        ax1.imshow(self.image, interpolation='nearest', extent=[self.xmin, self.xmax,
                   self.ymin, self.ymax], origin='lower')
        # plt.savefig('.png',dpi=dpi)
        # ax1.axes.get_xaxis().set_visible(False)
        # ax1.axes.get_xaxis().set_visible(False)
        ax1.set_axis_off()
        ax1.set_title(title)
        # fig.add_axes(ax1)
        plt.draw()
        if (savefig):
            fig.savefig(filename, dpi=dpi, bbox_inches='tight', pad_inches=0)

    def plot_wall_dist(self,
                       step=10, scale=10, scale_units='inches', id=1, title="",
                       savefig=False, filename='fig.png', dpi=150):
        """To plot the wall distances

        Parameters
        ----------

        id: integer
            Figure id (number)
        title: string
            Figure title
        savefig: boolean
            writes the figure as a png file if true
        filename: string
            png filename used to write the figure
        dpi: integer
            number of pixel per inch for the saved figure
        """
        fig = plt.figure(id)
        ax1 = fig.add_subplot(111)
        ax1.imshow(self.image, interpolation='nearest',
                   extent=[self.xmin, self.xmax, self.ymin, self.ymax], origin='lower')
        ax1.imshow(self.wall_distance, interpolation='nearest',
                   extent=[self.xmin, self.xmax, self.ymin, self.ymax], alpha=0.7,
                   origin='lower')
        ax1.quiver(self.X[::step, ::step], self.Y[::step, ::step],
                   self.wall_grad_X[::step, ::step],
                   self.wall_grad_Y[::step, ::step],
                   scale=scale, scale_units=scale_units)
        ax1.set_axis_off()
        ax1.set_title(title)
        # plt.savefig('.png',dpi=dpi)
        plt.draw()
        if (savefig):
            fig.savefig(filename, dpi=dpi, bbox_inches='tight', pad_inches=0)

    def plot_desired_velocity(self,
                              destination_name, step=10, scale=10, scale_units='inches', id=1,
                              title="", savefig=False, filename='fig.png', dpi=150):
        """To plot the desired velocity

        Parameters
        ----------
        destination_name: string
            name of the considered destination
        step: integer
            draw an arrow every step pixels
        scale: integer
            scaling for the quiver arrows
        scale_units: string
            unit name for quiver arrow scaling
        id: integer
            Figure id (number)
        title: string
            Figure title
        savefig: boolean
            writes the figure as a png file if true
        filename: string
            png filename used to write the figure
        dpi: integer
            number of pixel per inch for the saved figure
        """
        fig = plt.figure(id)
        ax1 = fig.add_subplot(111)
        ax1.imshow(self.image, interpolation='nearest',
                   extent=[self.xmin, self.xmax, self.ymin, self.ymax], origin='lower')
        ax1.imshow(self.destinations[destination_name].distance, interpolation='nearest',
                   extent=[self.xmin, self.xmax, self.ymin, self.ymax], alpha=0.7,
                   origin='lower')
        ax1.quiver(self.X[::step, ::step], self.Y[::step, ::step],
                   self.destinations[destination_name].desired_velocity_X[::step, ::step],
                   self.destinations[destination_name].desired_velocity_Y[::step, ::step],
                   scale=scale, scale_units=scale_units)
        ax1.set_title(title)
        # plt.savefig('.png',dpi=dpi)
        ax1.set_axis_off()
        plt.draw()
        if (savefig):
            fig.savefig(filename, dpi=dpi, bbox_inches='tight', pad_inches=0)

    def __str__(self):
        """To print the main caracteristics of a Domain object
        """
        return "--> "+self.name  \
            + ":\n    dimensions: ["+str(self.xmin)+","+str(self.xmax)+"]x[" \
            + str(self.ymin)+","+str(self.ymax)+"]" \
            + "\n    width: "+str(self.width)+" height: "+str(self.height) \
            + "\n    background image: "+str(self.__background) \
            + "\n    image of the domain: "+str(self.__image_filename) \
            + "\n    wall_colors: "+str(self.wall_colors) \
            + "\n    destinations: "+str(self.destinations)
