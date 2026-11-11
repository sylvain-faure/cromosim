# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
# License: GPL

"""Geometry, navigation fields, and plotting utilities for Cromosim domains."""

import logging
import sys

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Ellipse, Polygon, Rectangle
from PIL import Image
from PIL import ImageDraw
import skfmm


LOGGER = logging.getLogger(__name__)


def _rgb_mask(red, green, blue, colors):
    """Return a mask selecting pixels that match any RGB color in ``colors``."""
    mask = np.zeros_like(red)
    for color in colors:
        mask += (
            (red == color[0])
            * (green == color[1])
            * (blue == color[2])
        )
    return mask


def _grid_indices(xyr, xmin, ymin, pixel_size):
    """Return row and column indices for pedestrian coordinates."""
    rows = np.floor(
        (xyr[:, 1] - ymin - 0.5 * pixel_size) / pixel_size
    ).astype(int)
    cols = np.floor(
        (xyr[:, 0] - xmin - 0.5 * pixel_size) / pixel_size
    ).astype(int)
    return rows, cols


class Destination:
    """
    Represent a target used to guide pedestrians through a domain.

    A destination may correspond to a door, a staircase, another floor, or any
    other target encoded in the domain image. Desired velocities are obtained
    from a fast-marching travel-time field and may optionally be overridden on
    pixels of specified colors.

    Parameters and attributes
    -------------------------
    name : str
        Destination name.
    colors : list
        RGB colors ``[[r, g, b], ...]`` used to identify destination pixels.
    excluded_colors : list
        RGB colors representing regions that cannot be crossed on the way to
        this destination.
    desired_velocity_from_color : list
        Optional explicit desired velocities associated with colors, encoded as
        ``[r, g, b, vx, vy]``.
    fmm_speed : numpy.ndarray or None
        Speed field used by the fast-marching solver. If ``None``, a unit speed
        field is created when the destination is added to a domain.
    velocity_scale : float
        Multiplicative factor applied to the pedestrian desired speed.
    next_destination : object, optional
        Description of the destination(s) reached after this one.
    next_destination_choice : object, optional
        Optional data used to select the next destination.
    next_destination_weights : object, optional
        Weights associated with the possible next destinations.
    next_domain : object, optional
        Domain or domains containing the next destination(s).
    next_transit_box : object, optional
        Transit-region description used when pedestrians must be duplicated or
        transferred between domains.
    next_transit_box_xshift, next_transit_box_yshift : object, optional
        Coordinate shifts associated with transit regions.
    excluded_id : numpy.ndarray or None
        Mask identifying pixels excluded for this destination.
    distance : numpy.ndarray or None
        Fast-marching travel-time (or distance) field.
    desired_velocity_X, desired_velocity_Y : numpy.ndarray or None
        Components of the normalized desired-velocity field.

    Examples
    --------
    See ``cromosim/examples/domain/domain_room.py``,
    ``cromosim/examples/domain/domain_stadium.py``,
    ``cromosim/examples/domain/domain_shibuya_crossing.py`` and
    ``cromosim/examples/domain/domain_from_json.py``.
    """
    def __init__(
        self,
        name,
        colors,
        excluded_colors=[],
        desired_velocity_from_color=[],
        fmm_speed=None,
        velocity_scale=1,
        next_destination_choice=None,
        next_destination=None,
        next_domain=None,
        next_destination_weights=None,
        next_transit_box_xshift=None,
        next_transit_box_yshift=None,
        next_transit_box=None,
    ):
        """
        Initialize a destination.

        Parameters
        ----------
        name : str
            Destination name.
        colors : list
            RGB colors ``[[r, g, b], ...]`` identifying destination pixels.
        excluded_colors : list, optional
            RGB colors identifying regions that cannot be crossed.
        desired_velocity_from_color : list, optional
            Explicit velocity overrides encoded as ``[r, g, b, vx, vy]``.
        fmm_speed : numpy.ndarray or None, optional
            Fast-marching speed field.
        velocity_scale : float, optional
            Scale factor applied to desired pedestrian speeds.
        next_destination_choice : object, optional
            Optional next-destination selection data.
        next_destination : object, optional
            Next destination or destinations.
        next_domain : object, optional
            Domain or domains associated with the next destination(s).
        next_destination_weights : object, optional
            Selection weights associated with the next destination(s).
        next_transit_box_xshift, next_transit_box_yshift : object, optional
            Coordinate shifts associated with transit regions.
        next_transit_box : object, optional
            Transit-region description.
        """
        self.name = name
        self.colors = colors
        self.excluded_colors = excluded_colors
        self.desired_velocity_from_color = desired_velocity_from_color
        self.fmm_speed = fmm_speed
        self.velocity_scale = velocity_scale
        self.next_destination = next_destination
        self.next_destination_choice = next_destination_choice
        self.next_destination_weights = next_destination_weights
        self.next_transit_box = next_transit_box  # [xmin,xmax,ymin,ymax]
        self.next_transit_box_xshift = next_transit_box_xshift
        self.next_transit_box_yshift = next_transit_box_yshift
        self.next_domain = next_domain
        self.excluded_id = None
        self.distance = None
        self.desired_velocity_X = None
        self.desired_velocity_Y = None

    def __str__(self):
        """Return a human-readable summary of the destination."""
        return "--> Destination: " \
            + "\n    name: "+str(self.name) \
            + "\n    colors: "+str(self.colors) \
            + "\n    excluded_colors: "+str(self.excluded_colors) \
            + "\n    desired_velocity_from_color: "+str(self.desired_velocity_from_color) \
            + "\n    next_destination: "+str(self.next_destination) \
            + "\n    next_destination_choice: "+str(self.next_destination_choice) \
            + "\n    next_destination_weights: "+str(self.next_destination_weights) \
            + "\n    next_domain: "+str(self.next_domain) \
            + "\n    velocity_scale: "+str(self.velocity_scale) \
            + "\n    transit_box: "+str(self.next_transit_box) \
            + "\n    transit_box shift: ("+str(self.next_transit_box_xshift)+", "+str(self.next_transit_box_yshift)+")"


class Domain:
    """
    Represent a two-dimensional computational domain.

    A domain is defined by a raster image (or a white image created in memory),
    its physical pixel size, wall colors, and a set of destinations. Matplotlib
    shapes may be drawn onto the raster before the domain is built.

    The class computes wall-distance fields and destination travel-time fields
    with a fast-marching method. Their gradients are then used to obtain wall
    normals and pedestrian desired-velocity directions.

    Attributes
    ----------
    name : str
        Domain name.
    pixel_size : float
        Physical size of one image pixel.
    width, height : int
        Image dimensions in pixels.
    xmin, ymin : float
        Physical coordinates of the lower-left corner.
    xmax, ymax : float
        Physical coordinates of the upper-right corner.
    X, Y : numpy.ndarray
        Mesh-grid arrays containing pixel-center coordinates.
    image : PIL.Image.Image or numpy.ndarray
        Domain image. After :meth:`build_domain`, it is stored as a NumPy array
        with the vertical axis flipped so that the origin is at the lower left.
    image_red, image_green, image_blue : numpy.ndarray
        RGB channels created by :meth:`build_domain`.
    wall_colors : list
        RGB colors used to identify wall pixels.
    wall_mask : numpy.ndarray or None
        Mask of wall pixels.
    wall_id : tuple or None
        Indices of wall pixels.
    wall_distance : numpy.ndarray or None
        Distance to the nearest wall.
    wall_grad_X, wall_grad_Y : numpy.ndarray or None
        Components of the normalized wall-distance gradient.
    destinations : dict or None
        Destinations registered in the domain.

    Examples
    --------
    See ``cromosim/examples/domain/domain_room.py``,
    ``cromosim/examples/domain/domain_stadium.py``,
    ``cromosim/examples/domain/domain_shibuya_crossing.py`` and
    ``cromosim/examples/domain/domain_from_json.py``.
    """

    def __init__(
        self,
        name="Domain",
        background="White",
        pixel_size=1.0,
        xmin=0.0,
        width=100,
        ymin=0.0,
        height=100,
        wall_colors=[[0, 0, 0]],
        npzfile=None,
    ):
        """
        Initialize a computational domain.

        Parameters
        ----------
        name : str, optional
            Domain name.
        background : str, optional
            Background image filename. ``"White"`` creates a blank white image.
        pixel_size : float, optional
            Physical size of one pixel.
        xmin, ymin : float, optional
            Coordinates of the lower-left corner.
        width, height : int, optional
            Image dimensions used when ``background="White"``.
        wall_colors : list, optional
            RGB colors used to identify walls.
        npzfile : str or None, optional
            Existing ``.npz`` domain file to load instead of creating a new domain.
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

            self.people_counter = 0
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
        """
        Save the domain state to a NumPy ``.npz`` archive.

        Parameters
        ----------
        outfile : str or path-like
            Output filename. NumPy appends ``.npz`` when appropriate.
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

    def add_shape(
        self,
        shape,
        outline_color=[0, 0, 0],
        fill_color=[255, 255, 255],
    ):
        """
        Draw a Matplotlib shape onto the domain image.

        Supported shapes are :class:`~matplotlib.lines.Line2D`,
        :class:`~matplotlib.patches.Circle`, :class:`~matplotlib.patches.Ellipse`,
        :class:`~matplotlib.patches.Rectangle`, and
        :class:`~matplotlib.patches.Polygon`.

        Parameters
        ----------
        shape : matplotlib artist
            Shape to draw.
        outline_color : sequence of int, optional
            RGB outline color.
        fill_color : sequence of int, optional
            RGB fill color. It is ignored for ``Line2D`` objects.
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
        """
        Finalize the raster representation and compute wall geometry.

        The current image is saved to ``"<domain name>_domain.png"``, converted to a
        NumPy array with a lower-left origin, and split into RGB channels. Wall
        pixels are then identified and the wall-distance field is computed.
        """
        self.__image_filename = self.name+'_domain.png'
        self.image.save(self.__image_filename)
        # Easy way to convert a Pillow image to numpy arrays...
        # The origin of img is at the top (left) and flipud allows to put it down...
        self.image = np.flipud(np.array(self.image))
        self.image_red = self.image[:, :, 0]
        self.image_green = self.image[:, :, 1]
        self.image_blue = self.image[:, :, 2]
        self.wall_mask = _rgb_mask(
            self.image_red,
            self.image_green,
            self.image_blue,
            self.wall_colors,
        )
        self.wall_id = np.where(self.wall_mask > 0)
        # Compute wall distances: wall = "wall_colors" pixels
        if self.wall_id[0].size > 0:
            self.compute_wall_distance()
        else:
            LOGGER.error(
                "Failed to compute wall distance: no pixels match wall colors %s.",
                self.wall_colors,
            )
            LOGGER.error("Check the domain image and configured wall colors.")
            sys.exit()

    def compute_wall_distance(self):
        """
        Compute the distance to walls and its normalized gradient.

        The distance field is evaluated with a fast-marching method using the wall
        pixels as the zero level set.
        """
        phi = np.ones(self.image_red.shape)
        if len(self.wall_id[0]) > 0:
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
        """
        Compute and register the fields associated with a destination.

        A fast-marching travel-time field is computed outside the destination's
        excluded colors. The negative normalized gradient defines the desired
        velocity direction. Explicit color-based velocity overrides are applied
        before normalization.

        Parameters
        ----------
        dest : Destination
            Destination to add to the domain.
        """

        # Compute mask for dest.excluded_colors (which can be for example the
        # wall colors...)
        excluded_color_mask = _rgb_mask(
            self.image_red,
            self.image_green,
            self.image_blue,
            dest.excluded_colors,
        )
        excluded_color_id = np.where(excluded_color_mask > 0)
        dest.excluded_id = np.zeros_like(self.image_red)
        dest.excluded_id[ excluded_color_id ] = 1
        
        # Compute mask for dest.colors (which correspond which correspond to
        # the destination of the people)
        dest_mask = _rgb_mask(
            self.image_red,
            self.image_green,
            self.image_blue,
            dest.colors,
        )
        mask_id = np.where((dest_mask >= 1))

        # Define rhs of the Eikonal equation (i.e. the speed)
        if dest.fmm_speed is not None:
            if dest.fmm_speed.shape != self.image_red.shape:
                LOGGER.error(
                    "Cannot compute destination %r: fmm_speed has shape %s, expected %s.",
                    dest.name,
                    dest.fmm_speed.shape,
                    self.image_red.shape,
                )
                sys.exit()
        else:
            dest.fmm_speed = np.ones_like(self.image_red)

        if mask_id[0].size > 0:
            dest_mask[mask_id] = True
            phi = np.ones(self.image_red.shape)
            phi[mask_id] = 0

            phi = np.ma.MaskedArray(phi, mask=excluded_color_mask)
            dest.distance = skfmm.travel_time(phi, dest.fmm_speed, dx=self.pixel_size)
            # dest.distance = skfmm.distance(phi, dx=self.pixel_size)
            if excluded_color_id[0].size > 0:
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
        if self.destinations is None:
            self.destinations = {}
        self.destinations[dest.name] = dest




    def people_desired_velocity(self, xyr, people_dest, II=None, JJ=None):
        """
        Evaluate desired velocities for a collection of pedestrians.

        Parameters
        ----------
        xyr : numpy.ndarray
            Array with at least four columns ``[x, y, radius, speed]``.
        people_dest : array-like of str
            Destination name associated with each pedestrian.
        II, JJ : numpy.ndarray or None, optional
            Precomputed row and column indices in the domain grid.

        Returns
        -------
        II, JJ : numpy.ndarray
            Row and column indices used for interpolation by pixel lookup.
        Vd : numpy.ndarray
            Desired velocities, with shape ``(n_people, 2)``.
        """
        if II is None or JJ is None:
            II, JJ = _grid_indices(
                xyr, self.xmin, self.ymin, self.pixel_size
            )
        Vd = np.zeros((xyr.shape[0], 2))
        for dest_name in np.unique(people_dest):
            ind = np.where(np.array(people_dest) == dest_name)[0]
            scale = self.destinations[dest_name].velocity_scale
            Vd[ind, 0] = xyr[ind, 3]*scale*self.destinations[dest_name].desired_velocity_X[II[ind], JJ[ind]]
            Vd[ind, 1] = xyr[ind, 3]*scale*self.destinations[dest_name].desired_velocity_Y[II[ind], JJ[ind]]
        return II, JJ, Vd

    def people_target_distance(self, xyr, people_dest, II=None, JJ=None):
        """
        Evaluate the distance to each pedestrian's current destination.

        Parameters
        ----------
        xyr : numpy.ndarray
            Array containing at least ``[x, y, radius]`` for each pedestrian.
        people_dest : array-like of str
            Destination name associated with each pedestrian.
        II, JJ : numpy.ndarray or None, optional
            Precomputed row and column indices in the domain grid.

        Returns
        -------
        II, JJ : numpy.ndarray
            Row and column grid indices.
        D : numpy.ndarray
            Destination distance minus pedestrian radius.
        """

        if II is None or JJ is None:
            II, JJ = _grid_indices(
                xyr, self.xmin, self.ymin, self.pixel_size
            )
        D = np.zeros(xyr.shape[0])
        for dest_name in np.unique(people_dest):
            ind = np.where(np.array(people_dest) == dest_name)[0]
            D[ind] = self.destinations[dest_name].distance[II[ind], JJ[ind]]-xyr[ind, 2]
        return II, JJ, D

    def people_wall_distance(self, xyr, II=None, JJ=None):
        """
        Evaluate the distance from each pedestrian boundary to the nearest wall.

        Parameters
        ----------
        xyr : numpy.ndarray
            Array containing at least ``[x, y, radius]`` for each pedestrian.
        II, JJ : numpy.ndarray or None, optional
            Precomputed row and column indices in the domain grid.

        Returns
        -------
        II, JJ : numpy.ndarray
            Row and column grid indices.
        D : numpy.ndarray
            Wall distance minus pedestrian radius.
        """
        if II is None or JJ is None:
            II, JJ = _grid_indices(
                xyr, self.xmin, self.ymin, self.pixel_size
            )
        D = self.wall_distance[II, JJ]-xyr[:, 2]
        return II, JJ, D

    def plot(self, id=1, title="", savefig=False, filename='fig.png', dpi=150):
        """
        Plot the computational domain.

        Parameters
        ----------
        id : int, optional
            Matplotlib figure number.
        title : str, optional
            Figure title.
        savefig : bool, optional
            Save the figure when ``True``.
        filename : str, optional
            Output image filename.
        dpi : int, optional
            Resolution of the saved figure.
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
        if savefig:
            fig.savefig(filename, dpi=dpi, bbox_inches='tight', pad_inches=0)

    def plot_wall_dist(self,
                       step=10, scale=10, scale_units='inches', id=1, title="",
                       savefig=False, filename='fig.png', dpi=150):
        """
        Plot the wall-distance field and its gradient.

        Parameters
        ----------
        step : int, optional
            Plot one gradient arrow every ``step`` pixels.
        scale : float, optional
            Matplotlib quiver scale.
        scale_units : str, optional
            Matplotlib quiver scale units.
        id : int, optional
            Matplotlib figure number.
        title : str, optional
            Figure title.
        savefig : bool, optional
            Save the figure when ``True``.
        filename : str, optional
            Output image filename.
        dpi : int, optional
            Resolution of the saved figure.
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
        if savefig:
            fig.savefig(filename, dpi=dpi, bbox_inches='tight', pad_inches=0)

    def plot_desired_velocity(self,
                              destination_name, step=10, scale=10, scale_units='inches', id=1,
                              title="", savefig=False, filename='fig.png', dpi=150):
        """
        Plot a destination travel-time field and its desired-velocity directions.

        Parameters
        ----------
        destination_name : str
            Name of the destination to display.
        step : int, optional
            Plot one velocity arrow every ``step`` pixels.
        scale : float, optional
            Matplotlib quiver scale.
        scale_units : str, optional
            Matplotlib quiver scale units.
        id : int, optional
            Matplotlib figure number.
        title : str, optional
            Figure title.
        savefig : bool, optional
            Save the figure when ``True``.
        filename : str, optional
            Output image filename.
        dpi : int, optional
            Resolution of the saved figure.
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
        if savefig:
            fig.savefig(filename, dpi=dpi, bbox_inches='tight', pad_inches=0)

    def __str__(self):
        """Return a human-readable summary of the domain."""
        return "--> "+self.name  \
            + ":\n    dimensions: ["+str(self.xmin)+","+str(self.xmax)+"]x[" \
            + str(self.ymin)+","+str(self.ymax)+"]" \
            + "\n    width: "+str(self.width)+" height: "+str(self.height) \
            + "\n    background image: "+str(self.__background) \
            + "\n    image of the domain: "+str(self.__image_filename) \
            + "\n    wall_colors: "+str(self.wall_colors) \
            + "\n    destinations: "+str(self.destinations)
