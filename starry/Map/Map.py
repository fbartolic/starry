from ..ops import Ops, vectorize, STARRY_RECTANGULAR_PROJECTION, STARRY_ORTHOGRAPHIC_PROJECTION
from .indices import get_ylm_inds, get_ul_inds
from .utils import get_ortho_latitude_lines, get_ortho_longitude_lines
from .sht import image2map, healpix2map, array2map
import numpy as np
import theano.tensor as tt
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
radian = np.pi / 180.0
degree = 1.0 / radian


__all__ = ["Map"]



class YlmBase(object):
    """

    """

    def __init__(self, ydeg, udeg, fdeg, lazy=True, quiet=False):
        """

        """
        # Instantiate the Theano ops class
        self._lazy = lazy
        self.ops = Ops(ydeg, udeg, fdeg, lazy, quiet=quiet)
        self.cast = self.ops.cast

        # Dimensions
        self._ydeg = ydeg
        self._Ny = (ydeg + 1) ** 2 
        self._udeg = udeg
        self._Nu = udeg + 1
        self._fdeg = fdeg
        self._Nf = (fdeg + 1) ** 2
        self._deg = ydeg + udeg + fdeg
        self._N = (ydeg + udeg + fdeg + 1) ** 2

        # Initialize
        self.reset()

    @property
    def lazy(self):
        """

        """
        return self._lazy

    @property
    def ydeg(self):
        """
        
        """
        return self._ydeg
    
    @property
    def Ny(self):
        """
        
        """
        return self._Ny

    @property
    def udeg(self):
        """
        
        """
        return self._udeg
    
    @property
    def Nu(self):
        """
        
        """
        return self._Nu

    @property
    def fdeg(self):
        """
        
        """
        return self._fdeg
    
    @property
    def Nf(self):
        """
        
        """
        return self._Nf

    @property
    def deg(self):
        """
        
        """
        return self._deg
    
    @property
    def N(self):
        """
        
        """
        return self._N

    @property
    def y(self):
        """

        """
        return self._y

    @property
    def u(self):
        """

        """
        return self._u

    @property
    def inc(self):
        """

        """
        return self._inc

    @inc.setter
    def inc(self, value):
        self._inc = self.cast(value)

    @property
    def obl(self):
        """

        """
        return self._obl
    
    @obl.setter
    def obl(self, value):
        self._obl = self.cast(value)

    @property
    def axis(self):
        """

        """
        return self.ops.get_axis(self._inc * radian, self._obl * radian)

    @axis.setter
    def axis(self, axis):
        """

        """
        axis = self.cast(axis)
        inc_obl = self.ops.get_inc_obl(axis)
        self._inc = inc_obl[0] * degree
        self._obl = inc_obl[1] * degree
        
    def __getitem__(self, idx):
        """

        """
        if isinstance(idx, tuple) and len(idx) == 2:
            # User is accessing a Ylm index
            inds = get_ylm_inds(self.ydeg, idx[0], idx[1])
            return self._y[inds]
        elif isinstance(idx, (int, np.int)):
            # User is accessing a limb darkening index
            inds = get_ul_inds(self.udeg, idx)
            return self._u[inds]
        else:
            raise ValueError("Invalid map index.")

    def __setitem__(self, idx, val):
        """

        """
        if isinstance(idx, tuple) and len(idx) == 2:
            # User is accessing a Ylm index
            inds = get_ylm_inds(self.ydeg, idx[0], idx[1])
            if 0 in inds:
                raise ValueError("The Y_{0,0} coefficient cannot be set.")
            if self._lazy:
                self._y = self.ops.set_map_vector(self._y, inds, val)
            else:
                self._y[inds] = val
        elif isinstance(idx, (int, np.int, slice)):
            # User is accessing a limb darkening index
            inds = get_ul_inds(self.udeg, idx)
            if 0 in inds:
                raise ValueError("The u_0 coefficient cannot be set.")
            if self._lazy:
                self._u = self.ops.set_map_vector(self._u, inds, val)
            else:
                self._u[inds] = val
        else:
            raise ValueError("Invalid map index.")

    def reset(self):
        """

        """
        y = np.zeros(self.Ny)
        y[0] = 1.0
        self._y = self.cast(y)

        u = np.zeros(self.Nu)
        u[0] = -1.0
        self._u = self.cast(u)

        f = np.zeros(self.Nf)
        f[0] = np.pi
        self._f = self.cast(f)

        self._inc = self.cast(90.0)
        self._obl = self.cast(0.0)

    def X(self, **kwargs):
        """
        Compute and return the light curve design matrix.

        """
        # Orbital kwargs
        theta = kwargs.pop("theta", 0.0)
        xo = kwargs.pop("xo", 0.0)
        yo = kwargs.pop("yo", 0.0)
        zo = kwargs.pop("zo", 1.0)
        ro = kwargs.pop("ro", 0.0)
        theta, xo, yo, zo = vectorize(theta, xo, yo, zo)
        theta, xo, yo, zo, ro = self.cast(theta, xo, yo, zo, ro)

        # Convert angles radians
        inc = self._inc * radian
        obl = self._obl * radian
        theta *= radian

        # Compute & return
        return self.ops.X(theta, xo, yo, zo, ro, 
                          inc, obl, self._u, self._f)

    def flux(self, **kwargs):
        """
        Compute and return the light curve.
        
        """
        # Orbital kwargs
        theta = kwargs.pop("theta", 0.0)
        xo = kwargs.pop("xo", 0.0)
        yo = kwargs.pop("yo", 0.0)
        zo = kwargs.pop("zo", 1.0)
        ro = kwargs.pop("ro", 0.0)
        theta, xo, yo, zo = vectorize(theta, xo, yo, zo)
        theta, xo, yo, zo, ro = self.cast(theta, xo, yo, zo, ro)

        # Convert angles radians
        inc = self._inc * radian
        obl = self._obl * radian
        theta *= radian

        # Compute & return
        return self.ops.flux(theta, xo, yo, zo, ro, 
                             inc, obl, self._y, self._u, self._f)

    def intensity(self, **kwargs):
        """
        Compute and return the intensity of the map
        at a given ``(lat, lon)`` or ``(x, y, z)``
        point on the surface.
        
        """
        # Get the Cartesian points
        lat = kwargs.pop("lat", None)
        lon = kwargs.pop("lon", None)
        if lat is None and lon is None:
            x = kwargs.pop("x", 0.0)
            y = kwargs.pop("y", 0.0)
            z = kwargs.pop("z", 1.0)
            x, y, z = vectorize(*self.cast(x, y, z))
        else:
            lat, lon = vectorize(*self.cast(lat, lon))
            lat *= radian
            lon *= radian
            xyz = self.ops.latlon_to_xyz(self.axis, lat, lon)
            x = xyz[0]
            y = xyz[1]
            z = xyz[2]

        # Compute & return
        return self.ops.intensity(x, y, z, self._y, self._u, self._f)

    def render(self, **kwargs):
        """
        Compute and return the sky-projected intensity of 
        the map on a square Cartesian grid.
        
        """
        res = kwargs.pop("res", 300)
        projection = kwargs.pop("projection", "ortho")
        if projection.lower().startswith("rect"):
            projection = STARRY_RECTANGULAR_PROJECTION
        elif projection.lower().startswith("ortho"):
            projection = STARRY_ORTHOGRAPHIC_PROJECTION
        else:
            raise ValueError("Unknown map projection.")
        theta = kwargs.pop("theta", 0.0)
        theta = vectorize(theta)

        # Convert angles radians
        inc = self._inc * radian
        obl = self._obl * radian
        theta *= radian

        # Compute & return
        return self.ops.render(res, projection, theta, inc, obl, 
                               self._y, self._u, self._f)

    def show(self, **kwargs):
        """
        Display an image of the map, with optional animation.

        """
        # Get kwargs
        cmap = kwargs.pop("cmap", "plasma")
        projection = kwargs.get("projection", "ortho")
        if projection.lower().startswith("rect"):
            projection = STARRY_RECTANGULAR_PROJECTION
        elif projection.lower().startswith("ortho"):
            projection = STARRY_ORTHOGRAPHIC_PROJECTION
        else:
            raise ValueError("Unknown map projection.")
        grid = kwargs.pop("grid", True)
        interval = kwargs.pop("interval", 75)
        mp4 = kwargs.pop("mp4", None)

        # Get the map orientation
        if self.lazy:
            inc = self.inc.eval()
            obl = self.obl.eval()
        else:
            inc = self.inc
            obl = self.obl

        # Render the map if needed
        image = kwargs.pop("image", None)
        if image is None:

            # We need to evaluate the variables so we can
            # plot the map!
            if self.lazy:

                # Get kwargs
                res = kwargs.get("res", 300)
                theta = vectorize(kwargs.pop("theta", 0.0))

                # Evaluate the variables
                inc = self._inc.eval()
                obl = self._obl.eval()
                y = self._y.eval()
                u = self._u.eval()
                f = self._f.eval()

                # Explicitly call the compiled version of `render`
                image = self.ops._compiled_render(
                    res, projection, theta * radian, inc * radian, 
                    obl * radian, y, u, f
                )

            else:

                # Easy!
                image = self.render(**kwargs)

        if len(image.shape) == 3:
            nframes = image.shape[-1]
        else:
            nframes = 1
            image = np.reshape(image, image.shape + (1,))

        # Animation
        animated = (nframes > 1)

        if projection == STARRY_RECTANGULAR_PROJECTION:
            # Set up the plot
            fig, ax = plt.subplots(1, figsize=(7, 3.75))
            extent = (-180, 180, -90, 90)

            # Grid lines
            if grid:
                latlines = np.linspace(-90, 90, 7)[1:-1]
                lonlines = np.linspace(-180, 180, 13)
                for lat in latlines:
                    ax.axhline(lat, color="k", lw=0.5, alpha=0.5, zorder=100)
                for lon in lonlines:
                    ax.axvline(lon, color="k", lw=0.5, alpha=0.5, zorder=100)
            ax.set_xticks(lonlines)
            ax.set_yticks(latlines)
            ax.set_xlabel("Longitude [deg]")
            ax.set_ylabel("Latitude [deg]")

        else:
            # Set up the plot
            fig, ax = plt.subplots(1, figsize=(3, 3))
            ax.axis('off')
            ax.set_xlim(-1.05, 1.05)
            ax.set_ylim(-1.05, 1.05)
            extent = (-1, 1, -1, 1)

            # Grid lines
            if grid:
                x = np.linspace(-1, 1, 10000)
                y = np.sqrt(1 - x ** 2)
                ax.plot(x, y, 'k-', alpha=1, lw=1)
                ax.plot(x, -y, 'k-', alpha=1, lw=1)
                lat_lines = get_ortho_latitude_lines(inc=inc, obl=obl)
                lon_lines = get_ortho_longitude_lines(inc=inc, obl=obl)
                for x, y in lat_lines + lon_lines:
                    ax.plot(x, y, 'k-', lw=0.5, alpha=0.5, zorder=100)

        # Plot the first frame of the image
        img = ax.imshow(image[:, :, 0], origin="lower", 
                        extent=extent, cmap=cmap,
                        interpolation="none",
                        vmin=np.nanmin(image), vmax=np.nanmax(image), 
                        animated=animated)

        # Display or save the image / animation
        if animated:
            interval = kwargs.pop("interval", 75)
            mp4 = kwargs.pop("mp4", None)
            
            def updatefig(i):
                img.set_array(image[:, :, i])
                return img,

            ani = FuncAnimation(fig, updatefig, interval=interval,
                                blit=False, frames=image.shape[-1])

            # Business as usual
            if (mp4 is not None) and (mp4 != ""):
                if mp4.endswith(".mp4"):
                    mp4 = mp4[:-4]
                ani.save('%s.mp4' % mp4, writer='ffmpeg')
                plt.close()
            else:
                try:
                    if 'zmqshell' in str(type(get_ipython())):
                        plt.close()
                        display(HTML(ani.to_jshtml()))
                    else:
                        raise NameError("")
                except NameError:
                    plt.show()
                    plt.close()
        else:
            plt.show()

    def load(self, image, healpix=False, **kwargs):
        """
        Load an image, array, or ``healpix`` map. 
        
        This routine uses various routines in ``healpix`` to compute the spherical
        harmonic expansion of the input image and sets the map's :py:attr:`y`
        coefficients accordingly.

        Args:
            image: A path to an image file, a two-dimensional ``numpy`` 
                array, or a ``healpix`` map array (if ``healpix==True``).
        
        Keyword arguments:
            healpix (bool): Treat ``image`` as a ``healpix`` array? 
                Default ``False``.
            sampling_factor (int): Oversampling factor when computing the 
                ``healpix`` representation of an input image or array. 
                Default 8. Increasing this number may improve the fidelity of 
                the expanded map, but the calculation will take longer.
            sigma (float): If not ``None``, apply gaussian smoothing 
                with standard deviation ``sigma`` to smooth over 
                spurious ringing features. Smoothing is performed with 
                the ``healpix.sphtfunc.smoothalm`` method. 
                Default ``None``.
        """
        # Is this a file name?
        if type(image) is str:
            y = image2map(image, lmax=self.ydeg, **kwargs)
        # or is it an array?
        elif (type(image) is np.ndarray):
            if healpix:
                y = healpix2map(image, lmax=self.ydeg, **kwargs)
            else:
                y = array2map(image, lmax=self.ydeg, **kwargs)
        else:
            raise ValueError("Invalid `image` value.")
        
        # Ingest the coefficients
        self._y = self.cast(y)

    def rotate(self, theta, axis=None):
        """

        """
        # Get inc and obl
        if axis is None:
            inc = self._inc * radian
            obl = self._obl * radian
        else:
            axis = self.cast(axis)
            inc_obl = self.ops.get_inc_obl(axis)
            inc = inc_obl[0]
            obl = inc_obl[1]

        # Reshape theta & convert to radians
        theta = self.cast(theta)
        theta = vectorize(theta)
        theta *= radian

        # Rotate
        self._y = self.ops.rotate(self._y, theta, inc, obl)

    def align(self, source=None, dest=None):
        """
        Rotate the map to align the ``source`` point/axis with the
        ``dest`` point/axis.

        The standard way of rotating maps in ``starry`` is to
        provide the axis and angle of rotation, but this isn't always
        convenient. In some cases, it is easier to specify a source
        point/axis and a destination point/axis, and rotate the map such that the
        source aligns with the destination. This is particularly useful for
        changing map projections. For instance, to view the map pole-on,

        .. code-block:: python

            map.align(source=map.axis, dest=[0, 0, 1])

        This rotates the map axis to align with the z-axis, which points
        toward the observer.

        Args:
            source: A unit vector describing the source position. 
                This point will be rotated onto ``dest``. Default 
                is the current map axis.
            dest: A unit vector describing the destination position. 
                The ``source`` point will be rotated onto this point. Default 
                is the current map axis.

        """
        if source is None:
            source = self.axis
        if dest is None:
            dest = self.axis
        source = self.cast(source)
        dest = self.cast(dest)
        self._y = self.ops.align(self._y, source, dest)


class DopplerBase(object):
    """
    
    """

    def reset(self):
        """

        """
        super(DopplerBase, self).reset()
        self._alpha = self.cast(0.0)
        self._veq = self.cast(0.0)

    @property
    def alpha(self):
        """
        The rotational shear coefficient, a float in the range ``[0, 1]``.
        
        The parameter :math:`\\alpha` is used to model linear differential
        rotation. The angular velocity at a given latitude :math:`\\theta`
        is

        :math:`\\omega = \\omega_{eq}(1 - \\alpha \\sin^2\\theta)`

        where :math:`\\omega_{eq}` is the equatorial angular velocity of
        the object.
        """
        return self._alpha
    
    @alpha.setter
    def alpha(self, value):
        self._alpha = self.cast(value)

    @property
    def veq(self):
        """The equatorial velocity of the object in arbitrary units."""
        return self._veq
    
    @veq.setter
    def veq(self, value):
        self._veq = self.cast(value)

    def _unset_doppler_filter(self):
        f = np.zeros(self.Nf)
        f[0] = np.pi
        self._f = self.cast(f)

    def _set_doppler_filter(self):
        self._f = self.ops.compute_doppler_filter(
            self._inc * radian, self._obl * radian, self._veq, self._alpha
        )

    def rv(self, **kwargs):
        """
        Compute the net radial velocity one would measure from the object.

        The radial velocity is computed as the ratio

            :math:`\\Delta RV = \\frac{\\int Iv \\mathrm{d}A}{\\int I \\mathrm{d}A}`

        where both integrals are taken over the visible portion of the 
        projected disk. :math:`I` is the intensity field (described by the
        spherical harmonic and limb darkening coefficients) and :math:`v`
        is the radial velocity field (computed based on the equatorial velocity
        of the star, its orientation, etc.)
        """
        # Orbital kwargs
        theta = kwargs.pop("theta", 0.0)
        xo = kwargs.pop("xo", 0.0)
        yo = kwargs.pop("yo", 0.0)
        zo = kwargs.pop("zo", 1.0)
        ro = kwargs.pop("ro", 0.0)
        theta, xo, yo, zo = vectorize(theta, xo, yo, zo)
        theta, xo, yo, zo, ro = self.cast(theta, xo, yo, zo, ro)

        # Convert angles radians
        inc = self._inc * radian
        obl = self._obl * radian
        theta *= radian

        # Compute
        return self.ops.rv(
            theta, xo, yo, zo, ro, inc, obl, self._y, 
            self._u, self._veq, self._alpha
        )

    def intensity(self, **kwargs):
        # Compute the velocity-weighted intensity if `rv==True`
        rv = kwargs.pop("rv", True)
        if rv:
            self._set_doppler_filter()
        res = super(DopplerBase, self).intensity(**kwargs)
        if rv:
            self._unset_doppler_filter()
        return res

    def render(self, **kwargs):
        # Render the velocity map if `rv==True`
        rv = kwargs.pop("rv", True)
        if rv:
            self._set_doppler_filter()
        res = super(DopplerBase, self).render(**kwargs)
        if rv:
            self._unset_doppler_filter()
        return res

    def show(self, **kwargs):
        # Override the `projection` kwarg if we're
        # plotting the radial velocity.
        rv = kwargs.pop("rv", True)
        if rv:
            kwargs.pop("projection", None)
            self._set_doppler_filter()
        res = super(DopplerBase, self).show(**kwargs)
        if rv:
            self._unset_doppler_filter()
        return res


class ReflectedBase(object):
    """

    """

    def render(self, **kwargs):
        source = kwargs.pop("source", [-1, 0, 0])
        res = super(ReflectedBase, self).render(**kwargs)

        # TODO!

        '''
        Scalar b = -sz
        if (likely(abs(b) < 1.0)) {
            # Compute the terminator curve
            Scalar invsr = Scalar(1.0) / sqrt(sx * sx + sy * sy);
            Scalar cosw = sy * invsr;
            Scalar sinw = -sx * invsr;
            xrot = xv * cosw + yv * sinw;
            xrot2 = xrot.array().square();
            yrot = -xv * sinw + yv * cosw;
            yterm = b * (Ones - xrot2).cwiseSqrt();
            # Compute the illumination
            I = sqrt(Scalar(1.0) - Scalar(b * b)) * yrot - b * z;
            I = (yrot.array() > yterm.array()).select(I, 0.0);
        } else if (b < 0) {
            # Noon
            I = z;
        } else {
            # Midnight
            I = z * 0.0;
        }
        '''


def Map(ydeg=0, udeg=0, doppler=False, reflected=False, lazy=True, quiet=False):
    """

    """

    Bases = (YlmBase,)

    # Doppler mode?
    if doppler:
        Bases = (DopplerBase,) + Bases
        fdeg = 3
    else:
        fdeg = 0
    
    # Reflected light?
    if reflected:
        Bases = (ReflectedBase,) + Bases

    class Map(*Bases): 
        pass

    return Map(ydeg, udeg, fdeg, lazy=lazy, quiet=quiet)