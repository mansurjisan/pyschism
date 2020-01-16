import numpy as np
from collections import defaultdict
from collections.abc import Iterable, Mapping
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from functools import lru_cache
from pyschism.mesh import gr3
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyschism.mesh.gr3 import Gr3
import pyschism.mesh.figures as fig
from pyschism.mesh.friction import (
    Fgrid,
    ManningsN,
    DragCoefficient,
    RoughnessLength
)


class Hgrid(Gr3):
    """
    Class that represents the unstructured planar mesh used by SCHISM.
    """

    def __init__(
        self,
        nodes,
        elements,
        boundaries=None,
        crs=None,
        description=None,
    ):
        nodes = {id: (coord, -value) for id, (coord, value) in nodes.items()}
        super().__init__(nodes, elements, crs, description)
        self._boundaries = boundaries

    def set_friction(self, value, ftype='manning'):

        # certify ftype
        ftypes = {
            'manning': ManningsN,
            'drag': DragCoefficient,
            'rough': RoughnessLength
        }
        msg = f"ftype argument must be one of {ftypes.keys()}"
        assert ftype.lower() in ftypes, msg

        # certify value
        msg = "value argument must be an instance of type "
        msg += f"{int}, {float} or an iterable ."
        assert isinstance(value, (Iterable, int, float, Fgrid)), msg

        if isinstance(value, Fgrid):
            self._fgrid = value

        elif isinstance(value, (int, float)):
            if ftype == 'manning':
                self._fgrid = ftypes[ftype].constant(self, value)

        return self.fgrid

    def add_boundary_type(self, ibtype):
        if ibtype not in self.boundaries:
            self._boundaries[ibtype] = defaultdict()

    def set_boundary_data(self, ibtype, id, indexes, **properties):
        keys = set(self.nodes.keys())
        for idx in indexes:
            msg = "Indexes must be subset of node id's."
            if isinstance(idx, (str, int)):
                assert set(idx).issubset(keys), msg
            elif isinstance(idx, Iterable):
                for _idx in idx:
                    assert set(_idx).issubset(keys), msg
        self._boundaries[ibtype] = {
            'indexes': indexes,
            **properties
        }

    @fig._figure
    def make_plot(
        self,
        axes=None,
        vmin=None,
        vmax=None,
        show=False,
        title=None,
        figsize=rcParams["figure.figsize"],
        extent=None,
        cbar_label=None,
        **kwargs
    ):
        if vmin is None:
            vmin = np.min(self.values)
        if vmax is None:
            vmax = np.max(self.values)
        kwargs.update(**fig.get_topobathy_kwargs(self.values, vmin, vmax))
        kwargs.pop('col_val')
        levels = kwargs.pop('levels')
        if vmin != vmax:
            self.tricontourf(
                axes=axes,
                levels=levels,
                vmin=vmin,
                vmax=vmax,
                **kwargs
            )
        else:
            self.tripcolor(axes=axes, **kwargs)
        self.quadface(axes=axes, **kwargs)
        axes.axis('scaled')
        if extent is not None:
            axes.axis(extent)
        if title is not None:
            axes.set_title(title)
        mappable = ScalarMappable(cmap=kwargs['cmap'])
        mappable.set_array([])
        mappable.set_clim(vmin, vmax)
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("bottom", size="2%", pad=0.5)
        cbar = plt.colorbar(
            mappable,
            cax=cax,
            orientation='horizontal'
        )
        cbar.set_ticks([vmin, vmax])
        cbar.set_ticklabels([np.around(vmin, 2), np.around(vmax, 2)])
        if cbar_label is not None:
            cbar.set_label(cbar_label)
        return axes

    @fig._figure
    def plot_boundary(
        self,
        ibtype,
        id,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        boundary = [int(idx)-1 for idx in self.boundaries[ibtype][id]]
        axes.plot(self.x[boundary], self.y[boundary], **kwargs)
        return axes

    @fig._figure
    def plot_boundaries(
        self,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):
        kwargs.update({'axes': axes})
        for ibtype, bnds in self.boundaries.items():
            for bnd in bnds:
                axes = self.plot_boundary(ibtype, bnd, **kwargs)
                kwargs.update({'axes': axes})
        return kwargs['axes']

    @property
    @lru_cache
    def grd(self):
        grd = super().grd
        grd.update({"boundaries": self.boundaries})
        return grd

    @property
    def boundaries(self):
        return self._boundaries

    @property
    def fgrid(self):
        try:
            return self.__fgrid
        except AttributeError:
            self._fgrid = ManningsN.constant(self, 0.025)
            return self.__fgrid

    @property
    def _boundaries(self):
        return self.__boundaries

    @property
    def _fgrid(self):
        return self.__fgrid

    @_boundaries.setter
    def _boundaries(self, boundaries):
        """
        elements in boundaries should be a subset of the node keys.
        """
        keys = set(self.nodes)
        if boundaries is not None:
            msg = "elements argument must be a dictionary of the form "
            msg += "\\{element_id:  (e0, ..., en)\\} where n==2 or n==3."
            assert isinstance(boundaries, Mapping), msg
            msg = "(e0, ..., en) must be a subset of the node keys."
            for geom in boundaries.values():
                for bnd in geom.values():
                    assert set(bnd).issubset(keys), msg
        else:
            boundaries = {}
        # ocean boundaries
        if None not in boundaries:
            boundaries[None] = {}
        # land boundaries
        if 0 not in boundaries:
            boundaries[0] = {}
        # interior boundaries
        if 1 not in boundaries:
            boundaries[1] = {}

        self.__boundaries = boundaries

    @_fgrid.setter
    def _fgrid(self, fgrid):
        assert isinstance(fgrid, Fgrid)
        self.__fgrid = fgrid
