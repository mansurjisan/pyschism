import numpy as np
import pathlib
import logging
from collections import defaultdict
import fiona
from functools import lru_cache
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely.geometry import LineString, mapping
import matplotlib.pyplot as plt
from pyschism.mesh import grd
from pyschism import figures as fig
from pyschism.mesh.base import EuclideanMesh2D



class Gmesh(EuclideanMesh2D):
    """
    boundaries = {ibtype: {id: {'indexes': [i0, ..., in], 'properties': object }}
    """
    def __init__(
        self,
        coords,
        triangles=None,
        quads=None,
        values=None,
        crs=None,
        description=None,
        boundaries=None,
    ):
        super().__init__(coords, triangles, quads, values, crs, description)
        self._boundaries = boundaries

    @staticmethod
    def open_grd(path, crs=None):
        return Gmesh(**grd.to_gmesh(grd.reader(path)), crs=crs)

    def add_boundary_type(self, ibtype):
        if ibtype not in self.boundaries:
            self.__boundaries[ibtype] = defaultdict()
        else:
            msg = f"Cannot add boundary_type={ibtype}: boundary type already "
            msg += "exists."
            raise Exception(msg)

    def set_boundary_data(self, ibtype, id, indexes, **properties):
        msg = "Indexes must be subset of node id's."
        for idx in np.asarray(indexes).flatten():
            assert idx in self.node_index.keys(), msg
        self.__boundaries[ibtype][id] = {
            'indexes': indexes,
            **properties
        }

    def clear_boundaries(self):
        self.__boundaries = {}
        self.__boundaries[None] = {}

    def delete_boundary_type(self, ibtype):
        del self.__boundaries[ibtype]

    def delete_boundary_data(self, ibtype, id):
        del self.__boundaries[ibtype][id]

    def write_boundaries(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.exists() and not overwrite:
            msg = "Destination path exists and overwrite=False"
            raise IOError(msg)
        with fiona.open(
                    path.absolute(),
                    'w',
                    driver='ESRI Shapefile',
                    crs=self.crs.srs,
                    schema={
                        'geometry': 'LineString',
                        'properties': {
                            'id': 'int',
                            'ibtype': 'str',
                            'bnd_id': 'str'
                            }
                        }) as dst:
            _cnt = 0
            for ibtype, bnds in self.boundaries.items():
                for id, bnd in bnds.items():
                    idxs = list(map(self.get_node_index, bnd['indexes']))
                    linear_ring = LineString(self.xy[idxs].tolist())
                    dst.write({
                            "geometry": mapping(linear_ring),
                            "properties": {
                                "id": _cnt,
                                "ibtype": ibtype,
                                "bnd_id": f"{ibtype}:{id}"
                                }
                            })
                    _cnt += 1

    @fig._figure
    def make_plot(
        self,
        axes=None,
        vmin=None,
        vmax=None,
        show=False,
        title=None,
        # figsize=rcParams["figure.figsize"],
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
        tags=True,
        axes=None,
        show=False,
        figsize=None,
        **kwargs
    ):

        boundary = list(map(
            self.get_node_index, self.boundaries[ibtype][id]['indexes']))
        p = axes.plot(self.x[boundary], self.y[boundary], **kwargs)
        if tags:
            axes.text(
                self.x[boundary[len(boundary)//2]],
                self.y[boundary[len(boundary)//2]],
                f"ibtype={ibtype}\nid={id}",
                color=p[-1].get_color()
                )
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
            for id in bnds:
                axes = self.plot_boundary(ibtype, id, **kwargs)
                kwargs.update({'axes': axes})
        return kwargs['axes']

    @property
    @lru_cache
    def logger(self):
        return logging.getLogger(__name__ + '.' + self.__class__.__name__)

    @property
    def boundaries(self):
        return self._boundaries

    @property
    def _grd(self):
        grd = super()._grd
        grd.update({'boundaries': self.boundaries})
        return grd

    @property
    def _sms2dm(self):
        sms2dm = super()._sms2dm
        def nodestrings(boundaries):
            return
        sms2dm.update({'nodestrings': nodestrings(self.boundaries)})
        return sms2dm
    
    @property
    def _boundaries(self):
        return self.__boundaries

    @_boundaries.setter
    def _boundaries(self, boundaries):
        self.clear_boundaries() # clear
        if boundaries is not None:
            for ibtype, bnds in boundaries.items():
                if ibtype is not None:
                    self.add_boundary_type(ibtype)
                for id, bnd in bnds.items():
                    if 'properties' in bnd.keys():
                        properties = bnd['properties']
                    else:
                        properties = {}
                    self.set_boundary_data(
                        ibtype,
                        id,
                        bnd['indexes'],
                        **properties
                    )
