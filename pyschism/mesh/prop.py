import pathlib
from typing import Union

import numpy as np
from shapely.geometry import Polygon, MultiPolygon

from pyschism.mesh.base import Gr3


class Prop:

    def __init__(self, gr3: Gr3, element_values: np.array):

        if not isinstance(gr3, Gr3):
            raise TypeError(
                f'Argument gr3 must be an instance of type {Gr3}, not type '
                f'{type(gr3)}.')

        values = np.array(element_values).flatten()
        if len(gr3.elements) != values.shape[0]:
            raise ValueError(
                'Shape mismatch between element_values and hgrid.')

        self.elements = gr3.elements
        self.values = values

    def __str__(self):
        f = []
        for i, (iele, element) in enumerate(self.elements.elements.items()):
            f.append(f'{iele} {self.values[i]:G}')
        return '\n'.join(f)

    @classmethod
    def constant(cls, gr3: Gr3, value: np.array):
        return cls(gr3, np.full((len(gr3.elements),), value))

    def write(self, path, overwrite: bool = False):
        path = pathlib.Path(path)
        if path.exists() and not overwrite:
            raise IOError('path exists and overwrite is False')
        with open(path, 'w') as f:
            f.write(str(self))

    @classmethod
    def from_geometry(
            cls,
            gr3: Gr3,
            region: Union[Polygon, MultiPolygon],
            inner_value: float,
            outer_value: float,
    ):

        if not isinstance(gr3, Gr3):
            raise TypeError(
                f'Argument gr3 must be an instance of type {Gr3}, not type '
                f'{type(gr3)}.')

        if not isinstance(region, (Polygon, MultiPolygon)):
            raise TypeError(
                f'Argument region must be an instance of types {Polygon} or '
                f'{MultiPolygon}, not type {type(region)}.')

        if isinstance(region, Polygon):
            region = [region]

        obj = cls.constant(gr3, np.nan)
        elements_gdf = gr3.elements.gdf
        gdf_in = elements_gdf.geometry.intersects(region)
        obj.values[gdf_in] = inner_value
        outer_indexes = np.setdiff1d(elements_gdf.index, np.where(gdf_in))
        obj.values[outer_indexes] = outer_value
        return obj


class Fluxflag(Prop):
    """ Class for writing fluxflag.prop file, which is parameter for
        checking volume and salt conservation.
    """
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, -1)

    @classmethod
    def from_geometry(cls, gr3: Gr3, region: Union[Polygon, MultiPolygon],
                      value: int):
        if value not in [1, -1]:
            raise ValueError('Argument value must be 1 or -1.')
        return super(Fluxflag, cls).from_geometry(
            gr3, region, inner_value=value, outer_value=-value)


class Tvdflag(Prop):
    """Class for writing tvd.prop file, which specify horizontal regions 
       where upwind or TVD/TVD^2 is used based on the element property values
       (0: upwind; 1: TVD/TVD^2).
    """
    @classmethod
    def default(cls, hgrid):
        return cls.constant(hgrid, 1)

    @classmethod
    def from_geometry(cls, gr3: Gr3, region: Union[Polygon, MultiPolygon],
                      value: int = 1):
        if value not in [1, -1]:
            raise ValueError('Argument value must be 1 or -1.')
        return super(Tvdflag, cls).from_geometry(
            gr3, region, inner_value=value, outer_value=-value)


def reg2multipoly(file):
    with open(file) as f:
        f.readline()
        npoly = int(f.readline())
        polygons = []
        for _ in range(npoly):
            line = f.readline().split()
            if len(line) == 0:
                break
            nvrt = int(line[0])
            tflag = int(line[1])
            exterior = []
            for _ in range(nvrt):
                exterior.append(tuple(list(map(float, f.readline().split()))))
            interiors = []
            line = f.readline().split()
            if len(line) != 0:
                while tflag == 1:
                    nvrt = int(line[0])
                    tflag = int(line[1])
                    interior = []
                    for _ in range(nvrt):
                        interior.append(tuple(list(map(float, f.readline().split()))))
                    interiors.append(interior)
            polygons.append(Polygon(exterior, interiors))
    return MultiPolygon(polygons)
