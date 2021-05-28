from abc import ABC, abstractmethod
from enum import Enum
from functools import lru_cache
import pathlib
import subprocess
import tempfile

import numpy as np

from pyschism.mesh.hgrid import Hgrid


class VgridType(Enum):

    LSC2 = 1
    SZ = 2

    @classmethod
    def _missing_(cls, name):
        raise ValueError(f'ivcor={name} is not a valid vgrid type.')


class Vgrid(ABC):

    @abstractmethod
    def __str__(self):
        raise NotImplementedError

    @staticmethod
    def default():
        return SZ.default()

    @classmethod
    def from_binary(cls, hgrid, binary='gen_vqs'):
        _tmpdir = tempfile.TemporaryDirectory()
        tmpdir = pathlib.Path(_tmpdir.name)
        hgrid = Hgrid.open(hgrid, crs='EPSG:4326')
        hgrid.write(tmpdir / 'hgrid.gr3')
        subprocess.check_call([binary], cwd=tmpdir)
        return cls.open(tmpdir / 'vgrid.in')

    @staticmethod
    def open(path):
        '''
        Based on:
        https://github.com/wzhengui/pylibs/blob/master/Utility/schism_file.py
        '''
        with open(path) as f:
            return VgridTypeDispatch[VgridType(
                int(f.read().strip().split()[0])).name].value.open(path)

    def write(self, path, overwrite=False):
        path = pathlib.Path(path)
        if path.is_file() and not overwrite:
            raise Exception(
                'File exists, pass overwrite=True to allow overwrite.')

        with open(path, 'w') as f:
            f.write(str(self))

    @property
    def ivcor(self):
        return VgridType[self.__class__.__name__].value

    @property
    @abstractmethod
    def nvrt(self):
        raise NotImplementedError

    @lru_cache(maxsize=1)
    def is2D(self):
        if str(self) == str(SZ.default()):
            return True
        return False

    def is3D(self):
        return ~self.is2D()


class LSC2(Vgrid):

    def __init__(self, sigma):
        self.sigma = sigma

    def __str__(self):
        f = [
            f'{self.ivcor}',
            f'{self.nvrt}',
        ]
        for i, row in enumerate(self.sigma):
            kbp = int((row == -1).sum())
            line = [
                f'{i+1}'.rjust(11),
                f'{kbp}'.rjust(11),
                7*' ',
                '-1.000000',
            ]
            for value in row:
                if value != -1:
                    line.append(7*' ')
                    line.append(f'{value:6f}')

            f.append(' '.join(line))
        return '\n'.join(f)

    @classmethod
    def open(cls, path):

        path = pathlib.Path(path)

        with open(path) as f:
            lines = f.readlines()

        ivcor = int(lines[0].strip().split()[0])
        if ivcor != 1:
            raise TypeError(f'File {path} is not an LSC2 grid (ivcor != 1).')

        nvrt = int(lines[1].strip().split()[0])

        kbp = np.array([int(i.split()[1])-1 for i in lines[2:]])

        sigma = -np.ones([len(kbp), nvrt])

        for i, line in enumerate(lines[2:]):
            sigma[i, kbp[i]:] = np.array(
                line.strip().split()[2:]).astype('float')

        return cls(sigma)

    @property
    def nvrt(self):
        return self.sigma.shape[1]


class SZ(Vgrid):

    def __init__(self, h_s, ztot, h_c, theta_b, theta_f, sigma):
        self.h_s = h_s
        self.ztot = ztot
        self.h_c = h_c
        self.theta_b = theta_b
        self.theta_f = theta_f
        self.sigma = sigma

    def __str__(self):
        f = [
            f'{self.ivcor:d} !ivcor',
            f'{self.nvrt:d} {self.kz:d} {self.h_s:G} '
            '!nvrt, kz (# of Z-levels); h_s '
            ' (transition depth between S and Z)',
            'Z levels',
        ]
        # print(self.ztot)
        # exit()
        for row in enumerate(self.ztot):
            # f.append(f'{int(row[0]):d} {row[1]:G}')
            # line = [f'']
            for i, x in enumerate(row):
                print(i, x)
                # line.append(f' {i+1:d} {x:G}')
                f.append(f'{i+1:d} {x:G}')
        f.extend([
            'S levels',
            f'{self.h_c:G} {self.theta_b:G} {self.theta_f:G} '
            ' !h_c, theta_b, theta_f',
            ])
        for row in enumerate(self.sigma):
            # f.append(f'{int(row[0])} {row[1]:G}')
            for i, x in enumerate(row):
                # line.append(f' {i+1:d} {x:G}')
                f.append(f'{i+1:d} {x:G}')
        return '\n'.join(f)

    @classmethod
    def open(cls, path):

        path = pathlib.Path(path)

        with open(path) as f:
            lines = f.readlines()

        ivcor = int(lines[0].strip().split()[0])
        if ivcor != 2:
            raise TypeError(f'File {path} is not an SZ grid (ivcor != 2).')

        nvrt = int(lines[1].strip().split()[0])

        kz, h_s = lines[1].strip().split()[1:3]
        kz = int(kz)
        h_s = float(h_s)

        # read z grid
        ztot = []
        irec = 2
        for i in np.arange(kz):
            irec = irec+1
            ztot.append(lines[irec].strip().split()[1])
        print(ztot)
        ztot = np.array(ztot).astype('float')
        print(ztot)
        # read s grid
        sigma = []
        irec = irec+2
        nsigma = nvrt - kz+1
        h_c, theta_b, theta_f = np.array(
            lines[irec].strip().split()[:3]).astype('float')
        for i in np.arange(nsigma):
            irec = irec + 1
            sigma.append(lines[irec].strip().split()[1])
        sigma = np.array(sigma).astype('float')
        return cls(h_s, ztot, h_c, theta_b, theta_f, sigma)

    @classmethod
    def default(cls):
        ztot = np.array([-1.e6])
        sigma = np.array([-1, 0.])
        return cls(1.e6, ztot, 40., 1., 1.e-4, sigma)

    @property
    def kz(self):
        return self.ztot.shape[0]

    @property
    def nvrt(self):
        return self.sigma.shape[0]


class VgridTypeDispatch(Enum):

    LSC2 = LSC2
    SZ = SZ
