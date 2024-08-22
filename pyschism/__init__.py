from importlib import util

# from pyschism.param.param import Param
from pyschism.stations import Stations
from pyschism.driver import ModelDriver, ModelConfig
from pyschism.plotting import read_gr3_file, plot_gr3  # <-- Add this line

__all__ = ['ModelConfig', 'ModelDriver', 'Stations', 'read_gr3_file', 'plot_gr3']

if util.find_spec("colored_traceback") is not None:
    import colored_traceback  # type: ignore[import]
    colored_traceback.add_hook(always=True)
