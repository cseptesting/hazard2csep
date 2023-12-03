import os

from hazard2csep import region_lib
from hazard2csep import forecast_lib
from hazard2csep import sm_lib
from hazard2csep import logger

# Handle proj installation issue to point bin libraries
try:
    venv = os.environ['CONDA_PREFIX']
    proj_path = os.path.join(venv, 'share')

except KeyError:  # using python-venv
    venv = os.environ['VIRTUAL_ENV']
    proj_path = os.path.join(venv, 'share')

os.environ['PROJ_LIB'] = f'{proj_path}/proj'


__version__ = '0.1.1'
