import os
from os.path import join, dirname
from hazard2csep.cmd import main

dir_script = dirname(__file__)


if __name__ == '__main__':

    output_dir = join(dir_script, 'output')
    os.makedirs(output_dir, exist_ok=True)

    filepath = join(dir_script, 'fs_bg_source_model.xml')
    forecast_path = join(output_dir, 'forecast.csv')
    main.project(filepath,
                 dest=forecast_path,
                 plot=True)


