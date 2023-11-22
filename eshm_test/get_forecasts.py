import os
from os.path import join, dirname
from oq2csep.cmd import main
from oq2csep.sm_lib import parse_logictree_files

dir_script = dirname(__file__)


if __name__ == '__main__':

    region = join(dir_script, 'regions', 'region_final.txt')

    output_dir = join(dir_script, 'forecasts')
    os.makedirs(output_dir, exist_ok=True)

    eshm13_dir = join(dir_script, 'eshm13', 'SHARE_OQ_input_20140807')
    eshm13_lt = join(eshm13_dir, 'source_model_logic_tree.xml')
    eshm13_branches = parse_logictree_files(eshm13_lt)
    branch_names = ['eshm13_as', 'eshm13_fsbg', 'eshm13_seifa']

    for i, (key, branch_path) in enumerate(eshm13_branches[1].items()):
        branch = [join(eshm13_branches[0], i) for i in branch_path]
        main.project(branch,
                     region=region,
                     dest=join(output_dir, branch_names[i] + '.csv'),
                     plot=True)

    eshm20_dir = join(dir_script, 'eshm20', 'oq_computational',
                      'oq_configuration_eshm20_v12e_region_main')
    eshm20_lt = join(eshm20_dir,
                     'source_model_logic_tree_eshm20_model_v12e.xml')
    eshm20_dir, eshm20_branches = parse_logictree_files(eshm20_lt)

    for key, model_paths in eshm20_branches.items():
        branch = [join(eshm20_dir, i) for i in model_paths]
        main.project(branch,
                     region=region,
                     dest=join(output_dir, key + '.csv'),
                     plot=True)

