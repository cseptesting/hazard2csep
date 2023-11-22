import os
from os.path import join, dirname
from oq2csep.cmd import main


dir_script = dirname(__file__)

eshm13_sm = join(dir_script, 'eshm13', 'SHARE_OQ_input_20140807')
ESHM13_PATH = {'seifa': join(eshm13_sm, 'seifa_model.xml'),
               'fsbg': join(eshm13_sm, 'faults_backg_source_model.xml'),
               'as': join(eshm13_sm, 'area_source_model.xml')}


eshm20_sm = join(dir_script, 'eshm20', 'oq_computational',
                 'oq_configuration_eshm20_v12e_region_main/source_models')
eshm20_subd_interface = [join(eshm20_sm, 'interface_v12b', i)
                         for i in ['CaA_IF2222222_M40.xml',
                                   'CyA_IF2222222_M40.xml',
                                   'GiA_IF2222222_M40.xml',
                                   'HeA_IF2222222_M40.xml']]
ESHM20_PATH = {
    'as': [join(eshm20_sm, 'asm_v12e',
                'asm_ver12e_winGT_fs017_lo_abgrs_maxmag_low.xml'),
           join(eshm20_sm, 'asm_v12e',
                'asm_ver12e_winGT_fs017_twingr.xml'),
           *eshm20_subd_interface],
    'fsbg': [
           join(eshm20_sm, 'fsm_v09', 'fs_ver09e_model_aGR_SRL_ML_fMthr.xml'),
           join(eshm20_sm, 'ssm_v09',
                'seis_ver12b_fMthr_asm_ver12e_winGT_fs017_agbrs_point.xml'),
           *eshm20_subd_interface]}


if __name__ == '__main__':

    output_dir = join(dir_script, 'regions')
    os.makedirs(output_dir, exist_ok=True)
    #
    reg_eshm13_as = join(output_dir, 'region_eshm13_as.txt')
    reg_eshm13_fsbg = join(output_dir,'region_eshm13_fsbg.txt')
    reg_eshm13_seifa = join(output_dir, 'region_eshm13_seifa.txt')
    reg_eshm13 = join(output_dir, 'region_eshm2013.txt')
    #
    main.region(ESHM13_PATH['as'], dest=reg_eshm13_as,
                plot=True, fill=True)
    main.region(ESHM13_PATH['fsbg'], dest=reg_eshm13_fsbg,
                plot=True, fill=True)
    main.region(ESHM13_PATH['seifa'], dest=reg_eshm13_seifa,
                plot=True,fill=True)
    main.region([reg_eshm13_as, reg_eshm13_fsbg, reg_eshm13_seifa],
                dest=reg_eshm13, intersect=True, plot=True)

    reg_eshm20_as = join(output_dir, 'region_eshm20_as.txt')
    reg_eshm20_fsbg = join(output_dir, 'region_eshm20_fsbg.txt')
    reg_eshm20 = join(output_dir, 'region_eshm20.txt')
    main.region(ESHM20_PATH['as'], dest=reg_eshm20_as, plot=True,
                fill=True)
    main.region(ESHM20_PATH['fsbg'], dest=reg_eshm20_fsbg, plot=True,
                fill=True)
    main.region([reg_eshm20_as, reg_eshm20_fsbg],
                dest=reg_eshm20, intersect=True, plot=True)

    reg_final = join(output_dir, 'region_final.txt')
    main.region([reg_eshm13, reg_eshm20],
                dest=reg_final, intersect=True, plot=True)


