import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import *
import numpy as np
import matplotlib.pyplot as plt
import cmath
from matplotlib.collections import LineCollection
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from FermiLevelPinning_HAN import *

ZnO = Make_defective_model()
ZnO.set_bulk_property(715.19e-24,-286.79394896, 0.7586, 2.5677, 'DOSCAR_mode', 'DOSCAR_ZnO_PBEU') # volume & total energy & vbm & cbm & DOSCAR mode & DOSCAR filename or effective masses
chemical_potential_ZnO_Orich = [-4.054142985, -4.92616616] # Zn O
chemical_potential_ZnO_Znrich = [-0.334123815, -8.64618533] # Zn O
create_anhil_oxygen_vacancy = [0.0, -1.0]
create_anhil_zinc_vacancy = [-1.0, 0.0]
ZnO.set_defect_class('oxygen_vacancy', chemical_potential_ZnO_Orich, create_anhil_oxygen_vacancy) # defect class & chemical potential & create_anhil_info
ZnO.set_defect_class('zinc_vacancy', chemical_potential_ZnO_Orich, create_anhil_zinc_vacancy)

ZnO.set_defect_detail('oxygen_vacancy', -277.18429135, 0) # defect class & total energy & charge state
ZnO.set_defect_detail('oxygen_vacancy', -273.74553174, -1)
ZnO.set_defect_detail('oxygen_vacancy', -270.06781731, -2)
ZnO.set_defect_detail('oxygen_vacancy', -279.07155048, 1)
ZnO.set_defect_detail('oxygen_vacancy', -281.79138418, 2)

ZnO.set_defect_detail('zinc_vacancy', -279.82989041, 0)
ZnO.set_defect_detail('zinc_vacancy', -279.29412802, -1)
ZnO.set_defect_detail('zinc_vacancy', -278.59102532, -2)
ZnO.set_defect_detail('zinc_vacancy', -280.32489306, 1)
ZnO.set_defect_detail('zinc_vacancy', -280.79925758, 2)
#ZnO.print_info()


#################3 Draw theromdynamic charge transition level ##############
#ZnO.draw_thermodynamic_charge_transition_level()

############################################################################
cal_ZnO = calculate_Fermi_level(ZnO, 300) # model name & temperature
cal_ZnO.set_additional_info('oxygen_vacancy',  0, 4.47434e22, 1, 1) # defect class & charge state & N_site & spin_degeneracy & structure_degeneracy
cal_ZnO.set_additional_info('oxygen_vacancy', -1, 4.47434e22, 2, 1)
cal_ZnO.set_additional_info('oxygen_vacancy', -2, 4.47434e22, 1, 1)
cal_ZnO.set_additional_info('oxygen_vacancy',  1, 4.47434e22, 2, 1)
cal_ZnO.set_additional_info('oxygen_vacancy',  2, 4.47434e22, 1, 1)

cal_ZnO.set_additional_info('zinc_vacancy',  0, 4.47434e22, 1, 1)
cal_ZnO.set_additional_info('zinc_vacancy', -1, 4.47434e22, 2, 1)
cal_ZnO.set_additional_info('zinc_vacancy', -2, 4.47434e22, 1, 1)
cal_ZnO.set_additional_info('zinc_vacancy',  1, 4.47434e22, 2, 1)
cal_ZnO.set_additional_info('zinc_vacancy',  2, 4.47434e22, 1, 1)

#print cal_ZnO.defect_info
#print cal_ZnO.detail_info_defect
#cal_ZnO.check_everything_okay()
#print cal_ZnO.ready
cal_ZnO.main_routine(3.0, 1000, 0.01, 1e-15) # initial_Ef, total number of scf loops, delta_Ef, convergence_criteria
#cal_ZnO.calculate_free_carriers(1.2)
#print cal_ZnO.calculate_defect_concentration(0, 2.5)

