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
#from numeric import *

class Make_defective_model:
    def __init__(self, vbm, cbm):
        self.vbm = vbm  # vbm should be absolute value from DFT calculations not ralative values as compared to vbm, Fermi level ...
        self.cbm = cbm  # cbm should be absolute value from DFT calculations not ralative values as compared to vbm, Fermi level ...
        self.band_gap = cbm - vbm

    def set_bulk_property(self, Etotal, DOS_info, additional_info):
        self.E_bulk = Etotal

        if DOS_info == 'DOSCAR':
            print 'Use real DOS from bulk supercell'
            print 'Read DOSCAR = ' + str(additional_info)
            self.doscar_name = str(additional_info)
            self.read_DOSCAR()
        else:
            print 'Use parabolic DOS'
            if not len(additional_info) == 2:
                print 'Length of additional info is wrong ; please check it'
                return 0
            mass_elec, mass_hole = float(additional_info[0]), float(additional_info[1])
            print 'Electron mass = ' + str(mass_elec) + '  /  ' + 'Hole mass = ' + str(mass_hole)


    def read_DOSCAR(self):
        f = open(self.doscar_name, 'r')
        fbuffer = f.readlines()
        f.close()

        grid = int(fbuffer[5].split()[2])
        Ef = float(fbuffer[5].split()[3])

        energy = [] ; total_dos = []

        for i in range(grid):
            tempf = fbuffer[6+i].split()
            energy.append(float(tempf[0]))
            total_dos.append(float(tempf[1]))
        energy = np.array(energy) ; total_dos = np.array(total_dos)

        max_energy = np.max(energy) ; min_energy = np.min(energy) ; delta_energy = (max_energy - min_energy) / (grid - 1.0)

        self.bulk_energy = energy
        self.bulk_totalDOS = total_dos

        return 0

    def set_intrinsic_defect(self, defect_name, E_defect, chemical_potential):
        

    def draw_thermodynamic_charge_transition_level(self):
        return 0

             

if __name__ == "__main__": 
    ZnO = Make_defective_model(1.2, 2.0)
    ZnO.set_bulk_property(-100, 'DOSCAR', 'DOSCAR')
