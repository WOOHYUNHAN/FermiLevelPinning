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

#from numeric import *

class Make_defective_model:
    def __init__(self):
        self.defect_info = []

    def set_bulk_property(self, Etotal, vbm, cbm, DOS_info, additional_info):
        self.E_bulk = Etotal
        self.vbm = vbm # vbm should be absolute value from DFT calculations not ralative values as compared to vbm, Fermi level ...
        self.cbm = cbm  # cbm should be absolute value from DFT calculations not ralative values as compared to vbm, Fermi level ...
        self.band_gap = cbm - vbm

        if DOS_info == 'DOSCAR':
            print 'Use real DOS from bulk supercell'
            print 'Read DOSCAR & filename = ' + str(additional_info)
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

    def set_defect_class(self, defect_class, chemical_potential, create_anhil_info):
        if len(chemical_potential) != len(create_anhil_info):
            print 'ERROR: please check both chemical potential and create_anhil_info'
            return 0

        find_class = False
        for i in range(len(self.defect_info)):
            if defect_class == self.defect_info[i][0]:
                find_class = True
                index_class = i
        
        if find_class:
            print 'ERROR: You already have the same defect class; ' + str(self.defect_info[index_class][0])
            return 0

        temp_class = [defect_class, chemical_potential, create_anhil_info]
        self.defect_info.append(temp_class)


        return 0

    def clear_defect_class(self):

        return 0

    def set_defect_detail(self, defect_name, E_defect, charge_state, correct=0.0):
        '''
        defect formation energy:
        ref: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.67.2339

        defect name is very important !!!!
        Before draw charge transition level please run print info function
        '''

        if len(self.defect_info) == 0:
            print 'ERROR: defect info is empty ; set defect class first and then set defect detail'
            return 0
        else:
            temp_name = ''
            for i in range(len(self.defect_info)):
                temp_name += str(self.defect_info[i][0]) + '; '
            print 'You have defect classes as shown below; total number = ' + str(len(self.defect_info))
            print temp_name

        find_class = False
        for i in range(len(self.defect_info)):
            if defect_name == self.defect_info[i][0]:
                find_class = True
                index_class = i
        
        if not find_class:
            print 'ERROR: Your defect name do not concide with any class declared before'
            return 0

        temp_defect_info = self.defect_info[index_class]
        chemical_info = np.array(temp_defect_info[1]) * np.array(temp_defect_info[2])
        charge_info = charge_state * self.vbm
        diff_energy = E_defect - self.E_bulk - np.sum(chemical_info) #+ charge_info + correct
        defect_info = [charge_state, diff_energy]

        find_charge = False
        stored_info_num = len(self.defect_info[index_class]) - 3
        #print stored_info_num
        if not stored_info_num == 0:
            for i in range(stored_info_num):
                if float(charge_state) == self.defect_info[index_class][i+3][0]:
                    find_charge = True
                    index_charge = i

        if find_charge:
            print 'You already have the same charge state; So the new info will replace the old info'
            print 'Charge state: ' + str(charge_state)
            self.defect_info[index_class][index_charge+3][0] = charge_state
            self.defect_info[index_class][index_charge+3][1] = diff_energy
        else:
            print 'New charge state is introduced '
            print 'Charge state: ' + str(charge_state)
            self.defect_info[index_class].append(defect_info)

        self.class_num = len(self.defect_info)


        

    def draw_thermodynamic_charge_transition_level(self, min_E_from_VBM, max_E_from_CBM, Estep):
        energy_step = np.array([-1*min_E_from_VBM + self.vbm + (self.cbm - self.vbm + max_E_from_CBM + min_E_from_VBM) * i / float(Estep) for i in range(Estep)])
        #print len(energy_step)

        min_energy = []
        name = []
        color = ['black', 'red', 'green', 'blue', 'pink']
        
        for i in range(self.class_num):
            name.append(str(self.defect_info[i][0]))
            min_energy.append([])
            for j in range(len(energy_step)):
                temp = []
                for k in range(len(self.defect_info[i])-3):
                    charge_state = float(self.defect_info[i][k+3][0])
                    diff_energy = float(self.defect_info[i][k+3][1])
                    temp.append(diff_energy + charge_state * energy_step[j])
                min_energy[i].append(np.min(np.array(temp)))


        gs = gridspec.GridSpec(1, 1, width_ratios=[1], height_ratios=[1])
        ax1 = plt.subplot(gs[0,0])

        for i in range(self.class_num):
            ax1.plot(energy_step, min_energy[i], '-', color=color[i], label=r''+str(name[i])+'')

        ax1.axvline(x=self.vbm) ; ax1.axvline(x=self.cbm)

        handles, labels = ax1.get_legend_handles_labels()
        leg=ax1.legend(handles[::-1], labels[::-1], frameon=True ,loc='best',numpoints=1,handletextpad=0.5,borderpad=0.05, ncol=1, labelspacing=0.3, handlelength=2, prop={'size':10})

        ax1.set_xlabel(r'Energy (eV)')
        ax1.set_ylabel(r'Formation energy (eV)')
        ax1.set_xlim(self.vbm - min_E_from_VBM, self.cbm + max_E_from_CBM)
        #ax1.set_ylim(min_y-0.01,-5.0)
        plt.savefig("chargetransitionlevel.png")
        plt.show()

        return 0


    def print_info(self):
        print '==============print info=============='
        print 'VBM, CBM, band_gap = ' + str(self.vbm) + ', ' + str(self.cbm) + ', ' + str(self.band_gap) + ', '
        print 'defect information'
        print 'defect class, chemical potential, create_anhil_info, charge state, diff energy'
        for i in range(len(self.defect_info)):
            temp_info = self.defect_info[i]
            print str(temp_info[0]) + ' ' + str(temp_info[1]) + ' ' + str(temp_info[2]) + ' ' 
            for j in range(len(temp_info)-3):
                temp_line = str(temp_info[j+3][0]) + ' ' + str(temp_info[j+3][1])
                print temp_line
        print '==============print info==============' 

if __name__ == "__main__": 
    ZnO = Make_defective_model()
    ZnO.set_bulk_property(-286.79394896, 0.7586, 2.5677, 'DOSCAR', 'DOSCAR_ZnO_PBEU')
    chemical_potential_ZnO_Orich = [-4.054142985, -4.92616616] # Zn O
    chemical_potential_ZnO_Znrich = [-0.334123815, -8.64618533] # Zn O
    create_anhil_oxygen_vacancy = [0.0, -1.0]
    create_anhil_zinc_vacancy = [-1.0, 0.0]
    ZnO.set_defect_class('oxygen_vacancy', chemical_potential_ZnO_Orich, create_anhil_oxygen_vacancy)
    ZnO.set_defect_class('zinc_vacancy', chemical_potential_ZnO_Orich, create_anhil_zinc_vacancy)
    
    ZnO.set_defect_detail('oxygen_vacancy', -277.18429135, 0)
    ZnO.set_defect_detail('oxygen_vacancy', -273.74553174, -1)
    ZnO.set_defect_detail('oxygen_vacancy', -270.06781731, -2)
    ZnO.set_defect_detail('oxygen_vacancy', -279.07155048, 1)
    ZnO.set_defect_detail('oxygen_vacancy', -281.79138418, 2)
    
    ZnO.set_defect_detail('zinc_vacancy', -279.82989041, 0)
    ZnO.set_defect_detail('zinc_vacancy', -279.29412802, -1)
    ZnO.set_defect_detail('zinc_vacancy', -278.59102532, -2)
    ZnO.set_defect_detail('zinc_vacancy', -280.32489306, 1)
    ZnO.set_defect_detail('zinc_vacancy', -280.79925758, 2)

    ZnO.draw_thermodynamic_charge_transition_level(0.2, 0.2, 400)
    #print ZnO.class_num

    #ZnO.set_defect_detail('oxygen_vacancy', -20.0, -2)
    ZnO.print_info()