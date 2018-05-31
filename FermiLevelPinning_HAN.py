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

    def set_bulk_property(self, vol, Etotal, vbm, cbm, DOS_info, additional_info):
        self.E_bulk = Etotal
        self.vbm = vbm # vbm should be absolute value from DFT calculations not ralative values as compared to vbm, Fermi level ...
        self.cbm = cbm  # cbm should be absolute value from DFT calculations not ralative values as compared to vbm, Fermi level ...
        self.band_gap = cbm - vbm
        self.DOS_info = DOS_info
        self.bulk_vol = vol

        if self.DOS_info == 'DOSCAR_mode':
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

        self.bulk_dos_energy = energy
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

        ax1.axvline(x=self.vbm) ; ax1.axvline(x=self.cbm) ; ax1.axhline(y=0.0)

        handles, labels = ax1.get_legend_handles_labels()
        leg=ax1.legend(handles[::-1], labels[::-1], frameon=True ,loc='best',numpoints=1,handletextpad=0.5,borderpad=0.05, ncol=1, labelspacing=0.3, handlelength=2, prop={'size':10})

        ax1.set_xlabel(r'Energy (eV)')
        ax1.set_ylabel(r'Formation energy (eV)')
        ax1.set_xlim(self.vbm - min_E_from_VBM, self.cbm + max_E_from_CBM)
        #ax1.set_ylim(min_y-0.01,-5.0)
        plt.savefig("charge_transition_level.png")
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

class calculate_Fermi_level:
    def __init__(self, model, temperature):
        self.temperature = temperature
        self.boltzmann = 8.617343 * 1e-5
        self.defect_info = model.defect_info #### [class, chemical_potential, create_anil_info, [charge_state1, diff_energy1], [charge_state2, diff_energy2]...]
        self.class_num = model.class_num
        self.E_bulk = model.E_bulk
        self.vbm = model.vbm 
        self.cbm = model.cbm
        self.band_gap = model.band_gap
        self.bulk_dos_energy = model.bulk_dos_energy
        self.bulk_totalDOS = model.bulk_totalDOS
        self.bulk_vol = model.bulk_vol
        self.DOS_info = model.DOS_info
        self.detail_info_defect = []
        for i in range(self.class_num):
            self.detail_info_defect.append([])
            num_charge = len(self.defect_info[i])-3
            #print num_charge
            for j in range(num_charge):
                self.detail_info_defect[i].append([])
        self.ready = False


    def set_additional_info(self, defect_class_name, charge_state, N_site, spin_degeneracy, struc_degeneracy):
        find_class, find_charge = False, False

        for i in range(len(self.defect_info)):
            if defect_class_name == self.defect_info[i][0]:
                find_class = True
                index_class = i
        
        if not find_class:
            print 'ERROR: Your input for defect class is wrong; please check again ' + str(defect_class_name)
            return 0

        for i in range(len(self.defect_info[index_class])-3):
            if charge_state == self.defect_info[index_class][i+3][0]:
                find_charge = True
                index_charge = i

        if not find_charge:
            print 'ERROR: Your input for defect charge state is wrong; please check again ' + str(charge_state)
            return 0
        #print self.detail_info_defect
        #[N_site, spin_degeneracy, struc_degenery]
        ######## info ########
        # N_site = possible equivalent defect sites per volume 
        # spin_degeneracy = possible equivalent spin configurations (two-fold degerenated levels, one electron == 2, two electrons == 1)
        # struc_degeneracy = possible equivalent defect structures (point defct == 1)
        ######################
        self.detail_info_defect[index_class][index_charge].append(N_site)
        self.detail_info_defect[index_class][index_charge].append(spin_degeneracy)
        self.detail_info_defect[index_class][index_charge].append(struc_degeneracy)
        return 0

    def check_everything_okay(self):
        for i in range(self.class_num):
            for j in range(len(self.defect_info[i])-3):
                if len(self.detail_info_defect[i][j]) == 0:
                    print "ERROR: some part of detail defect info is empty"
                    print self.defect_info[i][0]
                    return 0
        print 'Everything is okay; you are ready to calculate Fermi level'
        self.ready = True
        return 0

    def calculate_defect_concentration(self, defect_class, Ef):
        # step 1: find most dominant charge state
        if defect_class > self.class_num - 1:
            print "ERROR: defect class is larger than the number of defect class"
            return 0

        num_charge = len(self.defect_info[defect_class]) - 3
        temp_energy = []
        for i in range(num_charge):
            charge_state = float(self.defect_info[defect_class][i+3][0])
            diff_energy = float(self.defect_info[defect_class][i+3][1])
            temp_energy.append(diff_energy + charge_state * Ef)

        min_index = np.argmin(np.array(temp_energy))
        #print min_index
        charge_state = self.defect_info[defect_class][min_index+3][0]
        formation_energy = float(temp_energy[min_index])
        #print charge_state, formation_energy

        # step 2: calculate defect concentrations

        N_site, spin_deg, struc_deg = self.detail_info_defect[defect_class][min_index]
        #print N_site, spin_deg, struc_deg

        temp_part = np.exp(-1.0 * formation_energy / (self.boltzmann * self.temperature))

        concentration = N_site * spin_deg * struc_deg * temp_part

        return charge_state, concentration

    def calculate_free_carriers(self, Ef):
        if self.DOS_info == 'DOSCAR_mode':
            ## READ DOSCAR
            temp_hole = 0
            temp_elec = 0
            max_energy = np.max(self.bulk_dos_energy) ; min_energy = np.min(self.bulk_dos_energy)  ; ngrid = len(self.bulk_dos_energy) ; dE = (max_energy - min_energy) / (ngrid - 1.0)
            DOS_ptype = np.zeros(ngrid) ; DOS_ntype = np.zeros(ngrid)
            E_ptype = Ef - self.bulk_dos_energy ; E_ntype = self.bulk_dos_energy - Ef
            for i in range(ngrid):
                if self.bulk_dos_energy[i] > self.cbm:
                    DOS_ntype[i] = self.bulk_totalDOS[i]
                elif self.bulk_dos_energy[i] < self.vbm:
                    DOS_ptype[i] = self.bulk_totalDOS[i]
                else:
                    pass 
        else:
            ## USE PARABOLIC DOS
            pass

        p_carrier = np.sum(DOS_ptype* dE / (np.exp(E_ptype/(self.boltzmann*self.temperature)) + 1.0)) / self.bulk_vol
        n_carrier = np.sum(DOS_ntype* dE / (np.exp(E_ntype/(self.boltzmann*self.temperature)) + 1.0)) / self.bulk_vol
        #print Ef, n_carrier, p_carrier

        return n_carrier, p_carrier

    def main_routine(self, initial_Ef, total_steps, delta, conv_criteria):
        self.check_everything_okay()
        if self.ready:
            print "ERROR: please check set addition info function"

        Ef = initial_Ef
        step = 0
        LRUN = True
        class_info = []
        temp_Ef = []
        print 'step' + '\t' + 'Fermi level' + '\t' + 'Electron' + '\t' + 'Hole'

        while step < total_steps and LRUN:
            electron = 0
            hole = 0
            ################# Defect concentration #########################
            for i in range(self.class_num):
                charge, concent = self.calculate_defect_concentration(i, Ef)
                class_info.append([charge, concent])
                if charge > 0:
                    ### donor ###
                    hole += abs(charge) * concent
                    #pass
                elif charge < 0:
                    ### acceptor ###
                    electron += abs(charge) * concent
                    #pass
                else:
                    ### netural ###
                    pass

            ################## free carriers from CB and VB ################
            n_carrier, p_carrier = self.calculate_free_carriers(Ef)

            electron += n_carrier
            hole += p_carrier


            ################## self-consistent loop ########################
            templine = str(step) + '\t' + str(Ef) + '\t' + str(electron) + '\t' + str(hole)  + '\t' + str(n_carrier)  + '\t' + str(p_carrier)
            for i in range(self.class_num):
                charge, concent = self.calculate_defect_concentration(i, Ef)
                templine += '\t' + str(concent)
            print templine

            temp_Ef.append(float(Ef))

            if abs(electron - hole) < conv_criteria:
                print "Converge"
                LRUN = False
            else:
                if electron > hole:
                    ### decrase Ef
                    Ef -= delta
                elif electron < hole:
                    ### increase Ef
                    Ef += delta
                else:
                    pass
            step += 1
        
        gs = gridspec.GridSpec(1, 1, width_ratios=[1], height_ratios=[1])
        ax1 = plt.subplot(gs[0,0])

        temp_step = [i for i in range(len(temp_Ef))]

        ax1.plot(temp_step, temp_Ef, '-', color='black')

        ax1.set_xlabel(r'Step')
        ax1.set_ylabel(r'Fermi level (eV)')
        ax1.set_xlim(0, len(temp_Ef))
        #ax1.set_ylim(min_y-0.01,-5.0)
        plt.savefig("FermiLevelPinning.png")
        #plt.show()


        return 0


if __name__ == "__main__": 
    ZnO = Make_defective_model()
    ZnO.set_bulk_property(715.19e-24,-286.79394896, 0.7586, 2.5677, 'DOSCAR_mode', 'DOSCAR_ZnO_PBEU')
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

    ###########################################################################
    #ZnO.draw_thermodynamic_charge_transition_level(0.2, 0.2, 400)
    #print ZnO.class_num

    #ZnO.set_defect_detail('oxygen_vacancy', -20.0, -2)
    #ZnO.print_info()
    ############################################################################
    cal_ZnO = calculate_Fermi_level(ZnO, 300)
    cal_ZnO.set_additional_info('oxygen_vacancy',  0, 4.47434e22, 1, 1)
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
    cal_ZnO.main_routine(3.0, 1000, 0.01, 1e-15)
    #cal_ZnO.calculate_free_carriers(1.2)
    #print cal_ZnO.calculate_defect_concentration(0, 2.5)

