import os
import pandas as pd
import numpy as np
import shutil
from scipy import optimize
import tensorflow as tf
from openpyxl import Workbook
from tensorflow.keras.layers import Dense, Activation, Lambda, Dropout
from sklearn.model_selection import train_test_split
from tensorflow.keras.regularizers import l2
from tensorflow.keras.optimizers import Adam, RMSprop, Adadelta, Adagrad, Adamax
from tensorflow.keras.models import Sequential
import keras_tuner as kt
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.metrics import MeanAbsoluteError
from tensorflow.keras.layers import LeakyReLU
LeakyReLU = LeakyReLU(negative_slope=0.1)

class IsotopomerAnalysis:

    def __init__(self):
        self.nmr_multiplets = pd.DataFrame()
        self.nmr_tocsy_multiplets = pd.DataFrame()
        self.gcms_data = pd.DataFrame()
        self.lcms_data = pd.DataFrame()
        self.nmr1d_data = pd.DataFrame()
        self.ver = '0.9.23'
        self.nat_abundance = 1.07  # [%]
        self.use_hsqc_multiplet_data = True
        self.use_gcms_data = True
        self.use_nmr1d_data = False
        self.gcms_scaling = 1.0
        self.hsqc_scaling = 1.0
        self.lcms_scaling = 1.0
        self.tocsy_scaling = 1.0
        self.metabolites = []
        self.n_exps = 0
        self.fit_isotopomers = {}
        self.isotopomer_percentages = {}
        self.nmr_isotopomers = {}
        self.nmr_isotopomer_percentages = {}
        self.gcms_percentages = {}
        self.nmr1d_percentages = {}
        self.exp_multiplets = {}
        self.exp_multiplet_percentages = {}
        self.exp_gcms = {}
        self.exp_nmr1d = {}
        self.hsqc = {}
        self.fitted_isotopomers = {}
        self.fitted_isotopomer_percentages = {}
        self.fitted_multiplets = {}
        self.fitted_multiplet_percentages = {}
        self.fitted_gcms_percentages = {}
        self.fitted_nmr1d_percentages = {}
        self.exp_hsqc_isos = {}
        self.hsqc_multiplets = {}
        self.hsqc_multiplets2 = {}
        self.hsqc_multiplet_percentages = {}
        self.n_bonds = {}
        self.chi2 = 0
        self.current_metabolite = ''
        self.current_experiment = 0
        self.results = []
        # end __init__

    def __str__(self):  # pragma: no cover
        r_string = '______________________________________________________________________________________\n'
        r_string += '\nMetaboLabPy Isotopomer Data Analysis (v. ' + self.ver + ')\n'
        r_string += '______________________________________________________________________________________\n\n'
        return r_string
        # end __str__

    def export_data(self, file_name=''):
        wb = Workbook()
        wb.remove(wb.active)
        for k in self.metabolites:
            wb.create_sheet(k)
            wb[k]["A1"] = "Experiment"
            wb[k]["B1"] = "Fitted Isotopomers"
            wb[k]["C1"] = "Fitted Isotopomer %"
            wb[k]["D1"] = "Multiplets"
            wb[k]["E1"] = "Exp. %"
            wb[k]["F1"] = "Sim. %"
            wb[k]["G1"] = "Exp. GC-MS %"
            wb[k]["H1"] = "Sim. GC-MS %"
            wb[k]["I1"] = "Exp. NMR1D %"
            wb[k]["J1"] = "Sim. NMR1D %"
            index = 2
            for l in range(self.n_exps):
                wb[k][f'A{int(index)}'] = str(l + 1)
                if k in self.fitted_isotopomers.keys():
                    for m in range(len(self.fitted_isotopomers[k][l])):
                        wb[k][f'B{int(index + m)}'] = str(self.fitted_isotopomers[k][l][m])
                        wb[k][f'C{int(index + m)}'] = str(f'{self.fitted_isotopomer_percentages[k][l][m]:4.2f}')

                if k in self.fitted_multiplets.keys():
                    for m in range(len(self.fitted_multiplets[k][l])):
                        wb[k][f'D{int(index + m)}'] = str(self.fitted_multiplets[k][l][m])
                        wb[k][f'E{int(index + m)}'] = str(f'{self.exp_multiplet_percentages[k][l][m]:4.2f}')
                        wb[k][f'F{int(index + m)}'] = str(f'{self.fitted_multiplet_percentages[k][l][m]:4.2f}')

                if k in self.exp_gcms.keys():
                    wb[k][f'G{int(index)}'] = str(np.round(np.array(self.exp_gcms[k][l])*100)/100)
                    wb[k][f'H{int(index)}'] = str(np.round(np.array(self.fitted_gcms_percentages[k][l])*100)/100)
                if k in self.exp_nmr1d.keys():
                    wb[k][f'I{int(index)}'] = str(np.round(np.array(self.exp_nmr1d[k][l])*100)/100)
                    wb[k][f'J{int(index)}'] = str(np.round(np.array(self.fitted_nmr1d_percentages[k][l])*100)/100)
                index += int(np.max(np.array([len(self.fitted_isotopomers[k][l]), len(self.fitted_multiplets[k][l]), 1.0])))

            wb.save(file_name)

    # end export_data

    def fct_data(self, fit_parameters):
        self.chi2 = 0
        isos = self.fit_isotopomers[self.current_metabolite]
        isos.pop(0)
        self.set_fit_isotopomers_simple(metabolite=self.current_metabolite, isotopomers=isos, percentages=fit_parameters)
        if self.use_hsqc_multiplet_data:
            self.set_hsqc_isotopomers(self.current_metabolite)
            self.set_fct_hsqc_data(exp_index=self.current_experiment, metabolite=self.current_metabolite)
            self.chi2 += self.fct_hsqc_data(exp_index=self.current_experiment, metabolite=self.current_metabolite)

        if self.use_gcms_data:
            self.set_gcms_percentages(metabolite=self.current_metabolite)
            self.chi2 += self.fct_gcms_data(exp_index=self.current_experiment, metabolite=self.current_metabolite)

        if self.use_nmr1d_data:
            self.set_nmr1d_percentages(metabolite=self.current_metabolite)
            self.chi2 += self.fct_nmr1d_data(exp_index=self.current_experiment, metabolite=self.current_metabolite)

        return self.chi2
    # end fct_data

    def fit_data(self, exp_index=0, metabolite='', fit_isotopomers=[]):
        if len(metabolite) == 0 or len(fit_isotopomers) == 0:
            return

        use_hsqc = self.use_hsqc_multiplet_data
        use_gcms = self.use_gcms_data
        use_nmr1d = self.use_nmr1d_data
        self.current_experiment = exp_index
        self.current_metabolite = metabolite
        fit_parameters = np.ones(len(fit_isotopomers))
        fit_parameters /= fit_parameters.sum()
        fit_parameters *= 100
        fit_parameters = list(fit_parameters)
        self.set_fit_isotopomers_simple(metabolite=metabolite, isotopomers=fit_isotopomers, percentages=fit_parameters)
        eval_parameters = optimize.minimize(self.fct_data, fit_parameters, method='Powell')
        e_pars = np.array(eval_parameters.x).tolist()
        #self.fct_data(e_pars)
        self.set_fit_isotopomers_simple(metabolite=metabolite, isotopomers=fit_isotopomers, percentages=e_pars)
        self.fitted_isotopomers[metabolite][exp_index] = self.fit_isotopomers[metabolite]
        self.fitted_isotopomer_percentages[metabolite][exp_index] = self.isotopomer_percentages[metabolite]
        self.fitted_multiplets[metabolite][exp_index] = []
        self.fitted_multiplet_percentages[metabolite][exp_index] = []
        self.fitted_multiplets[metabolite][exp_index] = self.exp_multiplets[metabolite][exp_index]
        self.fitted_multiplet_percentages[metabolite][exp_index] = np.zeros(len(self.exp_multiplets[metabolite][exp_index]))
        for k in range(len(self.hsqc[metabolite])):
            if self.hsqc[metabolite][k] == 1:
                for l in range(len(self.hsqc_multiplets2[metabolite][exp_index][k])):
                    idx = self.exp_multiplets[metabolite][exp_index].index(self.hsqc_multiplets2[metabolite][exp_index][k][l])
                    self.fitted_multiplet_percentages[metabolite][exp_index][idx] = self.hsqc_multiplet_percentages[metabolite][exp_index][k][l]

        self.fitted_gcms_percentages[metabolite][exp_index] = self.gcms_percentages[metabolite]
        self.fitted_nmr1d_percentages[metabolite][exp_index] = self.nmr1d_percentages[metabolite]
        self.use_hsqc_multiplet_data = use_hsqc
        self.use_gcms_data = use_gcms
        self.use_nmr1d_data = use_nmr1d
        return
    # end fit_data

    def fit_all_exps(self, metabolite='', fit_isotopomers=[]):
        if len(metabolite) == 0 or len(fit_isotopomers) == 0:
            return

        print(f'Fitting all experiments for {metabolite}...')
        for k in range(self.n_exps):
            self.fit_data(exp_index=k, metabolite=metabolite, fit_isotopomers=fit_isotopomers)

    def set_fct_hsqc_data(self, exp_index=0, metabolite=''):
        if len(metabolite) == 0:
            return

        d = self.exp_multiplets[metabolite][exp_index]
        n = self.nmr_isotopomers[metabolite]
        p = self.nmr_isotopomer_percentages[metabolite]
        nn = []
        pp = []
        num_carbons = len(self.hsqc[metabolite])
        n_bonds = self.n_bonds[metabolite]
        for k in range(num_carbons):
            nn.append([])
            pp.append([])
            for l in range(len(n)):
                if n[l][k] == 1:
                    nn[k].append(n[l])
                    pp[k].append(p[l])

        for k in range(num_carbons):
            min_idx = max(0, k - n_bonds)
            max_idx = min(num_carbons, k + n_bonds + 1)
            for l in range(len(nn[k])):
                nnn = np.array(nn[k][l])
                for a in range(0, min_idx):
                    nnn[a] = 0

                for b in range(max_idx, num_carbons):
                    nnn[b] = 0

                nn[k][l] = list(nnn)

        mm = []
        qq = []
        for k in range(num_carbons):
            mm.append([])
            qq.append([])
            kk = -1
            ppp = list.copy(pp[k])
            nnn = list.copy(nn[k])
            while len(nnn) > 0:
                mm[k].append(nnn[0])
                qq[k].append(ppp[0])
                temp = list.copy(nnn[0])
                del nnn[0]
                del ppp[0]
                while temp in nnn:
                    kk += 1
                    idx = nnn.index(temp)
                    qq[k][kk] += ppp[idx]
                    del nnn[idx]
                    del ppp[idx]

        self.hsqc_multiplets[metabolite][exp_index] = mm
        for k in range(len(qq)):
            qq[k] = list(np.array(qq[k]) * 100.0 / np.sum(np.array(qq[k])))

        self.hsqc_multiplet_percentages[metabolite][exp_index] = qq
        mm2 = []
        for k in range(len(mm)):
            mm2.append([])
            for l in range(len(mm[k])):
                ll = list(np.where(np.array(mm[k][l]) == 1)[0] + 1)
                ll.pop(ll.index(k + 1))
                ll.insert(0, k+1)
                mm2[k].append(ll)

        self.hsqc_multiplets2[metabolite][exp_index] = mm2
        return
    # end set_fct_hsqc_data

    def set_sim_hsqc_data(self, exp_index=0, metabolite=''):
        if len(metabolite) == 0:
            return
        
        n = self.nmr_isotopomers[metabolite]
        p = self.nmr_isotopomer_percentages[metabolite]
        nn = []
        pp = []
        num_carbons = len(self.hsqc[metabolite])
        n_bonds = self.n_bonds[metabolite]
        for k in range(num_carbons):
            nn.append([])
            pp.append([])
            for l in range(len(n)):
                if n[l][k] == 1:
                    nn[k].append(n[l])
                    pp[k].append(p[l])
                    
        for k in range(num_carbons):
            min_idx = max(0, k - n_bonds)
            max_idx = min(num_carbons, k + n_bonds + 1)
            for l in range(len(nn[k])):
                nnn = np.array(nn[k][l])
                for a in range(0, min_idx):
                    nnn[a] = 0
                    
                for b in range(max_idx, num_carbons):
                    nnn[b] = 0
                    
                nn[k][l] = list(nnn)
                
        mm = []
        qq = []
        for k in range(num_carbons):
            mm.append([])
            qq.append([])
            kk = -1
            ppp = list.copy(pp[k])
            nnn = list.copy(nn[k])
            while len(nnn) > 0:
                mm[k].append(nnn[0])
                qq[k].append(ppp[0])
                temp = list.copy(nnn[0])
                del nnn[0]
                del ppp[0]
                while temp in nnn:
                    kk += 1
                    idx = nnn.index(temp)
                    qq[k][kk] += ppp[idx]
                    del nnn[idx]
                    del ppp[idx]
                    
        self.hsqc_multiplets[metabolite][exp_index] = mm
        for k in range(len(qq)):
            qq[k] = list(np.array(qq[k]) * 100.0 / np.sum(np.array(qq[k])))
            
        self.hsqc_multiplet_percentages[metabolite][exp_index] = qq
        mm2 = []
        for k in range(len(mm)):
            mm2.append([])
            for l in range(len(mm[k])):
                ll = list(np.where(np.array(mm[k][l]) == 1)[0] + 1)
                ll.pop(ll.index(k + 1))
                ll.insert(0, k+1)
                mm2[k].append(ll)
                
        self.hsqc_multiplets2[metabolite][exp_index] = mm2
        exp_multiplets = []
        exp_multiplet_percentages = []
        for k in range(len(self.hsqc[metabolite])):
            if self.hsqc[metabolite][k] == 1:
                for l in range(len(self.hsqc_multiplets2[metabolite][exp_index][k])):
                    exp_multiplets.append(self.hsqc_multiplets2[metabolite][exp_index][k][l])
                    exp_multiplet_percentages.append(self.hsqc_multiplet_percentages[metabolite][exp_index][k][l])

        self.exp_multiplets[metabolite][exp_index] = exp_multiplets
        self.exp_multiplet_percentages[metabolite][exp_index] = exp_multiplet_percentages
        return
    # end set_sim_hsqc_data
    
    def fct_hsqc_data(self, exp_index=0, metabolite=''):
        if len(metabolite) == 0:
            return -1

        pp = self.hsqc_multiplet_percentages[metabolite][exp_index]
        mm = self.hsqc_multiplets2[metabolite][exp_index]
        perc = list(np.array(self.exp_multiplet_percentages[metabolite][exp_index]))
        mult = self.exp_multiplets[metabolite][exp_index]
        for k in range(len(self.hsqc[metabolite])):
            if self.hsqc[metabolite][k] == 1:
                for l in range(len(mm[k])):
                    idx = mult.index(mm[k][l])
                    perc[idx] -= pp[k][l]

        perc = np.array(perc)
        perc *= perc
        chi2 = perc.sum()
        return chi2
    # end fct_hsqc_data

    def fct_gcms_data(self, exp_index=0, metabolite=''):
        if len(metabolite) == 0:
            return -1

        perc = np.array(self.gcms_percentages[metabolite]) - np.array(self.exp_gcms[metabolite][exp_index])
        perc *= perc
        chi2 = perc.sum()
        return chi2
    # end fct_gcms_data

    def fct_nmr1d_data(self, exp_index=0, metabolite=''):
        if len(metabolite) == 0:
            return -1

        hsqc = np.array(self.hsqc[metabolite])
        perc = np.array(self.exp_nmr1d[metabolite][exp_index])
        hsqc[np.where(perc < 0)[0]] = np.zeros(len(np.where(perc < 0)[0]))
        perc -= np.array(self.nmr1d_percentages[metabolite])
        perc *= hsqc
        perc *= perc
        chi2 = perc.sum()
        return chi2
    # end fct_nmr1d_data

    def init_metabolite(self, metabolite='', hsqc=[]):
        if len(metabolite) == 0 or len(hsqc) == 0:
            return
        
        self.__init__()
        self.metabolites.append(metabolite)
        self.fit_isotopomers[metabolite] = []
        self.isotopomer_percentages[metabolite] = []
        self.nmr_isotopomers[metabolite] = []
        self.nmr_isotopomer_percentages[metabolite] = []
        self.gcms_percentages[metabolite] = []
        self.nmr1d_percentages[metabolite] = []
        self.n_bonds[metabolite] = []
        self.hsqc_multiplets[metabolite] = []
        self.hsqc_multiplets2[metabolite] = []
        self.hsqc_multiplet_percentages[metabolite] = []
        self.fitted_isotopomers[metabolite] = []
        self.fitted_isotopomer_percentages[metabolite] = []
        self.fitted_multiplets[metabolite] = []
        self.fitted_multiplet_percentages[metabolite] = []
        self.fitted_gcms_percentages[metabolite] = []
        self.fitted_nmr1d_percentages[metabolite] = []
        self.exp_multiplets[metabolite] = []
        self.exp_multiplet_percentages[metabolite] = []
        self.exp_hsqc_isos[metabolite] = []
        self.n_bonds[metabolite] = 1
        self.hsqc_multiplets[metabolite] = []
        self.hsqc_multiplets2[metabolite] = []
        self.hsqc_multiplet_percentages[metabolite] = []
        self.exp_multiplets[metabolite].append([])
        self.exp_multiplet_percentages[metabolite].append([])
        self.exp_hsqc_isos[metabolite].append([])
        self.hsqc_multiplets[metabolite].append([])
        self.hsqc_multiplets2[metabolite].append([])
        self.hsqc_multiplet_percentages[metabolite].append([])
        self.fitted_isotopomers[metabolite].append([])
        self.fitted_isotopomer_percentages[metabolite].append([])
        self.fitted_multiplets[metabolite].append([])
        self.fitted_multiplet_percentages[metabolite].append([])
        self.fitted_gcms_percentages[metabolite].append([])
        self.fitted_nmr1d_percentages[metabolite].append([])
        self.hsqc[metabolite] = hsqc
        self.exp_gcms[metabolite] = []
        self.exp_gcms[metabolite].append([])
        self.n_exps = 1

        self.X_train = None
        self.X_val = None
        self.y_train = None
        self.y_val = None
        return

    def metabolite(self, metabolite=''):
        if len(metabolite) == 0 or metabolite not in self.metabolites:
            return

        r_string = '______________________________________________________________________________________\n'
        r_string += f'\nMetaboLabPy Isotopomer Data Analysis (v. {self.ver}) for {metabolite}\n'
        r_string += '______________________________________________________________________________________\n\n'
        r_string += 'HSQC NMR multiplet data:\n'
        print(r_string)
        print(self.nmr_multiplets[metabolite])
        print('\n\nGC-MS data:\n')
        print(self.gcms_data[metabolite])
        print('\n\nIsotopomer data:\n')
        print(f'isotopomers: {self.fit_isotopomers[metabolite]}\npercentages: {self.isotopomer_percentages[metabolite]}')
    # end metabolite

    def multiplet_fct(self):
        return
    # end multiplet_fct

    def read_hsqc_multiplets(self, file_name=''):
        if len(file_name) == 0:
            return

        self.nmr_multiplets = pd.read_excel(file_name, sheet_name=None, keep_default_na=False)
        #self.metabolites = []
        for k in self.nmr_multiplets.keys():
            self.metabolites.append(k)
            if k not in self.fit_isotopomers.keys():
                self.fit_isotopomers[k] = []
                self.isotopomer_percentages[k] = []
                self.nmr_isotopomers[k] = []
                self.nmr_isotopomer_percentages[k] = []
                self.gcms_percentages[k] = []
                self.nmr1d_percentages[k] = []
                self.n_bonds[k] = []
                self.hsqc_multiplets[k] = []
                self.hsqc_multiplets2[k] = []
                self.hsqc_multiplet_percentages[k] = []
                self.fitted_isotopomers[k] = []
                self.fitted_isotopomer_percentages[k] = []
                self.fitted_multiplets[k] = []
                self.fitted_multiplet_percentages[k] = []
                self.fitted_gcms_percentages[k] = []
                self.fitted_nmr1d_percentages[k] = []

        self.metabolites = sorted(list(set(self.metabolites)))
        self.n_exps = int(len(self.nmr_multiplets[self.metabolites[0]].keys())/6)
        for k in self.metabolites:
            self.hsqc[k] = list(map(int, self.nmr_multiplets[k]['HSQC.0'][0].split()))
            self.exp_multiplets[k] = []
            self.exp_multiplet_percentages[k] = []
            self.exp_hsqc_isos[k] = []
            self.n_bonds[k] = int(self.nmr_multiplets[k]['HSQC.0'][1].replace('n_bonds: ', ''))
            zero_iso = np.zeros(len(self.hsqc[k]))
            self.hsqc_multiplets[k] = []
            self.hsqc_multiplets2[k] = []
            self.hsqc_multiplet_percentages[k] = []
            for l in range(self.n_exps):
                multiplet_string  = f'Multiplet.{l}'
                percentages_string = f'Percentages.{l}'
                multiplets = self.nmr_multiplets[k][multiplet_string]
                percentages = self.nmr_multiplets[k][percentages_string]
                exp_multiplets = []
                exp_percentages = []
                exp_hsqc_isos = []
                for m in range(len(multiplets)):
                    exp_m = list(map(int, multiplets[m].replace(',', '').split()))
                    exp_multiplets.append(exp_m)
                    exp_percentages.append(percentages[m])

                self.exp_multiplets[k].append(exp_multiplets)
                self.exp_multiplet_percentages[k].append(exp_percentages)
                self.hsqc_multiplets[k].append([])
                self.hsqc_multiplets2[k].append([])
                self.hsqc_multiplet_percentages[k].append([])

            if len(self.fitted_isotopomers[k]) == 0:
                for l in range(self.n_exps):
                    self.fitted_isotopomers[k].append([])
                    self.fitted_isotopomer_percentages[k].append([])
                    self.fitted_multiplets[k].append([])
                    self.fitted_multiplet_percentages[k].append([])
                    self.fitted_gcms_percentages[k].append([])
                    self.fitted_nmr1d_percentages[k].append([])



        return
    # end read_hsqc_multiplets

    def read_nmr1d_data(self, file_name=''):
        if len(file_name) == 0:
            return

        self.nmr1d_data = pd.read_excel(file_name, sheet_name=None, keep_default_na=False)
        for k in self.nmr1d_data.keys():
            self.metabolites.append(k)
            if k not in self.fit_isotopomers.keys():
                self.fit_isotopomers[k] = []
                self.isotopomer_percentages[k] = []
                self.nmr_isotopomers[k] = []
                self.nmr_isotopomer_percentages[k] = []
                self.gcms_percentages[k] = []
                self.nmr1d_percentages[k] = []
                self.fitted_isotopomers[k] = []
                self.fitted_isotopomer_percentages[k] = []
                self.fitted_multiplets[k] = []
                self.fitted_multiplet_percentages[k] = []
                self.fitted_gcms_percentages[k] = []
                self.fitted_nmr1d_percentages[k] = []

            self.metabolites = sorted(list(set(self.metabolites)))
            self.n_exps = int(len(self.nmr1d_data[self.metabolites[0]].keys())/4)
            self.exp_nmr1d[k] = []
            for l in range(self.n_exps):
                percentages_string = f'Percentages.{l}'
                percentages = self.nmr1d_data[k][percentages_string]
                exp_percentages = []
                for m in range(len(percentages)):
                    exp_percentages.append(percentages[m])

                self.exp_nmr1d[k].append(exp_percentages)
                if len(self.fitted_isotopomers[k]) == 0:
                    self.fitted_isotopomers[k].append([])
                    self.fitted_isotopomer_percentages[k].append([])
                    self.fitted_multiplets[k].append([])
                    self.fitted_multiplet_percentages[k].append([])
                    self.fitted_gcms_percentages[k].append([])
                    self.fitted_nmr1d_percentages[k].append([])

        return
    # end read_nmr1d_data

    def read_gcms_data(self, file_name=''):
        if len(file_name) == 0:
            return

        self.gcms_data = pd.read_excel(file_name, sheet_name=None, keep_default_na=False)
        for k in self.gcms_data.keys():
            self.metabolites.append(k)
            if k not in self.fit_isotopomers.keys():
                self.fit_isotopomers[k] = []
                self.isotopomer_percentages[k] = []
                self.nmr_isotopomers[k] = []
                self.nmr_isotopomer_percentages[k] = []
                self.gcms_percentages[k] = []
                self.nmr1d_percentages[k] = []
                self.fitted_isotopomers[k] = []
                self.fitted_isotopomer_percentages[k] = []
                self.fitted_multiplets[k] = []
                self.fitted_multiplet_percentages[k] = []
                self.fitted_gcms_percentages[k] = []
                self.fitted_nmr1d_percentages[k] = []


            self.metabolites = sorted(list(set(self.metabolites)))
            self.n_exps = int(len(self.gcms_data[self.metabolites[0]].keys())/4)
            self.exp_gcms[k] = []
            for l in range(self.n_exps):
                percentages_string = f'Percentages.{l}'
                percentages = self.gcms_data[k][percentages_string]
                exp_percentages = []
                for m in range(len(percentages)):
                    exp_percentages.append(percentages[m])

                self.exp_gcms[k].append(exp_percentages)
                if len(self.fitted_isotopomers[k]) == 0:
                    self.fitted_isotopomers[k].append([])
                    self.fitted_isotopomer_percentages[k].append([])
                    self.fitted_multiplets[k].append([])
                    self.fitted_multiplet_percentages[k].append([])
                    self.fitted_gcms_percentages[k].append([])
                    self.fitted_nmr1d_percentages[k].append([])

        return
    # end read_gcms_data

    def reset_all_fit_isotopomers(self):
        for k in self.metabolites:
            self.reset_fit_isotopomers(k)
    # end reset_all_fit_isotopomers

    def reset_fit_isotopomers(self, metabolite=''):
        if len(metabolite) == 0 or metabolite not in self.metabolites:
            return

        self.fit_isotopomers[metabolite] = []
        self.isotopomer_percentages[metabolite] = []

    #end reset_fit_isotopomers

    def set_fit_isotopomers(self, metabolite='', isotopomers=[], percentages=[], exp_index=0):
        if len(metabolite) == 0 or metabolite not in self.metabolites or len(percentages) == 0 or len(isotopomers) == 0:
            print(
                'Usage:\nself.set_fit_isotopomers(metabolite="L-LacticAcid", isotopomers=[[0, 0, 1], [0, 1, 1]], percentages=[3, 5]')
            return

        if len(isotopomers) != len(percentages):
            print('Length of percentages vector does not match number of isotopomers')
            return

        self.fit_isotopomers[metabolite] = isotopomers
        self.isotopomer_percentages[metabolite] = percentages

        zero_isotopomer = list(np.zeros(len(self.fit_isotopomers[metabolite][0]), dtype=int))
        if zero_isotopomer not in self.fit_isotopomers[metabolite]:
            self.fit_isotopomers[metabolite].append(zero_isotopomer)
            self.isotopomer_percentages[metabolite].append(0.0)

        p_sum = sum(self.isotopomer_percentages[metabolite])
        idx = self.fit_isotopomers[metabolite].index(zero_isotopomer)
        if p_sum < 100.0:
            self.isotopomer_percentages[metabolite][idx] = 100.0 - p_sum + self.isotopomer_percentages[metabolite][idx]

        p_sum = sum(self.isotopomer_percentages[metabolite])
        for k in range(len(self.isotopomer_percentages[metabolite])):
            self.isotopomer_percentages[metabolite][k] *= 100.0 / p_sum

        new_isotopomer_list = []
        new_percentages_list = []
        new_isotopomer_list.append(self.fit_isotopomers[metabolite][idx])
        new_percentages_list.append(self.isotopomer_percentages[metabolite][idx])
        self.fit_isotopomers[metabolite].pop(idx)
        self.isotopomer_percentages[metabolite].pop(idx)
        while len(self.fit_isotopomers[metabolite]) > 0:
            new_isotopomer_list.append(self.fit_isotopomers[metabolite].pop(0))
            new_percentages_list.append(self.isotopomer_percentages[metabolite].pop(0))

        self.fit_isotopomers[metabolite] = new_isotopomer_list.copy()
        self.isotopomer_percentages[metabolite] = new_percentages_list.copy()

    # end set_fit_isotopomers

    def set_hsqc_isotopomers(self, metabolite=''):
        if len(metabolite) == 0 or metabolite not in self.metabolites:
            return

        self.nmr_isotopomers[metabolite] = []
        self.nmr_isotopomer_percentages[metabolite] = []
        for k in range(len(self.fit_isotopomers[metabolite])):
            n_zeros = len(self.fit_isotopomers[metabolite][k]) - sum(self.fit_isotopomers[metabolite][k])
            self.nmr_isotopomers[metabolite].append(self.fit_isotopomers[metabolite][k].copy())

            # Ensure percentage is correctly extracted
            percentage = self.isotopomer_percentages[metabolite][0][k] if isinstance(
                self.isotopomer_percentages[metabolite][0], list) else self.isotopomer_percentages[metabolite][k]


            pp = percentage * (1.0 - n_zeros * self.nat_abundance / 100.0)
            self.nmr_isotopomer_percentages[metabolite].append(pp)
            idx1 = 0
            for l in range(n_zeros):
                d2 = self.fit_isotopomers[metabolite][k].copy()
                idx2 = d2.index(0, idx1)
                d2[idx2] = 1
                idx1 = idx2 + 1
                self.nmr_isotopomers[metabolite].append(d2)
                pp = percentage * self.nat_abundance / 100.0
                self.nmr_isotopomer_percentages[metabolite].append(pp)

        new_nmr_isotopomers = []
        new_isotopomer_percentages = []
        for k in range(len(self.nmr_isotopomers[metabolite])):
            if self.nmr_isotopomers[metabolite][k] not in new_nmr_isotopomers:
                new_nmr_isotopomers.append(self.nmr_isotopomers[metabolite][k].copy())
                new_isotopomer_percentages.append(self.nmr_isotopomer_percentages[metabolite][k])
            else:
                idx = new_nmr_isotopomers.index(self.nmr_isotopomers[metabolite][k])
                new_isotopomer_percentages[idx] += self.nmr_isotopomer_percentages[metabolite][k]

        self.nmr_isotopomers[metabolite] = new_nmr_isotopomers.copy()
        self.nmr_isotopomer_percentages[metabolite] = new_isotopomer_percentages.copy()

    # end set_hsqc_isotopomers

    def set_gcms_percentages(self, metabolite=''):
        if len(metabolite) == 0 or metabolite not in self.metabolites:
            return

        d_sums = []
        for k in range(len(self.fit_isotopomers[metabolite])):
            d_sums.append(sum(self.fit_isotopomers[metabolite][k]))

        d_sums = np.array(d_sums)
        if metabolite in self.nmr_multiplets.keys():
            gcms_data = list(np.zeros(len(self.nmr_multiplets[metabolite]['HSQC.0'][0].split()) + 1, dtype=int))
        elif metabolite in self.gcms_data.keys():
            gcms_data = list(np.zeros(len(self.gcms_data[metabolite]), dtype=int))
        else:
            gcms_data = list(np.zeros(len(self.hsqc[metabolite]) + 1, dtype=int))

        # introduce flattened percentages
        flattened_percentages = [item for sublist in self.isotopomer_percentages[metabolite] for item in sublist]
        percentages = np.array(flattened_percentages)

        for k in range(len(gcms_data)):
            gcms_data[k] = percentages[np.where(d_sums == k)].sum()

        self.gcms_percentages[metabolite] = gcms_data.copy()
    # end set_gcms_isotopomers

    def sim_hsqc_data(self, metabolite='', exp_index=0, isotopomers=[], percentages=[]):
        if len(metabolite) == 0 or metabolite not in self.metabolites or len(isotopomers) == 0 or len(percentages) == 0:
            return

        self.set_fit_isotopomers_simple(metabolite, isotopomers, percentages)
        self.set_hsqc_isotopomers(metabolite)
        # sim effect of read_hsqc_data
        self.set_sim_hsqc_data(exp_index, metabolite)
        return
        # end sim_hsqc_data

    def sim_gcms_data(self, metabolite='', exp_index=0):
        if len(metabolite) == 0 or metabolite not in self.metabolites:
            return

        self.set_gcms_percentages(metabolite)
        self.exp_gcms[metabolite][exp_index] = list(np.copy(self.gcms_percentages[metabolite]))
    # end sim_gcms_data

    def set_nmr1d_percentages(self, metabolite=''):
        if len(metabolite) == 0 or metabolite not in self.metabolites:
            return

        self.nmr1d_percentages[metabolite] = np.zeros(len(self.nmr_isotopomers[metabolite][0]))
        for k in range(len(self.nmr_isotopomers[metabolite])):
            self.nmr1d_percentages[metabolite] += np.array(self.nmr_isotopomers[metabolite][k])*self.nmr_isotopomer_percentages[metabolite][k]

        self.nmr1d_percentages[metabolite] = list(self.nmr1d_percentages[metabolite])
    # end set_nmr1d_isotopomers

# Neural Network Addition - RA

    def generate_isotopomer_percentages(self, metabolite_type):
        if metabolite_type == 'three-carbon':
            num_isotopomers = 7
        elif metabolite_type == 'four_carbon':
            num_isotopomers = 15
        elif metabolite_type == 'five_carbon':
            num_isotopomers = 31
        elif metabolite_type == 'six_carbon':
            num_isotopomers = 63
        else:
            raise ValueError("Unknown metabolite type")

        unlabelled_percentage = np.random.uniform(20, 80)
        remaining_percentage = 100 - unlabelled_percentage

        isotopomer_presence = np.random.rand(num_isotopomers) < 0.5
        present_isotopomers = np.sum(isotopomer_presence)

        if present_isotopomers == 0:
            isotopomer_presence[np.random.randint(0, num_isotopomers)] = True
            present_isotopomers = 1

        random_values = np.random.rand(present_isotopomers)
        random_percentages = (random_values / random_values.sum()) * remaining_percentage

        isotopomer_percentages = [0] * num_isotopomers
        random_idx = 0
        for i in range(num_isotopomers):
            if isotopomer_presence[i]:
                isotopomer_percentages[i] = random_percentages[random_idx]
                random_idx += 1

        percentages = [unlabelled_percentage] + isotopomer_percentages

        return percentages

    def add_noise(self, data, noise_level):
        data = np.array(data)
        noise = np.random.normal(0, noise_level, data.shape)
        return data + noise

    def add_noise_to_hsqc_gcms(self, metabolite, num_samples, hsqc_noise_level, gcms_noise_level):
        for exp_index in range(num_samples):
            self.exp_multiplet_percentages[metabolite][exp_index] = self.add_noise(
                self.exp_multiplet_percentages[metabolite][exp_index], hsqc_noise_level
            ).tolist()
            self.exp_gcms[metabolite][exp_index] = self.add_noise(
                self.exp_gcms[metabolite][exp_index], gcms_noise_level
            ).tolist()

    def get_training_data(self, metabolite='', exp_index=0):
        hsqc_data = self.exp_multiplet_percentages[metabolite][exp_index]
        gcms_data = self.exp_gcms[metabolite][exp_index]

        flattened_gcms = np.array(gcms_data).flatten()
        features = np.hstack([hsqc_data, flattened_gcms])
        return np.array(features).reshape(1, -1)

    def get_training_labels(self, metabolite='', exp_index=0, percentages=[], fit_isotopomers=None):
        labels = percentages[exp_index]
        return np.array(labels).reshape(1, -1)

    def init_metabolite_multiple_samples(self, metabolites=[], hsqc=[], num_samples=1000):
        for metabolite in metabolites:
            if len(metabolite) == 0 or len(hsqc) == 0:
                continue

            if metabolite not in self.metabolites:
                self.metabolites.append(metabolite)

            self.fit_isotopomers[metabolite] = []
            self.isotopomer_percentages[metabolite] = [[] for _ in range(num_samples)]
            self.nmr_isotopomers[metabolite] = []
            self.nmr_isotopomer_percentages[metabolite] = []
            self.gcms_percentages[metabolite] = []
            self.nmr1d_percentages[metabolite] = []
            self.n_bonds[metabolite] = 1
            self.hsqc[metabolite] = hsqc

            self.exp_multiplets[metabolite] = [[] for _ in range(num_samples)]
            self.exp_multiplet_percentages[metabolite] = [[] for _ in range(num_samples)]
            self.exp_gcms[metabolite] = [[] for _ in range(num_samples)]
            self.exp_hsqc_isos[metabolite] = [[] for _ in range(num_samples)]
            self.hsqc_multiplets[metabolite] = [[] for _ in range(num_samples)]
            self.hsqc_multiplets2[metabolite] = [[] for _ in range(num_samples)]
            self.hsqc_multiplet_percentages[metabolite] = [[] for _ in range(num_samples)]
            self.fitted_isotopomers[metabolite] = [[] for _ in range(num_samples)]
            self.fitted_isotopomer_percentages[metabolite] = [[] for _ in range(num_samples)]
            self.fitted_multiplets[metabolite] = [[] for _ in range(num_samples)]
            self.fitted_multiplet_percentages[metabolite] = [[] for _ in range(num_samples)]
            self.fitted_gcms_percentages[metabolite] = [[] for _ in range(num_samples)]
            self.fitted_nmr1d_percentages[metabolite] = [[] for _ in range(num_samples)]

            self.n_exps = num_samples

    def set_fit_isotopomers_simple(self, metabolite='', isotopomers=[], percentages=[], exp_index=0):
        if len(metabolite) == 0 or metabolite not in self.metabolites or len(percentages) == 0 or len(isotopomers) == 0:
            print(
                'Usage:\nself.set_fit_isotopomers_simple(metabolite="L-LacticAcid", isotopomers=[[0, 0, 1], [0, 1, 1]], percentages=[3, 5]')
            return

        if len(isotopomers) != len(percentages):
            print('Length of percentages vector does not match number of isotopomers')
            return

        while len(self.isotopomer_percentages[metabolite]) <= exp_index:
            self.isotopomer_percentages[metabolite].append([])

        self.isotopomer_percentages[metabolite][exp_index] = percentages
        self.fit_isotopomers[metabolite] = isotopomers

    def create_nn_model(self, hp, input_dim, output_dim):
        l2_lambda = hp.Float('l2_lambda', min_value=0.001, max_value=0.1, step=0.001)
        learning_rate = hp.Float('learning_rate', min_value=0.001, max_value=0.1, step=0.001)
        num_neurons = hp.Int('num_neurons', min_value=64, max_value=256, step=32)
        num_layers = hp.Int('num_layers', min_value=2, max_value=4)
        dropout_rate = hp.Float('dropout_rate', min_value=0.0, max_value=0.5, step=0.1)

        model = Sequential()
        model.add(Dense(num_neurons, activation='relu', kernel_regularizer=l2(l2_lambda), input_dim=input_dim))
        for _ in range(num_layers - 1):
            model.add(Dense(num_neurons, activation='relu', kernel_regularizer=l2(l2_lambda)))
            if dropout_rate > 0:
                model.add(Dropout(dropout_rate))
        model.add(Dense(output_dim, activation='relu'))
        model.add(Lambda(lambda x: 100 * x / tf.reduce_sum(x, axis=1, keepdims=True)))

        optimizer = Adam(learning_rate=learning_rate)
        model.compile(optimizer=optimizer, loss='mean_squared_error')
        return model

    def fit_data_nn(self, metabolite='', fit_isotopomers=None, num_samples=1000, percentages=[], hsqc=None,
                    tuner_project_name=''):
        if len(metabolite) == 0 or fit_isotopomers is None:
            print("Metabolite or isotopomers are not provided")
            return

        use_hsqc = self.use_hsqc_multiplet_data
        use_gcms = self.use_gcms_data
        use_nmr1d = self.use_nmr1d_data

        self.use_hsqc_multiplet_data = True
        self.use_gcms_data = True
        self.use_nmr1d_data = False

        X = []
        y = []

        for exp_index in range(num_samples):
            self.current_experiment = exp_index
            self.current_metabolite = metabolite

            X_train = self.get_training_data(metabolite=metabolite, exp_index=exp_index)
            y_train = self.get_training_labels(metabolite=metabolite, exp_index=exp_index, percentages=percentages)

            if X_train is None or y_train is None:
                print(f"Training data or labels could not be retrieved for sample {exp_index}.")
                continue

            X.append(X_train)
            y.append(y_train)

        X = np.array(X).reshape(num_samples, -1)
        y = np.array(y).reshape(num_samples, -1)

        X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.1, random_state=42)

        input_dim = X_train.shape[1]
        output_dim = y_train.shape[1]

        tuner_dir = f'my_dir/{tuner_project_name}'
        if os.path.exists(tuner_dir):
            shutil.rmtree(tuner_dir)

        tuner = kt.RandomSearch(
            lambda hp: self.create_nn_model(hp, input_dim, output_dim),
            objective='val_loss',
            max_trials=200,
            executions_per_trial=1,
            directory='my_dir',
            project_name=tuner_project_name
        )

        early_stopping = EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)

        tuner.search(X_train, y_train, epochs=100, validation_data=(X_val, y_val), callbacks=[early_stopping])

        best_hps = tuner.get_best_hyperparameters(num_trials=1)[0]

        model = self.create_nn_model(best_hps, input_dim, output_dim)
        model.fit(X_train, y_train, epochs=100, batch_size=32, validation_data=(X_val, y_val),
                  callbacks=[early_stopping])

        y_pred = model.predict(X_val)

        mae_metric = MeanAbsoluteError()
        mae_metric.update_state(y_val, y_pred)
        mae = mae_metric.result().numpy()

        self.results.append({
            'Metabolite': metabolite,
            'HSQC': hsqc,
            'Best Hyperparameters': best_hps.values,
            'Best Validation Loss': tuner.oracle.get_best_trials(num_trials=1)[0].score,
            'Mean Absolute Percentage Error': mae
        })

        # Access the isotopomer distribution for the test sample (Last sample in validation set)
        test_sample_index = 99  # Index of the test sample in validation set
        test_sample_pred = y_pred[test_sample_index]
        test_sample_true = y_val[test_sample_index]

        print(f'Test Sample True Isotopomer Distribution: {test_sample_true}')
        print(f'Test Sample Predicted Isotopomer Distribution: {test_sample_pred}')

        self.use_hsqc_multiplet_data = use_hsqc
        self.use_gcms_data = use_gcms
        self.use_nmr1d_data = use_nmr1d

    def save_results(self, filename='results.xlsx'):
        df = pd.DataFrame(self.results)
        df.to_excel(filename, index=False)