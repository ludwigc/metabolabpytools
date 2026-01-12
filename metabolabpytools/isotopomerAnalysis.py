import os
import pandas as pd
import numpy as np
import scipy as sp
import os
import math
import time
from scipy import optimize
from openpyxl import Workbook  # pragma: no cover
import tensorflow as tf
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import layers, models
from keras_tuner import BayesianOptimization
import numpy as np
import pandas as pd
import os
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from collections import defaultdict  # Add this import statement




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
                    wb[k][f'G{int(index)}'] = str(np.round(np.array(self.exp_gcms[k][l]) * 100) / 100)
                    wb[k][f'H{int(index)}'] = str(np.round(np.array(self.fitted_gcms_percentages[k][l]) * 100) / 100)
                if k in self.exp_nmr1d.keys():
                    wb[k][f'I{int(index)}'] = str(np.round(np.array(self.exp_nmr1d[k][l]) * 100) / 100)
                    wb[k][f'J{int(index)}'] = str(np.round(np.array(self.fitted_nmr1d_percentages[k][l]) * 100) / 100)
                index += int(
                    np.max(np.array([len(self.fitted_isotopomers[k][l]), len(self.fitted_multiplets[k][l]), 1.0])))

            wb.save(file_name)

    # end export_data

    def fct_data(self, fit_parameters):
        self.chi2 = 0
        isos = self.fit_isotopomers[self.current_metabolite]
        isos.pop(0)
        self.set_fit_isotopomers(metabolite=self.current_metabolite, isotopomers=isos, percentages=fit_parameters)
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
        self.set_fit_isotopomers(metabolite=metabolite, isotopomers=fit_isotopomers, percentages=fit_parameters)
        eval_parameters = optimize.minimize(self.fct_data, fit_parameters, method='Powell')
        e_pars = np.array(eval_parameters.x).tolist()
        # self.fct_data(e_pars)
        self.set_fit_isotopomers(metabolite=metabolite, isotopomers=fit_isotopomers, percentages=e_pars)
        self.fitted_isotopomers[metabolite][exp_index] = self.fit_isotopomers[metabolite]
        self.fitted_isotopomer_percentages[metabolite][exp_index] = self.isotopomer_percentages[metabolite]
        self.fitted_multiplets[metabolite][exp_index] = []
        self.fitted_multiplet_percentages[metabolite][exp_index] = []
        self.fitted_multiplets[metabolite][exp_index] = self.exp_multiplets[metabolite][exp_index]
        self.fitted_multiplet_percentages[metabolite][exp_index] = np.zeros(
            len(self.exp_multiplets[metabolite][exp_index]))
        for k in range(len(self.hsqc[metabolite])):
            if self.hsqc[metabolite][k] == 1:
                for l in range(len(self.hsqc_multiplets2[metabolite][exp_index][k])):
                    idx = self.exp_multiplets[metabolite][exp_index].index(
                        self.hsqc_multiplets2[metabolite][exp_index][k][l])
                    self.fitted_multiplet_percentages[metabolite][exp_index][idx] = \
                    self.hsqc_multiplet_percentages[metabolite][exp_index][k][l]

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
                    if kk >= len(qq[k]):
                        qq[k].append(0)  # Ensure qq[k] is long enough
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
                ll.insert(0, k + 1)
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

                nn[k][l] = nnn.tolist()

        nnn = []
        ppp = []
        for k in range(len(nn)):
            nnn.append([])
            ppp.append([])
            for l in range(len(nn[k])):
                if nn[k][l] in nnn[k]:
                    ppp[k][nnn[k].index(nn[k][l])] += pp[k][l]
                else:
                    nnn[k].append(nn[k][l])
                    ppp[k].append(pp[k][l])

        nn = nnn
        pp = ppp

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
                    # Ensure that qq[k] has enough elements
                    if kk >= len(qq[k]):
                        qq[k].append(0)  # Append a zero if the list is not long enough
                    qq[k][kk] += ppp[idx]
                    del nnn[idx]
                    del ppp[idx]

        self.hsqc_multiplets[metabolite][exp_index] = mm
        for k in range(len(qq)):
            qq[k] = (np.array(qq[k]) * 100.0 / np.sum(np.array(qq[k]))).tolist()

        self.hsqc_multiplet_percentages[metabolite][exp_index] = qq
        mm2 = []
        for k in range(len(mm)):
            mm2.append([])
            for l in range(len(mm[k])):
                ll = (np.where(np.array(mm[k][l]) == 1)[0] + 1).tolist()
                ll.pop(ll.index(k + 1))
                ll.insert(0, k + 1)
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
        print(
            f'isotopomers: {self.fit_isotopomers[metabolite]}\npercentages: {self.isotopomer_percentages[metabolite]}')

    # end metabolite

    def multiplet_fct(self):
        return

    # end multiplet_fct

    def read_hsqc_multiplets(self, file_name=''):
        if len(file_name) == 0:
            return

        self.nmr_multiplets = pd.read_excel(file_name, sheet_name=None, keep_default_na=False)
        # self.metabolites = []
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
        self.n_exps = int(len(self.nmr_multiplets[self.metabolites[0]].keys()) / 6)
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
                multiplet_string = f'Multiplet.{l}'
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
            self.n_exps = int(len(self.nmr1d_data[self.metabolites[0]].keys()) / 4)
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
            self.n_exps = int(len(self.gcms_data[self.metabolites[0]].keys()) / 4)
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

    # end reset_fit_isotopomers

    def set_fit_isotopomers(self, metabolite='', isotopomers=[], percentages=[]):
        if len(metabolite) == 0 or metabolite not in self.metabolites or len(isotopomers) == 0 or len(percentages) == 0:
            print(
                'Usage:\nself.set_fit_isotopomers(metabolite="L-LacticAcid", isotopomers=[[0, 0, 1], [0, 1, 1]], percentages=[3, 5]')
            return

        if len(isotopomers) != len(percentages):
            print('length of percentages vector does not match number of isotopomers')
            return

        self.reset_fit_isotopomers(metabolite)
        for k in range(len(isotopomers)):
            self.fit_isotopomers[metabolite].append(isotopomers[k])
            self.isotopomer_percentages[metabolite].append(percentages[k])

        zero_isotopomer = np.zeros(len(self.fit_isotopomers[metabolite][0]), dtype=int).tolist()
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
            pp = self.isotopomer_percentages[metabolite][k] * (1.0 - n_zeros * self.nat_abundance / 100.0)
            self.nmr_isotopomer_percentages[metabolite].append(pp)
            idx1 = 0
            for l in range(n_zeros):
                d2 = self.fit_isotopomers[metabolite][k].copy()
                idx2 = d2.index(0, idx1)
                d2[idx2] = 1
                idx1 = idx2 + 1
                self.nmr_isotopomers[metabolite].append(d2)
                pp = self.isotopomer_percentages[metabolite][k] * self.nat_abundance / 100.0
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
            gcms_data = np.zeros(len(self.nmr_multiplets[metabolite]['HSQC.0'][0].split()) + 1, dtype=int).tolist()
        elif metabolite in self.gcms_data.keys():
            gcms_data = np.zeros(len(self.gcms_data[metabolite]), dtype=int).tolist()
        else:
            gcms_data = np.zeros(len(self.hsqc[metabolite]) + 1, dtype=int).tolist()

        percentages = np.array(self.isotopomer_percentages[metabolite].copy())
        for k in range(len(gcms_data)):
            gcms_data[k] = percentages[np.where(d_sums == k)].sum()

        self.gcms_percentages[metabolite] = gcms_data.copy()

    # end set_gcms_isotopomers

    def sim_hsqc_data(self, metabolite='', exp_index=0, isotopomers=[], percentages=[]):
        if len(metabolite) == 0 or metabolite not in self.metabolites or len(isotopomers) == 0 or len(percentages) == 0:
            return

        self.set_fit_isotopomers(metabolite, isotopomers, percentages)
        self.set_hsqc_isotopomers(metabolite)
        # sim effect of read_hsqc_data
        self.set_sim_hsqc_data(exp_index, metabolite)
        return
        # end sim_hsqc_data

    def sim_gcms_data(self, metabolite='', exp_index=0):
        if len(metabolite) == 0 or metabolite not in self.metabolites:
            return

        self.set_gcms_percentages(metabolite)
        self.exp_gcms[metabolite][exp_index] = np.copy(self.gcms_percentages[metabolite]).tolist()

    # end sim_gcms_data

    def set_nmr1d_percentages(self, metabolite=''):
        if len(metabolite) == 0 or metabolite not in self.metabolites:
            return

        self.nmr1d_percentages[metabolite] = np.zeros(len(self.nmr_isotopomers[metabolite][0]))
        for k in range(len(self.nmr_isotopomers[metabolite])):
            self.nmr1d_percentages[metabolite] += np.array(self.nmr_isotopomers[metabolite][k]) * \
                                                  self.nmr_isotopomer_percentages[metabolite][k]

        self.nmr1d_percentages[metabolite] = self.nmr1d_percentages[metabolite].tolist()
        # end set_nmr1d_isotopomers

    pass


class IsotopomerAnalysisNN(IsotopomerAnalysis):

    def __init__(self):
        super().__init__()

    def simulate_hsqc_gcms(self, distributions, hsqc_vector):
        all_isotopomer_data = []
        all_hsqc_data = []
        all_gcms_data = []

        for i, distribution in enumerate(distributions):
            isotopomers = distribution['isotopomers']
            percentages = distribution['percentages']

            # Use existing logic to simulate HSQC and GC-MS data
            metabolite_name = f'TestMetabolite_{i + 1}'
            self.init_metabolite(metabolite_name, hsqc=hsqc_vector)  # HSQC vector defined externally

            # Simulate the isotopomer data
            self.sim_hsqc_data(metabolite=metabolite_name, exp_index=0, isotopomers=isotopomers,
                               percentages=percentages)
            self.sim_gcms_data(metabolite=metabolite_name, exp_index=0)

            # Retrieve simulated data
            hsqc_multiplets = self.exp_multiplets[metabolite_name][0]
            hsqc_percentages = self.exp_multiplet_percentages[metabolite_name][0]
            gcms_percentages = self.exp_gcms[metabolite_name][0]

            # Create separate DataFrames for each sample
            isotopomer_df = pd.DataFrame({
                'isotopomer': [str(iso) for iso in isotopomers],
                'isotopomer%': percentages,
                'sample': i + 1
            })

            hsqc_df = pd.DataFrame({
                'hsqc_multiplet': [str(hsqc) for hsqc in hsqc_multiplets],
                'hsqc%': hsqc_percentages,
                'sample': i + 1
            })

            gcms_df = pd.DataFrame({
                'GCMS%': gcms_percentages,
                'sample': i + 1
            })

            # Append to the list of all data
            all_isotopomer_data.append(isotopomer_df)
            all_hsqc_data.append(hsqc_df)
            all_gcms_data.append(gcms_df)

        # Concatenate all DataFrames into a single DataFrame for each type
        combined_isotopomer_data = pd.concat(all_isotopomer_data, ignore_index=True)
        combined_hsqc_data = pd.concat(all_hsqc_data, ignore_index=True)
        combined_gcms_data = pd.concat(all_gcms_data, ignore_index=True)

        return combined_isotopomer_data, combined_hsqc_data, combined_gcms_data

    def save_simulation_data(self, combined_isotopomer_data, combined_hsqc_data, combined_gcms_data, hsqc_vector):
        # Ensure the directory exists
        os.makedirs('sim_data', exist_ok=True)

        # Convert the HSQC vector to a string for the file name
        hsqc_vector_str = ''.join(map(str, hsqc_vector))
        file_name = f'sim_data/sim_{hsqc_vector_str}.xlsx'

        # Save to Excel with multiple sheets
        with pd.ExcelWriter(file_name) as writer:
            combined_isotopomer_data.to_excel(writer, sheet_name='Isotopomer_Data', index=False)
            combined_hsqc_data.to_excel(writer, sheet_name='HSQC_Data', index=False)
            combined_gcms_data.to_excel(writer, sheet_name='GCMS_Data', index=False)

        print(f"Data successfully saved to {file_name}")



    def generate_isotopomer_distributions(self, n_distributions=10000, n_carbons=3):
        distributions = []
        for _ in range(n_distributions):
            isotopomers = []
            percentages = []

            # Generate all possible isotopomer patterns for n_carbons
            possible_isotopomers = [map(int, bin(i)[2:].zfill(n_carbons)).tolist() for i in range(2 ** n_carbons)]

            # Coin flip to determine inclusion for each isotopomer, always include [0, 0, 0]
            for isotopomer in possible_isotopomers:
                if np.random.rand() > 0.5 or isotopomer == [0] * n_carbons:
                    isotopomers.append(isotopomer)

            # Assign the unlabeled isotopomer percentage within 10-30%
            unlabeled_percentage = np.random.uniform(30, 100)
            percentages.append(unlabeled_percentage)

            # Calculate the remaining percentage
            remaining_percentage = 100 - unlabeled_percentage

            # Generate random proportions for the remaining isotopomers
            num_other_isotopomers = len(isotopomers) - 1  # Excluding the unlabeled one
            if num_other_isotopomers > 0:
                random_proportions = np.random.rand(num_other_isotopomers)
                random_proportions /= random_proportions.sum()  # Normalize to sum to 1

                # Assign the remaining percentage based on the random proportions
                for proportion in random_proportions:
                    percentages.append(proportion * remaining_percentage)

            # Final check to ensure percentages sum exactly to 100
            percentages = [p * (100 / sum(percentages)) for p in percentages]

            distributions.append({'isotopomers': isotopomers, 'percentages': percentages})

        return distributions


    def load_spreadsheet_by_hsqc_vector(self, hsqc_vector):
        # Convert HSQC vector to string
        hsqc_vector_str = ''.join(map(str, hsqc_vector))
        file_name = f'sim_data/sim_{hsqc_vector_str}.xlsx'

        if os.path.exists(file_name):
            with pd.ExcelFile(file_name) as xls:
                isotopomer_data = pd.read_excel(xls, sheet_name='Isotopomer_Data')
                hsqc_data = pd.read_excel(xls, sheet_name='HSQC_Data')
                gcms_data = pd.read_excel(xls, sheet_name='GCMS_Data')
            return isotopomer_data, hsqc_data, gcms_data
        else:
            raise FileNotFoundError(f"No spreadsheet found for HSQC vector {hsqc_vector}")

    def collate_y_labels(self, isotopomer_data, n_carbons):
        # Define the 8 possible isotopomers for a 3-carbon metabolite
        possible_isotopomers = [map(int, bin(i)[2:].zfill(n_carbons)).tolist() for i in range(2 ** n_carbons)]

        # Convert the possible isotopomers to strings to match the data format
        possible_isotopomers_str = [str(iso) for iso in possible_isotopomers]

        Y = []

        for sample in isotopomer_data['sample'].unique():
            # Initialize Y_sample with zeros
            Y_sample = np.zeros(len(possible_isotopomers_str))

            # Filter data for the current sample
            sample_data = isotopomer_data[isotopomer_data['sample'] == sample]

            # Populate Y_sample
            for i, iso_str in enumerate(possible_isotopomers_str):
                match = sample_data[sample_data['isotopomer'] == iso_str]
                if not match.empty:
                    Y_sample[i] = match['isotopomer%'].values[0]

            Y.append(Y_sample)

        return np.array(Y)

    def generate_possible_hsqc_multiplets(self, hsqc_vector):
        active_positions = [i + 1 for i, x in enumerate(hsqc_vector) if x == 1]
        possible_multiplets = []
        max_position = len(hsqc_vector)

        # Define possible multiplets for each position, ensuring connections do not exceed max_position
        all_multiplets = {
            1: [[1], [1, 2]] if max_position >= 2 else [[1]],
            2: [[2], [2, 1], [2, 3], [2, 1, 3]] if max_position >= 3 else [[2], [2, 1]],
            3: [[3], [3, 4], [3, 2], [3, 2, 4]] if max_position >= 4 else [[3], [3, 2]],
            4: [[4], [4, 5], [4, 3], [4, 3, 5]] if max_position >= 5 else [[4], [4, 3]],
            5: [[5], [5, 6], [5, 4], [5, 4, 6]] if max_position >= 6 else [[5], [5, 4]],
            6: [[6], [6, 5]] if max_position >= 6 else [[6]],
            # Extend if needed for more carbons
        }

        # Loop through active positions and gather their possible multiplets
        for pos in active_positions:
            if pos in all_multiplets:
                possible_multiplets.extend(all_multiplets[pos])

        return possible_multiplets

    def collate_x_labels_with_noise(self, hsqc_data, gcms_data, all_possible_hsqc_multiplets, hsqc_noise_level=0.025,
                                    gcms_noise_level=0.075):
        X = []

        for sample in hsqc_data['sample'].unique():
            # Filter data for the current sample and create copies to avoid SettingWithCopyWarning
            sample_hsqc_data = hsqc_data[hsqc_data['sample'] == sample].copy()
            sample_gcms_data = gcms_data[gcms_data['sample'] == sample].copy()

            # Add noise to HSQC data
            sample_hsqc_data.loc[:, 'hsqc%'] += np.random.normal(0, hsqc_noise_level, sample_hsqc_data['hsqc%'].shape)

            # Initialize a dictionary to hold hsqc percentages, filling with zeros initially
            hsqc_dict = {str(multiplet): 0 for multiplet in all_possible_hsqc_multiplets}

            # Fill in the actual hsqc percentages from the sample data
            for _, row in sample_hsqc_data.iterrows():
                hsqc_dict[str(row['hsqc_multiplet'])] = row['hsqc%']

            # Extract hsqc percentages in the order of all_possible_hsqc_multiplets
            hsqc_percentages_ordered = [hsqc_dict[str(multiplet)] for multiplet in all_possible_hsqc_multiplets]

            # Add noise to GC-MS data
            sample_gcms_data.loc[:, 'GCMS%'] += np.random.normal(0, gcms_noise_level, sample_gcms_data['GCMS%'].shape)

            # Combine HSQC multiplet percentages and GC-MS percentages
            X_sample = np.hstack([hsqc_percentages_ordered, sample_gcms_data['GCMS%'].values])

            X.append(X_sample)

        # Convert list to array
        return np.array(X)

    def collate_x_labels_without_noise(self, hsqc_data, gcms_data, all_possible_hsqc_multiplets):
        X = []

        for sample in hsqc_data['sample'].unique():
            # Filter data for the current sample and create copies to avoid SettingWithCopyWarning
            sample_hsqc_data = hsqc_data[hsqc_data['sample'] == sample].copy()
            sample_gcms_data = gcms_data[gcms_data['sample'] == sample].copy()

            # Initialize a dictionary to hold hsqc percentages, filling with zeros initially
            hsqc_dict = {str(multiplet): 0 for multiplet in all_possible_hsqc_multiplets}

            # Fill in the actual hsqc percentages from the sample data
            for _, row in sample_hsqc_data.iterrows():
                hsqc_dict[str(row['hsqc_multiplet'])] = row['hsqc%']

            # Extract hsqc percentages in the order of all_possible_hsqc_multiplets
            hsqc_percentages_ordered = [hsqc_dict[str(multiplet)] for multiplet in all_possible_hsqc_multiplets]

            # Combine HSQC multiplet percentages and GC-MS percentages without adding noise
            X_sample = np.hstack([hsqc_percentages_ordered, sample_gcms_data['GCMS%'].values])

            X.append(X_sample)

        # Convert list to array
        return np.array(X)

    def create_dynamic_nn_model(self, input_dim, output_dim):
        model = models.Sequential()
        model.add(layers.Input(shape=(input_dim,)))  # Define input shape using Input layer
        model.add(layers.Dense(128, activation='relu'))
        model.add(layers.Dense(64, activation='relu'))
        model.add(layers.Dense(32, activation='relu'))
        # Use ReLU activation in the output layer to ensure non-negative outputs
        model.add(layers.Dense(output_dim, activation='relu'))
        return model

    def train_neural_network(self, X_noisy, Y, plot_filename="trainingloss0110.png", epochs=100, batch_size=32, validation_split=0.2):
        input_dim = X_noisy.shape[1]  # Number of features in X
        output_dim = Y.shape[1]  # Number of isotopomers in Y

        # Create the model
        model = self.create_dynamic_nn_model(input_dim, output_dim)

        # Compile the model
        model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mae'])

        # Display the model architecture
        model.summary()

        # Split the data into training and validation sets
        X_train, X_val, Y_train, Y_val = train_test_split(X_noisy, Y, test_size=validation_split, random_state=42)

        # Train the model
        history = model.fit(X_train, Y_train, epochs=epochs, batch_size=batch_size, validation_data=(X_val, Y_val))

       # self.plot_and_save_training_history(history, plot_filename, "Training and Validation Loss, L-AsparticAcid [0, 1, 1, 0]")

        # Evaluate the model on the validation set
        val_loss, val_mae = model.evaluate(X_val, Y_val)
        print(f"Validation Loss: {val_loss}, Validation MAE: {val_mae}")

        # Make predictions
        predictions = model.predict(X_val)

        # Example: Comparing predictions with actual Y values
        for i in range(5):
            print(f"Predicted: {predictions[i]}, Actual: {Y_val[i]}")

        return model, history

    def tune_model(self, X, Y, hsqc_vector, plot_filename="TUNEDtrainingloss0110.png"):
        input_dim = X.shape[1]
        output_dim = Y.shape[1]

        hypermodel = self.MetaboliteHyperModel(input_dim, output_dim)

        tuner = BayesianOptimization(
            hypermodel.build,
            objective="val_loss",
            max_trials=200,
            executions_per_trial=5,
            directory="tuning_dir",
            project_name=f"metabolite_tuning_{'_'.join(map(str, hsqc_vector))}"
        )

        tuner.search_space_summary()

        early_stopping = EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)

        tuner.search(X, Y, epochs=100, validation_split=0.2, verbose=1, callbacks=[early_stopping])

        tuner.results_summary()

        # Save the best model for this specific HSQC vector
        best_model = tuner.get_best_models(num_models=1)[0]

        # Split the data into training and validation sets
        X_train, X_val, Y_train, Y_val = train_test_split(X, Y, test_size=0.2, random_state=42)

        # history = best_model.fit(X_train, Y_train, epochs=100, validation_data=(X_val, Y_val), verbose=1,
        #                          callbacks=[early_stopping])
        #
        # self.plot_and_save_training_history(history, plot_filename,
        #                                     "Tuned Model: Training and Validation Loss")

        # Evaluate the best model on the validation set
        val_loss, val_mae = best_model.evaluate(X_val, Y_val)
        print(f"Validation Loss: {val_loss}, Validation MAE: {val_mae}")

        # Save the best model and its summary
        self.save_model(best_model, hsqc_vector)
        self.save_model_summary(best_model, val_loss, val_mae, tuner, hsqc_vector)

        # Perform Monte Carlo Dropout predictions after hyperparameter tuning
        mean_pred, std_dev_pred = self.mc_dropout_predict(best_model, X_val, n_iter=10000)

        # Example: Comparing normalized predictions with actual Y values
        for i in range(5):
            print(f"Sample {i+1} - Predicted Mean: {mean_pred[i]}, Standard Deviation: {std_dev_pred[i]}")



        return best_model, X_val, Y_val, mean_pred, std_dev_pred

    class NormalizationLayer(layers.Layer):
        def call(self, inputs):
            # Normalize the inputs to sum to 1
            normalized = tf.math.divide_no_nan(inputs, tf.reduce_sum(inputs, axis=-1, keepdims=True))
            # Scale to sum to 100
            return normalized * 100

    class MetaboliteHyperModel:
        def __init__(self, input_dim, output_dim):
            self.input_dim = input_dim
            self.output_dim = output_dim

        def build(self, hp):
            model = models.Sequential()
            model.add(layers.Input(shape=(self.input_dim,)))
            for i in range(hp.Int('num_layers', 1, 6)):
                model.add(layers.Dense(
                    units=hp.Int(f'units_{i}', min_value=32, max_value=512, step=32, default=64),
                    activation='relu',
                    kernel_regularizer=tf.keras.regularizers.L2(hp.Float('l2_lambda', 1e-5, 1e-2, sampling='log'))
                ))
                model.add(layers.Dropout(rate=hp.Float('dropout_rate', 0.1, 0.5, step=0.05)))
            model.add(layers.Dense(self.output_dim, activation='relu'))
            model.add(IsotopomerAnalysisNN.NormalizationLayer())  # Add the custom normalization layer
            model.compile(
                optimizer=tf.keras.optimizers.Adam(
                    learning_rate=hp.Float('learning_rate', 1e-4, 1e-2, sampling='log', default=0.001)),
                loss='mean_squared_error',
                metrics=['mae']
            )
            return model

    def generate_model_filename(self, hsqc_vector):
        hsqc_str = '_'.join(map(str, hsqc_vector))
        filename = f"model_hsqc_{hsqc_str}.keras"
        return filename

    def save_model(self, model, hsqc_vector, directory="saved_models"):
        if not os.path.exists(directory):
            os.makedirs(directory)
        filename = self.generate_model_filename(hsqc_vector)
        model.save(os.path.join(directory, filename))
        print(f"Model saved as {filename} in {directory}")

    def save_model_summary(self, model, val_loss, val_mae, tuner, hsqc_vector, directory="model_summaries"):
        if not os.path.exists(directory):
            os.makedirs(directory)

        hyperparameters = tuner.get_best_hyperparameters(1)[0]
        summary_data = {
            'Model Name': [self.generate_model_filename(hsqc_vector)],
            'MSE': [val_loss],
            'MAE': [val_mae],
            'num_layers': [hyperparameters.get('num_layers')],
            'learning_rate': [hyperparameters.get('learning_rate')],
            'l2_lambda': [hyperparameters.get('l2_lambda')],
            'dropout_rate': [hyperparameters.get('dropout_rate')]
        }

        summary_df = pd.DataFrame(summary_data)
        summary_filename = os.path.join(directory, f"model_summary_{self.generate_model_filename(hsqc_vector)}.csv")
        summary_df.to_csv(summary_filename, index=False)
        print(f"Model summary saved as {summary_filename}")

    # New method to load data
    def load_hsqc_and_gcms_data(self, hsqc_data_file, gcms_data_file):
        """Loads HSQC and GC-MS data from provided file paths."""
        self.read_hsqc_multiplets(hsqc_data_file)
        self.read_gcms_data(gcms_data_file)
        print("HSQC and GC-MS data loaded successfully.")

    # New method to inspect HSQC and GC-MS data
    def inspect_metabolite_data(self, metabolite_name):
        """Inspects HSQC and GC-MS data for a given metabolite."""
        if metabolite_name in self.exp_multiplets:
            print(f"HSQC Multiplets for {metabolite_name}:")
            for exp_idx, multiplets in enumerate(self.exp_multiplets[metabolite_name]):
                print(f"Experiment {exp_idx + 1}: {multiplets}")

            print(f"\nHSQC Multiplet Percentages for {metabolite_name}:")
            for exp_idx, percentages in enumerate(self.exp_multiplet_percentages[metabolite_name]):
                print(f"Experiment {exp_idx + 1}: {percentages}")

            print(f"\nGC-MS Percentages for {metabolite_name}:")
            for exp_idx, gcms_percentages in enumerate(self.exp_gcms[metabolite_name]):
                print(f"Experiment {exp_idx + 1}: {gcms_percentages}")
        else:
            print(f"No data found for metabolite {metabolite_name}")

    def create_feature_vectors(self, metabolite_name):
        """Combines HSQC and GC-MS data into feature vectors for a given metabolite, handling duplicate multiplets."""
        if metabolite_name not in self.exp_multiplets:
            print(f"No data found for metabolite {metabolite_name}")
            return None

        X_real_data = []

        # Iterate over each experiment
        for exp_idx in range(self.n_exps):
            # Create a dictionary to accumulate percentages for duplicate multiplets
            hsqc_dict = defaultdict(list)

            # Get the HSQC multiplets and their percentages for the current experiment
            hsqc_multiplets = self.exp_multiplets[metabolite_name][exp_idx]
            hsqc_percentages = self.exp_multiplet_percentages[metabolite_name][exp_idx]

            # Fill the dictionary with the multiplet as key and percentages as values
            for multiplet, percentage in zip(hsqc_multiplets, hsqc_percentages):
                hsqc_dict[str(multiplet)].append(percentage)

            # Average the percentages for any duplicate multiplets
            averaged_percentages = [np.mean(hsqc_dict[multiplet]) for multiplet in hsqc_dict]

            # Get the GC-MS percentages for the current experiment
            gcms_percentages = self.exp_gcms[metabolite_name][exp_idx]

            # Combine averaged HSQC percentages and GC-MS data into a single feature vector
            X_sample = np.hstack([averaged_percentages, gcms_percentages])

            # Append the feature vector to the list
            X_real_data.append(X_sample)

        # Convert the list to a numpy array
        X_real_data = np.array(X_real_data)

        return X_real_data

    def load_model_and_predict(self, model_path, X_real_data, n_carbons, n_iter=10000):
        """Loads a trained model, makes predictions on real data with Monte Carlo Dropout, and simulates HSQC/GC-MS data."""

        @tf.keras.utils.register_keras_serializable()
        class NormalizationLayer(tf.keras.layers.Layer):
            def call(self, inputs):
                normalized = tf.math.divide_no_nan(inputs, tf.reduce_sum(inputs, axis=-1, keepdims=True))
                return normalized * 100

        # Load the trained model with custom objects
        best_model = tf.keras.models.load_model(
            model_path,
            custom_objects={'NormalizationLayer': NormalizationLayer}
        )

        # Define a function to run the model with training=True
        def mc_dropout_predict(model, X, n_iter=10000):
            predictions = []
            for _ in range(n_iter):
                predictions.append(model(X, training=True))  # Directly call the model with training=True
            predictions = np.array(predictions)
            return np.mean(predictions, axis=0), np.std(predictions, axis=0)

        # Obtain mean and standard deviation using Monte Carlo Dropout
        mean_predictions, std_dev_predictions = mc_dropout_predict(best_model, X_real_data, n_iter=n_iter)

        print("Mean Predictions:", mean_predictions)
        print("Standard Deviation of Predictions:", std_dev_predictions)

        # Convert mean predictions into isotopomer distributions
        predicted_distributions = []
        possible_isotopomers = [map(int, bin(i)[2:].zfill(n_carbons)).tolist() for i in range(2 ** n_carbons)]

        for prediction in mean_predictions:
            # Filter out isotopomers with zero percentages
            filtered_isotopomers = []
            filtered_percentages = []
            for iso, perc in zip(possible_isotopomers, prediction):
                if perc > 0:
                    filtered_isotopomers.append(iso)
                    filtered_percentages.append(perc)
            predicted_distributions.append({
                'isotopomers': filtered_isotopomers,
                'percentages': filtered_percentages
            })

        return mean_predictions, std_dev_predictions, predicted_distributions

    def simulate_from_predictions(self, predicted_distributions, hsqc_vector):
        """Simulates HSQC and GC-MS data from predicted distributions."""
        predicted_hsqc_data, predicted_gcms_data = [], []

        all_possible_multiplets = self.generate_possible_hsqc_multiplets(hsqc_vector)

        for distribution in predicted_distributions:
            self.init_metabolite('PredictedMetabolite', hsqc_vector)
            self.sim_hsqc_data(metabolite='PredictedMetabolite', exp_index=0, isotopomers=distribution['isotopomers'],
                               percentages=distribution['percentages'])
            self.sim_gcms_data(metabolite='PredictedMetabolite', exp_index=0)

            hsqc_multiplets = self.exp_multiplets['PredictedMetabolite'][0]
            hsqc_percentages = self.exp_multiplet_percentages['PredictedMetabolite'][0]
            gcms_percentages = self.exp_gcms['PredictedMetabolite'][0]

            # Ensure all possible HSQC multiplets are represented
            hsqc_dict = {str(multiplet): 0 for multiplet in all_possible_multiplets}
            for multiplet, perc in zip(hsqc_multiplets, hsqc_percentages):
                hsqc_dict[str(multiplet)] = perc

            uniform_hsqc_percentages = [hsqc_dict[str(multiplet)] for multiplet in all_possible_multiplets]

            predicted_hsqc_data.append({
                'multiplets': all_possible_multiplets,
                'percentages': uniform_hsqc_percentages
            })
            predicted_gcms_data.append(gcms_percentages)

        return predicted_hsqc_data, predicted_gcms_data


    def mc_dropout_predict(self, model, X, n_iter=10000):

        predictions = []

        for _ in range(n_iter):
            # Ensure dropout is active during prediction
            pred = model(X, training=True)  # Pass training=True directly to the model
            predictions.append(pred)

        predictions = np.array(predictions)
        mean_prediction = np.mean(predictions, axis=0)
        std_dev_prediction = np.std(predictions, axis=0)

        return mean_prediction, std_dev_prediction

    def combine_hsqc_gcms(self, predicted_hsqc_data, predicted_gcms_data):
        combined_vector = []
        for hsqc, gcms in zip(predicted_hsqc_data, predicted_gcms_data):
            combined_vector.append(np.hstack([hsqc['percentages'], gcms]))
        return np.array(combined_vector)

    def save_results_summary(self, X_real_data, predicted_distributions, std_dev_predictions,
                             predicted_hsqc_data, predicted_gcms_data, hsqc_vector):
        # Collect results in a dictionary format
        results_data = {
            "Sample": [],
            "Isotopomers": [],
            "Predicted Isotopomer Distribution": [],
            "Standard Deviation": [],
            "HSQC Multiplets": [],
            "Real HSQC %": [],
            "Back Calculated Sim HSQC %": [],
            "Real GC-MS %": [],
            "Back Calculated Sim GC-MS %": [],
        }

        # Determine correct length to split HSQC and GC-MS data
        gcms_length = len(hsqc_vector) + 1  # Adding 1 for the unlabelled carbon
        slice_length = len(X_real_data[0]) - gcms_length
        np.set_printoptions(suppress=True, precision=3)

        # Capture the data
        for i, (pred_dist, std_dev, hsqc, gcms, back_hsqc, back_gcms) in enumerate(zip(
                predicted_distributions, std_dev_predictions, predicted_hsqc_data,
                predicted_gcms_data, X_real_data[:, :slice_length], X_real_data[:, slice_length:]
        )):
            results_data["Sample"].append(i + 1)
            results_data["Isotopomers"].append(pred_dist['isotopomers'])
            results_data["Predicted Isotopomer Distribution"].append(np.round(pred_dist['percentages'], 3))
            results_data["Standard Deviation"].append(np.round(std_dev, 3))
            results_data["HSQC Multiplets"].append(hsqc['multiplets'])
            results_data["Real HSQC %"].append(np.round(back_hsqc, 3))  # Correctly placing the back-calculated data
            results_data["Back Calculated Sim HSQC %"].append(
                np.round(hsqc['percentages'], 3))  # Correctly placing the real data
            results_data["Real GC-MS %"].append(np.round(back_gcms, 3))  # Correctly placing the back-calculated data
            results_data["Back Calculated Sim GC-MS %"].append(np.round(gcms, 3))  # Correctly placing the real data

        # Convert the dictionary to a DataFrame
        results_df = pd.DataFrame(results_data)

        # Define the directory and ensure it exists
        results_directory = "nn_analysis_results"
        if not os.path.exists(results_directory):
            os.makedirs(results_directory)

        # Convert HSQC vector to string for the filename
        hsqc_vector_str = "_".join(map(str, hsqc_vector))
        results_filename = f"results_summary_hsqc_{hsqc_vector_str}.xlsx"
        results_filepath = os.path.join(results_directory, results_filename)

        # Save the results to an Excel file
        results_df.to_excel(results_filepath, index=False)
        print(f"Results successfully saved to {results_filepath}")

    import matplotlib.pyplot as plt

    # Plotting training history with improvements
    def plot_and_save_training_history(self, history, filename, title):
        plt.figure(figsize=(10, 6))
        plt.plot(history.history['loss'], label='Training Loss', color='blue', linewidth=2)
        plt.plot(history.history['val_loss'], label='Validation Loss', color='orange', linestyle='--', linewidth=2)
        plt.xlabel('Epochs', fontsize=14)
        plt.ylabel('Loss', fontsize=14)
        plt.title(title, fontsize=16)
        plt.legend(fontsize=12)
        plt.grid(True)
        plt.savefig(filename, dpi=300)  # Save the figure as an image file
        plt.show()


