import os
import pandas as pd
import numpy as np
import scipy as sp
import os
import math
import time
import metabolabpy.__init__ as ml_version


class IsotopomerAnalysis:

    def __init__(self):
        self.nmr_multiplets = pd.DataFrame()
        self.nmr_tocsy_multiplets = pd.DataFrame()
        self.gcms_data = pd.DataFrame()
        self.lcms_data = pd.DataFrame()
        self.nmr_1d_data = pd.DataFrame()
        self.ver = ml_version.__version__
        self.nat_abundance = 1.07  # [%]
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
        # end __init__

    def __str__(self):  # pragma: no cover
        r_string = '______________________________________________________________________________________\n'
        r_string += '\nMetaboLabPy Isotopomer Data Analysis (v. ' + self.ver + ')\n'
        r_string += '______________________________________________________________________________________\n\n'
        return r_string
        # end __str__

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

        self.metabolites = sorted(list(set(self.metabolites)))
        self.n_exps = int(len(self.nmr_multiplets[self.metabolites[0]].keys()    )/6)
        return
    # end read_hsqc_multiplets

    def read_nmr1d_data(self, file_name=''):
        if len(file_name) == 0:
            return

        self.nmr_1d_data = pd.read_excel(file_name, sheet_name=None, keep_default_na=False)
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


        self.metabolites = sorted(list(set(self.metabolites)))
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
            print('Usage:\nia.set_fit_isotopomers(metabolite="L-LacticAcid", isotopomers=[[0, 0, 1], [0, 1, 1]], percentages=[3, 5]')
            return

        if len(isotopomers) != len(percentages):
            print('length of percentages vector does not match number of isotopomers')
            return

        self.reset_fit_isotopomers(metabolite)
        for k in range(len(isotopomers)):
            self.fit_isotopomers[metabolite].append(isotopomers[k])
            self.isotopomer_percentages[metabolite].append(percentages[k])

        zero_isotopomer = list(np.zeros(len(self.fit_isotopomers[metabolite][0]), dtype=int))
        if zero_isotopomer not in self.fit_isotopomers[metabolite]:
            self.fit_isotopomers[metabolite].append(zero_isotopomer)
            self.isotopomer_percentages[metabolite].append(0.0)

        p_sum = sum(self.isotopomer_percentages[metabolite])
        idx = self.fit_isotopomers[metabolite].index(zero_isotopomer)
        if  p_sum < 100.0:
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
            pp = self.isotopomer_percentages[metabolite][k] * (1.0 - n_zeros*self.nat_abundance / 100.0)
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

    def set_gcms_isotopomers(self):
        if len(self.gcms_data.keys()) == 0:
            return

    # end set_gcms_isotopomers
