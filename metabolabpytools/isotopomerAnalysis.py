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
        self.fit_isotopomers = np.array([[]], dtype=int)
        self.nat_abundance = 1.07  # [%]
        self.gcms_scaling = 1.0
        self.hsqc_scaling = 1.0
        self.lcms_scaling = 1.0
        self.tocsy_scaling = 1.0
        self.metabolites = []
        self.n_exps = 0
        # end __init__

    def __str__(self):  # pragma: no cover
        r_string = '______________________________________________________________________________________\n'
        r_string += '\nMetaboLabPy Isotopomer Data Analysis (v. ' + self.ver + ')\n'
        r_string += '______________________________________________________________________________________\n\n'
        return r_string
        # end __str__

    def read_hsqc_multiplets(self, file_name=''):
        if len(file_name) == 0:
            return

        self.nmr_multiplets = pd.read_excel(file_name, sheet_name=None, keep_default_na=False)
        self.metabolites = []
        for k in self.nmr_multiplets.keys():
            self.metabolites.append(k)

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
        return
    # end read_gcms_data

    def set_fit_isotopomers(self):
        if len(nmr_tocsy_multiplets.keys()) == 0:
            self.set_hsqc_isotopomers()

        if len(lcms_data.keys()) == 0:
            self.set_gcms_isotopomers()

    # end set_fit_isotopomers

    def set_hsqc_isotopomers(self):
        if len(self.nmr.multiplets.keys()) == 0:
            return

    # end set_hsqc_isotopomers

    def set_gcms_isotopomers(self):
        if len(self.gcms_data.keys()) == 0:
            return

    # end set_gcms_isotopomers
