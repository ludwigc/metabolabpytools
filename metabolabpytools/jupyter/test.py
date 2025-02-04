from metabolabpytools import isotopomerAnalysis1
ia = isotopomerAnalysis1.IsotopomerAnalysis1()
isotopomers = [[0,1,1],[1,0,1],[0,0,1]]
percentages = [18.0, 10.0,2.0]
hsqc = [0,1,1]
metabolite = 'L-LacticAcid'
exp_index = 0
ia.init_metabolite(metabolite, hsqc)
ia.set_fit_isotopomers(metabolite=metabolite, isotopomers=isotopomers, percentages=percentages)
print(f'Isotopomers : {ia.fit_isotopomers[metabolite]}\nIsotopomer %: {ia.isotopomer_percentages[metabolite]}')
ia.sim_hsqc_data(metabolite=metabolite, exp_index=exp_index, isotopomers=ia.fit_isotopomers[metabolite], percentages=ia.isotopomer_percentages[metabolite])

