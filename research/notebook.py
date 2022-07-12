# <markdowncell>
# # Figures for Martijn / RSF / HCO3/Co2


# <codecell>
# set working directory to main directory of package
import codecs
import imp

from matplotlib.pyplot import ylabel
from set_cwd_to_project_root import project_root
from pathlib import Path

import pandas as pd
import numpy as np
import plotting

from chemistry_tools.formulae.formula import Formula
from phreeqpython.phreeqpython import PhreeqPython
from phreeqc_databases_and_definitions.cdmusic_surface import CDMusicSurface
from phreeqc_databases_and_definitions.common_functions import run_over_ph_range

current_file_dir = Path(__file__).parent
pd.options.plotting.backend = "plotly"
%load_ext autoreload
%autoreload 2

# <markdowncell>
# ## Ferrihydrite surface
# <codecell>
pp=PhreeqPython()#project_root / 'phreeqc_databases_and_definitions' / 'phreeqc.dat')
sol_dict = {
  'Temp':      12,
  'pH':        8,
  'units':     'mg/kgw',
  'Amm': '3.3 mg/kgw',
  'Ca': '121.6 mg/kgw',
  'Cl': '32.8 mg/kgw',
  '[Fe+2]': '11.9 mg/kgw',
  'K':         "2 mg/kgW",
  'Mg': '10.0 mg/kgw',
  'Na': '17.9 mg/kgw',
  'P': "0.98 mg/kgw",
  'Si': '17.15 mg/kgw as SiO3',
  'S(6)': '0.002 mg/kgw as SO4',
  'Alkalinity': "437.0 mg/kgW as HCO3 ",
#   'Alkalinity': "437.0 mg/kgW as HCO3 ",
}
sol_dict_tank=sol_dict.copy(   )
sol_dict_tank['Temp'] = 20
sol_dict_tank['Alkalinity'] = '290 mg/kgW as HCO3 '
sol_dict_tank['Ca'] = 127
sol_dict_tank['pH'] = 7.3

sol_ruw=pp.add_solution(sol_dict)
sol_tank=pp.add_solution(sol_dict_tank)

# pp.ip.run_string("USE SOLUTION 0\n" + hfo.phreeqc_str + "\n SAVE SURFACE 0\nSAVE SOLUTION 0")
# pp.ip.run_string("USE SOLUTION 0\nUSE SURFACE 0\nEND\nSAVE SOLUTION 1")
plotting.createSpeciesPlot(solutions=[sol_ruw, sol_tank], render_plot=False, )#only_contain_species='F' )
# plotting.createSIPlot(solution=[sol_ruw, sol_tank], labels=['ruw','tank'], render_plot=False, )#only_contain_species='F' )

# <markdowncell>
# ## Check for charge balance (crude way)
# <codecell>
sol_charge_bal_dict = sol_dict.copy()
sol_charge_bal_dict['Na'] = f"{sol_charge_bal_dict['Na']} charge"
sol_charge_bal = pp.add_solution(sol_charge_bal_dict)
(sol_charge_bal.elements['Na'] - sol_tank.elements['Na']) / Formula.from_string('Na').mass
# <markdowncell>
# ## Let calcite precipitate
# <codecell>

sol_tank_precip = sol_tank.copy()
sol_tank_precip.saturate('calcite')
plotting.createSIPlot(solutions=[sol_ruw, sol_tank, sol_tank_precip], labels=['ruw','tank', 'tank after calcite precip'], render_plot=False)#, min_concentration=-2)#only_contain_species='F' )
# <markdowncell>
# species plot
# <codecell>
plotting.createSpeciesPlot(solutions=[sol_ruw, sol_tank, sol_tank_precip], render_plot=False, )#only_contain_species='F' )
# <markdowncell>
# change ph
sol_tank_co2_equil = sol_tank_precip.copy()
air     = pp.add_gas({'CO2(g)': 0.00316}, volume=1000, pressure=1, fixed_pressure=False, fixed_volume=True, ) # atmosphere
sol_tank_co2_equil.interact(air.copy())
plotting.createSpeciesPlot(solutions=[sol_ruw, sol_tank, sol_tank_precip, sol_tank_co2_equil],
                           labels=['ruw', 'tank', 'precip. calc',
                                   'precip calc+co2 eq.'],
                           render_plot=False, )  # only_contain_species='F' )
# <codecell>
plotting.createSIPlot(solutions=[sol_ruw, sol_tank, sol_tank_precip, sol_tank_co2_equil],
                           labels=['ruw', 'tank', 'precip. calc',
                                   'precip calc+co2 eq.'],
                           render_plot=False, )  # only_contain_species='F' )

# <codecell>
sol_changing_ph = sol_tank.copy()
phs = [7.3, 7, 6.7, 6.4]
si_calcites = []
for ph in phs:
    sol_changing_ph.change_ph(ph)
    si_calcites.append(sol_changing_ph.si('Calcite'))

sol_tank_ph67 = sol_tank.copy().change_ph(6.7).interact(air.copy())
sol_tank_precip.I
plotting.createSpeciesPlot(solutions=[sol_ruw, sol_tank, sol_tank_ph67], render_plot=False, )#only_contain_species='F' )
# <codecell>
plotting.createSIPlot(solutions=[sol_ruw, sol_tank, sol_tank_ph67], render_plot=False, )#only_contain_species='F' )

# <codecell>
sol_eq_ca_co2 = {}
# sol_eq_ca_co2.saturate(['Calcite', 'CO2(g)'], [0, -3])
for fracion in [0.10, 0.20, 0.30, 1.00]:
    key = f'{fracion*100}% air'
    sol_eq_ca_co2[key] = sol_tank.copy()
    sol_eq_ca_co2[key].equalize(['Calcite', 'CO2(g)', 'O2(g)'],
                       [0, np.log10(fracion*0.039), np.log10(fracion*0.2)],
                       [0, 10, 10])
plotting.createSIPlot(solutions=[sol_tank] + list(sol_eq_ca_co2.values()),
                           labels=['tank'] + list(sol_eq_ca_co2.keys()),
                           render_plot=False,
                           min_concentration=-5 )  # only_contain_species='F' )

# <codecell>
plotting.createSpeciesPlot(solutions=[sol_tank] + list(sol_eq_ca_co2.values()),
                           labels=['tank'] + list(sol_eq_ca_co2.keys()),
                           render_plot=False,
                           min_concentration=0,
                           only_contain_species='C'
                           )
# <codecell>
pd.DataFrame({k: [v.si('Calcite'), v.pH, v.total('HCO3-'), v.total('CO2'), v.total('Ca+2')]
             for k, v in {'tank': sol_tank, **sol_eq_ca_co2}.items()},
             index=['SI', 'pH', 'HCO3-','CO2', 'Ca'])


# <markdowncell>
# Amount of calcium removed
# <codecell>
change_in_ca = sol_tank.elements['Ca'] - sol_eq_ca_co2['30.0% air'].elements['Ca']
change_in_ca * Formula.from_string('Ca').mass

# <markdowncell>
# SI of CO2 and Calcite
# <codecell>
pd.DataFrame({'SI': [sol_tank.si('Calcite'), sol_eq_ca_co2['30.0% air'].si('Calcite')],
              'pH': [sol_tank.pH, sol_eq_ca_co2['30.0% air'].pH]},
              index=['tank', 'precip calc+co2 eq.'])


# <markdowncell>
# # Check how much NaHCO3 needs to be added to get the same pH and IS
# create solution
# <codecell>
sol_eq_co2 = {}
start_ph = 6.0
# sol_eq_co2.saturate(['Calcite', 'CO2(g)'], [0, -3])
for fracion in [0.10, 0.20, 0.30, 1.00]:
    key = f'{fracion*100}% air'
    sol_eq_co2[key] = sol_tank.copy()
    sol_eq_co2[key].change_ph(start_ph)
    sol_eq_co2[key].equalize(['Calcite', 'CO2(g)', 'O2(g)'],
                       [0, np.log10(fracion*0.039), np.log10(fracion*0.2)],
                       [0, 10, 10])
plotting.createSIPlot(solutions=[sol_tank] + list(sol_eq_co2.values()),
                           labels=['tank'] + list(sol_eq_co2.keys()),
                           render_plot=False,
                           min_concentration=-8,
                           # only_contain_species='F'
                           )
# <markdowncell>
# SI, CO2, HCO3, CO2, Ca after first setting it to pH=6.5 and then
# equilibrate with % of air and calcite
# <codecell>
pd.DataFrame({k: [v.si('Calcite'), v.pH, v.total('HCO3-'), v.total('CO2'), v.total('Ca+2'), v.I]
             for k, v in {'tank': sol_tank, **sol_eq_co2}.items()},
             index=['SI', 'pH', 'HCO3-','CO2', 'Ca', 'IS'])

# <markdowncell>
# Change back to original pH of tank (7.2)
# <codecell>

sol_eq_co2_then_ph72 = {}
target_ph = 7.2
# sol_eq_co2_then_ph72.saturate(['Calcite', 'CO2(g)'], [0, -3])
for fracion in [0.10, 0.20, 0.30, 1.00]:
    key = f'{fracion*100}% air'
    sol_eq_co2_then_ph72[key] = sol_eq_co2[key].copy()
    sol_eq_co2_then_ph72[key].change_ph(target_ph)
    # sol_eq_co2_then_ph72[key].equalize(['Calcite', 'CO2(g)', 'O2(g)'],
    #                    [0, np.log10(fracion*0.039), np.log10(fracion*0.2)],
    #                    [0, 10, 10])
plotting.createSIPlot(solutions=[sol_tank] + list(sol_eq_co2_then_ph72.values()),
                           labels=['tank'] + list(sol_eq_co2_then_ph72.keys()),
                           render_plot=False,
                           min_concentration=-8,
                           # only_contain_species='F'
                           )

# <markdowncell>
# SI, CO2, HCO3, CO2, Ca After bring pH to 7.2
# <codecell>
pd.DataFrame({k: [v.si('Calcite'), v.pH, v.total('HCO3-'), v.total('CO2'), v.total('Ca+2'), v.I]
             for k, v in {'tank': sol_tank, **sol_eq_co2_then_ph72}.items()},
             index=['SI', 'pH', 'HCO3-','CO2', 'Ca', 'IS'])
# <markdowncell>
# Add NaCl to get the same IS for 30% and 100% air:
# <codecell>
sol_eq_co2_then_ph72['30.0% air with nacl'] = sol_eq_co2_then_ph72['30.0% air'].copy()
sol = sol_eq_co2_then_ph72['30.0% air with nacl']
dose = 0.01
total_dose = 0
mw_NaCl = Formula.from_string('NaCl').mass
while sol.I < sol_eq_co2_then_ph72['100.0% air'].I:
    total_dose += dose
    sol.add('NaCl', 0.01, 'mmol')

relative_difference = 100 * (sol_eq_co2_then_ph72['30.0% air with nacl'].I - sol_eq_co2_then_ph72['100.0% air'].I)/sol.I
print(f"Relative difference in IS of 100% air and 30% air after adding {total_dose:.2f} mmol (={total_dose*mw_NaCl:.2f} mg/L) is "+
      f"{relative_difference:.2f}% and should be less than 1%")

# <markdowncell>
# SI, CO2, HCO3, CO2, Ca After bring pH to 7.2 and adding NaCl
# <codecell>
pd.DataFrame({k: [v.si('Calcite'), v.pH, v.total('HCO3-'), v.total('CO2'), v.total('Ca+2'), v.I]
             for k, v in {'tank': sol_tank, **sol_eq_co2_then_ph72}.items()},
             index=['SI', 'pH', 'HCO3-','CO2', 'Ca', 'IS'])

# <markdowncell>
# ### Try to get higher pH with only adding NaHCO3 when using 100% air
# <codecell>
target_ph = 7.2
dose = 10. / Formula.from_string('NaHCO3').mass # 10 mg/L in mmol/L
sol = sol_eq_co2['100.0% air'].copy()
phs = []
hco3s = []
dose_steps = np.arange(0, 200)

for dose_i in dose_steps:
    key = f'{dose_i} NaHCO3'
    sol.add('NaHCO3', dose, units='mmol' )
    phs.append( sol.pH)
    hco3s.append(sol.total('HCO3-'))
hco3s = np.array(hco3s)

from plotly import express as px
px.scatter(x=dose*dose_steps * Formula.from_string('NaHCO3').mass, y=phs, labels=dict(y='pH', x='Dose (mg NaHCO3/L)') )
# <codecell>
px.scatter(x=dose*dose_steps * Formula.from_string('NaHCO3').mass, y=hco3s* Formula.from_string('NaHCO3').mass, labels=dict(y='hco3 [mol/L]', x='Dose (mg NaHCO3/L)') )

# <markdowncell>
# - bij ph7 .2 weten we de CO2 en HCO3- verhouding, en we weten HCO3 concentratie dus kunnen total C berekenen. Dus
# verschil dat totaal en hoeveen C-totaal er in het vat zit, weten we hoveel C we moeten toevpegen (maakt niet uit in welke vorm)
# dan stellen ph (en zorgen dat CO2 niet ontsnapt) geeft de gewenste pH en HCO3 (en IS?)

# <markdowncell>
# - with dose of 700-1000 mg/L NaHCO3 kom je op een pH van 7.2-7.3 with the 100% air solution
# - with dose of 30-110 mg/L NaHCO3 kom je op een pH van 7.2-7.3 with the 30% air solution
# recreate those solutions (in a dict)
mw_nahco3 = Formula.from_string('NaHCO3').mass
na_hco3_dosing = {'tank': sol_tank.copy(),
                  '30.0% air': sol_eq_ca_co2['30.0% air'].copy().add('NaHCO3', 30 / mw_nahco3, 'mmol'),
                  '100.0% air': sol_eq_ca_co2['100.0% air'].copy().add('NaHCO3', 900 / mw_nahco3, 'mmol')}
# <markdowncell>
# SI, CO2, HCO3, CO2, Ca After bring pH to 7.2 with only HCO3
# <codecell>
pd.DataFrame({k: [v.si('Calcite'), v.pH, v.total('HCO3-'), v.total('CO2'), v.total('Ca+2'), v.I]
             for k, v in na_hco3_dosing.items()},
             index=['SI', 'pH', 'HCO3-','CO2', 'Ca', 'IS'])


# <markdowncell>
# # Synthetic water
# Repeat what we have done above for synthetic water instead of the real water from Spannenburg

# <markdowncell>
# # Create the original synthetic water at pH=3.5
# <codecell>
sol_synth_dict = {
  'Temp':      20,
  'pH':        '7 charge',
  'units':     'mmol/kgw',
  'Amm': '0.18',
  'Ca': '3.03',
  'Cl': '7.4',
  #'[Fe+2]': '0.21',
  'K':         '0.05',
  'Mg': '0.42',
  'Mn': '0.009',
  'Na': '0.52 ',
  'P': "0.01 ",
  'Si': '0.23',
  'S(6)': '0.02',
  'N(3)': '0.154'

#   'Alkalinity': "437.0 mg/kgW as HCO3 ",
}
sol_synt = pp.add_solution(sol_synth_dict)
sol_synths = dict(original= sol_synt)
sol_synt_after_hcl = sol_synt.copy()
volume_reactor = 60# liter
hcl_dose = 50 # mmol , originally 50 mmol
na_hco3_dose = 0.1 # mmol / L
target_inorg_C = 5.41 # mmol / L

sol_synths['after hcl'] = sol_synt.copy().add('HCl', hcl_dose/volume_reactor, 'mmol') # 50 mmol HCl toegevoegd aan reactor met volume van 60 Liter
sol_synths['after hcl and bubble CO2'] = sol_synths['after hcl'].copy()
sol_synths['after hcl and bubble CO2'].equalize(['Calcite', 'CO2(g)', 'O2(g)'],
                    [0, np.log10(fracion*0.039), np.log10(fracion*0.2)],
                    [0, 10, 10])
sol_synths['after hcl, bubble CO2 add HCO3'] = sol_synths['after hcl and bubble CO2'].copy()
inorg_C = sol_synths['after hcl, bubble CO2 add HCO3'].elements['C(4)'] * 1000 # mmol / L
extra_inorg_C = target_inorg_C - inorg_C
sol_synths['after hcl, bubble CO2 add HCO3'].add('NaHCO3', extra_inorg_C, 'mmol')
sol_synths['after hcl, bubble CO2 add HCO3'].change_ph(7.2)


mw_hco3 = Formula.from_string('HCO3').mass
mw_co2 = Formula.from_string('CO2').mass
mw_ca = Formula.from_string('Ca').mass
mw_na = Formula.from_string('Na').mass
mw_cl = Formula.from_string('Cl').mass
pd.DataFrame({k: [v.si('Calcite'), v.pH, v.total('HCO3', 'mmol') * mw_hco3, v.total('CO2', 'mmol')*mw_co2, v.total('Ca+2','mmol')*mw_ca,
                  v.I, v.total('Na', 'mmol') * mw_na, v.total('Cl', 'mmol') * mw_cl]
             for k, v in sol_synths.items()},
             index=['SI', 'pH', 'HCO3-[mg/L]','CO2[mg/L]', 'Ca[mg/L]', 'IS', 'Na[mg/L]', 'Cl [mg/L]'])

# <markdowncell>
# Amount of NaHCO3 added:
# <codecell>
print(f"{na_hco3_dose:.2f} mmol/L (={na_hco3_dose*mw_nahco3:.2f} mg/L) NaHCO3 added to the reactor")
# %%
