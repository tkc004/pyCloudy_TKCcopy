#!/usr/bin/env python
# coding: utf-8

# # In this example we use the MdB class to access a database of models. 

# The dabase is 3MdB, described here: https://sites.google.com/site/mexicanmillionmodels/the-different-projects/hii_chim



import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

script_dir = Path(__file__).resolve().parent
temp_model_dir = script_dir / 'temp_models'
temp_model_dir.mkdir(exist_ok=True)
fig_dir = script_dir / 'figures'
fig_dir.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(temp_model_dir / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(temp_model_dir / ".cache"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _example_utils import save_fig




# Defining the connection parameters.
required_env = ('MdB_HOST', 'MdB_USER', 'MdB_PASSWD', 'MdB_DB_17')
missing_env = [name for name in required_env if not os.environ.get(name)]
if missing_env:
    raise RuntimeError(
        'MdB database configuration is missing: ' + ', '.join(missing_env)
    )
host, user, passwd, db = (os.environ[name] for name in required_env)

try:
    from sqlalchemy import create_engine
except ImportError as exc:
    raise RuntimeError(
        'The MdB example requires SQLAlchemy and the PyMySQL driver'
    ) from exc




request = """SELECT
12+oxygen AS OH, 
nitrogen-oxygen AS NO, 
lumi AS logU, 
BLND_372700A/H__1_486133A AS O2, 
BLND_436300A/H__1_486133A AS O3_4363, 
O__3_500684A/H__1_486133A AS O3, 
N__2_658345A/H__1_486133A AS N2, 
(S__2_671644A + S__2_673082A)/H__1_486133A AS S2
FROM tab_17
WHERE ref = 'BOND'
"""
sqlEngine = create_engine(f'mysql+pymysql://{user}:{passwd}@{host}:{3306}/{db}')

with sqlEngine.connect() as db_con:
    res = pd.read_sql(request, con=db_con)




print(len(res))




res




plt.figure(figsize=(10, 8))
plt.scatter(np.log10(res['N2']), np.log10(res['O3']), c=res['logU'], edgecolor = 'none')
plt.xlabel('log [NII]/Ha')
plt.ylabel('log [OIII]/Hb')
cb = plt.colorbar()
cb.set_label('logU');
save_fig(plt.gcf(), fig_dir / 'bond_bpt_logu.png')




plt.figure(figsize=(10, 8))
plt.scatter(np.log10(res['N2']), np.log10(res['O3']), c=res['OH'], edgecolor = 'none')
plt.xlabel('log [NII]/Ha')
plt.ylabel('log [OIII]/Hb')
cb = plt.colorbar()
cb.set_label('O/H');
save_fig(plt.gcf(), fig_dir / 'bond_bpt_oh.png')




with sqlEngine.connect() as db_con:
    res = pd.read_sql("SELECT count(*) as N FROM tab_17 WHERE ref like 'PNe_2020'", con=db_con)
print("Total number of models with ref='PNe_2020': {}".format(res.N.values[0]))




# Query the database
com1 = 'BB' # Blackbody
com2 = 'C' # Constant density
com4 = 'S' # Solar metallicity
com5 = 'N' # No dust
com6 = 1 # selected models
request = f"""SELECT
    A_HYDROGEN_vol_1, A_HELIUM_vol_1, A_HELIUM_vol_2, A_CARBON_vol_2, A_NITROGEN_vol_1, A_OXYGEN_vol_1,A_OXYGEN_vol_2,
    A_NEON_vol_2, A_NEON_vol_4, A_SULPHUR_vol_1, A_SULPHUR_vol_2, A_CHLORINE_vol_1, A_CHLORINE_vol_2, A_CHLORINE_vol_3,
    A_ARGON_vol_2, A_ZINC_vol_3, A_IRON_vol_2, A_NICKEL_vol_2, MassFrac, atm1
FROM tab_17, abion_17
WHERE tab_17.ref like 'PNe_2020'
    AND tab_17.N = abion_17.N
    AND com1 = '{com1}'
    AND com2 = '{com2}'
    AND com4 = '{com4}'
    AND com5 = '{com5}' 
    AND com6 = {com6}
    """
with sqlEngine.connect() as db_con:
    res = pd.read_sql(request, con=db_con)




print(request)




print(len(res))




plt.figure(figsize=(10, 8))
plt.scatter(res['A_OXYGEN_vol_2']/(res['A_OXYGEN_vol_1']+res['A_OXYGEN_vol_2']), 
                np.log10(res['A_OXYGEN_vol_1']/res['A_NITROGEN_vol_1']), c=res['atm1'])
plt.xlabel(r'O$^{++}$/(O$^+$+O$^{++}$)')
plt.ylabel(r'log ICF$_{th}$(N$^+$/O$^+$)')
cb = plt.colorbar()
cb.set_label('Stellar Temperature')
save_fig(plt.gcf(), fig_dir / 'pne_ionization_correction.png')
