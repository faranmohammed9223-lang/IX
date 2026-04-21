import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
import sys


fn1 = "/Users/mohammedfaran/Desktop/Research project & material/ix_github_files/Github_full/Water_Treatment_Models/PSDM"  ### psdm
if fn1 not in sys.path:
    sys.path.insert(0, fn1)

import PSDM.PSDM as PSDM

# import psdm

data_in = "experimental_2_PSDM.xlsx"  # file path

chem_data = PSDM.process_input_data(data_in, sheet_name="Properties")
k_data = pd.read_excel(
    data_in, sheet_name="Kdata", index_col=0
)  # extra sheets I added to the Excel file to help save time
kf_val = pd.read_excel(data_in, sheet_name="kf_val", index_col=0)
raw_data, column_info, compounds, carbons = PSDM.process_input_file(
    data_in, data_sheet="data", column_sheet="columnSpecs"
)
COMPOUND = "PMPA"  ####  Targeted species
carbon = carbons[0]
col = column_info[carbon].copy()
print(compounds)
print(col)

###################
# FIX weight (wt)  USING TARGET ebed & SETUP BED VOLUME
target_ebed = 0.381
area = np.pi * float(col["diam"]) ** 2 / 4.0
bedvol = area * float(col["L"])
col["wt"] = (1 - target_ebed) * bedvol * float(col["rhop"])
flrt = float(col["flrt"])
Dl = float(
    np.asarray(kf_val.loc[COMPOUND, "DL (cm^2/s)"]).ravel()[0]
)  # liquid diffusivity
kf = float(np.asarray(kf_val.loc[COMPOUND, "kf (cm/s)"]).ravel()[0])
K_init = float(np.asarray(k_data.loc["K", COMPOUND]).ravel()[0])
xn_init = float(np.asarray(k_data.loc["1/n", COMPOUND]).ravel()[0])
Dp = Dl / float(col["tortu"])  # pore diffusion
bv_factor = (
    flrt / bedvol
)  # We will need this factor to convert time (minutes) to Bed Volumes (BVs)

#######
influent_name = col["influentID"]
effluent_name = carbon
infl_series = raw_data[influent_name][COMPOUND].dropna()
c0 = float(infl_series.max())  # Max safely grabs the spiking concentration
if 0 in raw_data.index:
    raw_data.loc[0, (influent_name, COMPOUND)] = c0
    raw_data.loc[0, (effluent_name, COMPOUND)] = 0.0
eff_series = raw_data[effluent_name][COMPOUND].dropna()

# Extract time and normalized concentration
x_exp_time = np.asarray(eff_series.index, dtype=float)
y_exp = np.asarray(eff_series.values, dtype=float) / c0

# Convert experimental time to Bed Volumes
x_exp_bv = x_exp_time * bv_factor


print("K_init =", K_init)
print("xn_init =", xn_init)
print("Dl =", Dl)
print("kf =", kf)
print("Dp =", Dp)
print("Corrected wt =", col["wt"])
print("ebed =", 1 - float(col["wt"]) / (bedvol * float(col["rhop"])))
print("BV Conversion Factor (BV/min) =", bv_factor)
print("True Experimental C0 =", c0)
print("BVs:", x_exp_bv[:])


#     More coming up!
