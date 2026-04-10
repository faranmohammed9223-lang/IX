import sys
from pathlib import Path
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

warnings.simplefilter("ignore")

fn2 = "/Users/mohammedfaran/Desktop/Research project & material/ix_github_files/Github_full/Water_Treatment_Models/IonExchangeModel/PSDM"
sys.path.append(fn2)

import PSDM as psdm_psdm
sys.modules["PSDM"] = psdm_psdm
from psdm import PSDM_tools
from psdm import PSDM
from psdm.PSDM_functions import process_input_file, process_input_data

# -----------------------------
# Load data
fn = Path("/Users/mohammedfaran/Desktop/Research Data/Experimental_Data(PSDM)/experimental_1_PSDM.xlsx")
rawdata_df, column_data, compounds, carbons = process_input_file(
    fn,
    data_sheet="data",
    column_sheet="columnSpecs"
)
comp_data = process_input_data(fn, sheet_name="Properties")

# -----------------------------
# PSDM model
column = PSDM.PSDM(
    column_data[column_data.columns[0]],
    comp_data,
    rawdata_df,
    xn=1.0,
    optimize=False,
    chem_type="PFAS",
    water_type="Organic Free"
)

# -----------------------------
# Show k_data
# -----------------------------
kdata_df = column.k_data.T.reset_index().rename(columns={"index": "compound"})
print("\nPSDM-generated k_data:")
print(kdata_df)

influent_key = column_data.loc["influentID"].values[0]
effluent_key = column_data.columns[0]


# Loop through compounds

for compound in compounds:

    try:
        print("\n==============================")
        print(compound)
        cin_series = rawdata_df[(influent_key, compound)].dropna()
        eff_series = rawdata_df[(effluent_key, compound)].dropna()

        if cin_series.empty or eff_series.empty:
            print("No data found → skipping")
            continue

        # Build dataframe
        df = pd.DataFrame({
            "time": eff_series.index.astype(float),
            "Ct": eff_series.values.astype(float)
        })

        cin = float(cin_series.mean())

        df = df.sort_values("time")

        df["Cin"] = cin
        df["Ce"] = df["Ct"]

# -----------------------------
        # Compute q (FIXED)
        t = df["time"].to_numpy()
        Ct = df["Ct"].to_numpy()

        delta_c = cin - Ct

        cum_int = np.zeros_like(delta_c)   # ✅ FIX

        for i in range(1, len(delta_c)):
            cum_int[i] = (
                cum_int[i-1]
                + 0.5 * (delta_c[i] + delta_c[i-1]) * (t[i] - t[i-1])
            )

        df["q"] = (column.flrt * cum_int) / column.wt

        fit_df = df[(df["Ce"] > 0) & (df["q"] >= 0)].copy()

        if fit_df.empty:
            print("No valid data for fitting → skipping")
            continue

        # -----------------------------
        # Plot using PSDM_tools
        # -----------------------------
        print("\nFreundlich fit:")
        PSDM_tools.isotherm_fit(
            fit_df[["Ce", "q"]],
            isotherm="freundlich",
            plot=True,
            save_plot=False
        )

        print("\nLangmuir fit:")
        PSDM_tools.isotherm_fit(
            fit_df[["Ce", "q"]],
            isotherm="langmuir",
            plot=True,
            save_plot=False
        )


    except Exception as e:
        print(f"{compound} failed:", e)