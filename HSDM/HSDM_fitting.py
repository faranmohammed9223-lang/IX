

get_ipython().system('python3 -m pip install "IXPY files Path"') ######
get_ipython().system('{sys.executable} -m pip install --upgrade openpyxl')
if not hasattr(DataFrame, "map"):
    DataFrame.map = DataFrame.applymap 
print("pandas version:", pd.__version__)
print("DataFrame.map exists:", hasattr(DataFrame, "map"))



import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from ixpy.hsdmix import HSDMIX
from ixpy.colloc import build_collocation


fn = "/Users/(note:file path)/Actual_exp_Data_Provided_1.xlsx"

# Load all sheets
sheets = pd.read_excel(fn, sheet_name=None)

# Breakthrough data
brth_df = sheets["PFAS_BrTh"].copy()

# PFAS columns
skip_cols = {"time", "BV-C/Co"}
pfas_cols = [c for c in brth_df.columns if c not in skip_cols]
print("PFAS found:", pfas_cols)

##### functions needed for the optimization 
def rmse(yhat, y):
    yhat = np.asarray(yhat, float)
    y    = np.asarray(y, float)
    m = np.isfinite(yhat) & np.isfinite(y)
    if m.sum() < 5:
        return np.inf
    return float(np.sqrt(np.mean((yhat[m] - y[m])**2)))

def normalize_breakthrough_ratio(y):
    y = np.asarray(y, float)
    y = np.where(np.isfinite(y), y, np.nan)
    return np.clip(y, 0.0, 2.0)

def clean_exp_series(brth_df, pfas):
    t  = brth_df["time"].to_numpy(float)
    bv = brth_df["BV-C/Co"].to_numpy(float)
    y  = brth_df[pfas].to_numpy(float)

    m = np.isfinite(t) & np.isfinite(bv) & np.isfinite(y)
    if m.sum() < 8:
        return None

    df = pd.DataFrame({"t": t[m], "bv": bv[m], "y": y[m]}).sort_values("t")
    df = df.groupby("t", as_index=False).mean()

    t_u = df["t"].to_numpy(float)
    if t_u.size < 8 or not np.all(np.diff(t_u) > 0):
        return None

    return { 
        "t_s": t_u,
        "BVx": df["bv"].to_numpy(float) / 1000.0,
        "y": normalize_breakthrough_ratio(df["y"].to_numpy(float))
    }

def extend_influent_to_time(model, t_end_s):
    Cin_t = model.Cin_t.copy()
    t_last = float(Cin_t.index.values[-1])

    if t_end_s <= t_last + 1e-9:
        return

    last_row = Cin_t.iloc[-1].copy()
    Cin_t.loc[float(t_end_s)] = last_row
    model.Cin_t = Cin_t.sort_index()

def simulate_outlet_ratio(model, pfas, t_eval_seconds):
    time_mult = float(model.params["time"])
    t_eval_frac = np.asarray(t_eval_seconds, float) / time_mult

    t_s, u = model.solve(t_eval=t_eval_frac, const_Cin=False, quiet=True)

    i = list(model.ions.index).index(pfas)

    Cout = u[0, i, -1, :].astype(float)
    C0 = float(model.Cin_t.iloc[0][pfas])
    y_ratio = Cout / C0 if C0 > 0 else np.full_like(Cout, np.nan)

    v = float(model.params["v"])
    L = float(model.params["L"])
    BVx = (t_s * v / L) / 1000.0

    return BVx, y_ratio



def fit_kxc_ds_robust(pfas):
    exp = clean_exp_series(brth_df, pfas)
    if exp is None or np.nanmax(exp["y"]) < 0.05:
        return None

    t_exp = exp["t_s"]
    BVx_exp = exp["BVx"]
    y_exp = exp["y"]

    base = HSDMIX(fn)

    if pfas not in base.ions.index:
        return None

    extend_influent_to_time(base, float(t_exp[-1]))

    ds0 = float(base.ions.loc[pfas, "Ds"]) if "Ds" in base.ions.columns else 1e-9
    k0  = float(base.ions.loc[pfas, "Kxc"]) if "Kxc" in base.ions.columns else 1.0

    x0 = np.log10([ds0, k0])

    bounds = [(-11, -6), (-3, 5)]
    BIG = 1e6

    def obj(x):
        ds, kxc = 10**x[0], 10**x[1]

        m = HSDMIX(fn)
        extend_influent_to_time(m, float(t_exp[-1]))

        m.ions.loc[pfas, "Ds"]  = ds
        m.ions.loc[pfas, "Kxc"] = kxc

        try:
            _, y_m = simulate_outlet_ratio(m, pfas, t_exp)
            val = rmse(y_m, y_exp)
            return val if np.isfinite(val) else BIG
        except:
            return BIG

    res = minimize(obj, x0, method="L-BFGS-B", bounds=bounds)

    if not res.success:
        return None

    ds_star, kxc_star = 10**res.x[0], 10**res.x[1]

    m = HSDMIX(fn)                                # HSDM 
    extend_influent_to_time(m, float(t_exp[-1]))
    m.ions.loc[pfas, "Ds"] = ds_star
    m.ions.loc[pfas, "Kxc"] = kxc_star

    BVx_m, y_m = simulate_outlet_ratio(m, pfas, t_exp)

    return {
        "PFAS": pfas,
        "Ds_cm2_s": ds_star,
        "Kxc": kxc_star,
        "RMSE": rmse(y_m, y_exp),
        "BVx_exp": BVx_exp,
        "y_exp": y_exp,
        "BVx_model": BVx_m,
        "y_model": y_m
    }


def plot_fit(r):
    plt.figure(figsize=(10,5))
    plt.plot(r["BVx_exp"], r["y_exp"], "o", alpha=0.4, label="exp")
    plt.plot(r["BVx_model"], r["y_model"], "-", label="model")
    plt.axhline(1.0, alpha=0.3)

    plt.xlabel("Bed Volumes")
    plt.ylabel("C/C0")
    plt.title(f"{r['PFAS']} | Kxc={r['Kxc']:.3g}, Ds={r['Ds_cm2_s']:.2e}")
    plt.legend()
    plt.grid()
    plt.show()

results = []

for pfas in pfas_cols:
    r = fit_kxc_ds_robust(pfas)

    if r is None:
        print(f"[SKIP] {pfas}")
        continue

    print(f"{pfas} | Kxc={r['Kxc']:.3g} | Ds={r['Ds_cm2_s']:.2e} | RMSE={r['RMSE']:.4f}")

    plot_fit(r)

    results.append({
        "PFAS": r["PFAS"],
        "Kxc": r["Kxc"],
        "Ds_cm2_s": r["Ds_cm2_s"],
        "RMSE": r["RMSE"]
    })

df = pd.DataFrame(results).sort_values("RMSE")
print(df)


def rmse(yhat, y):
    yhat, y = np.array(yhat), np.array(y)
    m = np.isfinite(yhat) & np.isfinite(y)
    return np.sqrt(np.mean((yhat[m] - y[m])**2))

def clean_data(df, pfas):
    t = df["time"].values
    y = df[pfas].values
    m = np.isfinite(t) & np.isfinite(y)
    return t[m], y[m]

def fit_pfas(pfas):
    t_exp, y_exp = clean_data(brth_df, pfas)
    def obj(x):
        ds, kxc = 10**x[0], 10**x[1]
        m = HSDMIX(fn)
        m.ions.loc[pfas, "Ds"] = ds
        m.ions.loc[pfas, "Kxc"] = kxc

        try:
            t_s, u = m.solve(const_Cin=False, quiet=True)
            i = list(m.ions.index).index(pfas)
            Cout = u[0, i, -1, :]
            C0 = m.Cin_t.iloc[0][pfas]
            y_model = Cout / C0
            return rmse(y_model[:len(y_exp)], y_exp)
        except:
            return 1e6

    res = minimize(obj, [-9, 0], bounds=[(-11, -6), (-3, 5)])

    if not res.success:
        return None
    ds, kxc = 10**res.x[0], 10**res.x[1]

    # Final model
    m = HSDMIX(fn)
    m.ions.loc[pfas, "Ds"] = ds
    m.ions.loc[pfas, "Kxc"] = kxc
    t_s, u = m.solve(const_Cin=False, quiet=True)
    i = list(m.ions.index).index(pfas)

    # ---- equalibrium outputs ----
    Ce = float(u[0, i, -1, -1])          # outlet concentration
    qe = float(u[1:, i, -1, -1].mean())  # average loading (simple approx)
    return {"PFAS": pfas,
            "Kxc": kxc,
            "Ds": ds,
            "Ce": Ce,
            "qe": qe}
    
results = []
for pfas in pfas_cols:
    print("Running:", pfas)
    r = fit_pfas(pfas)
    if r is None:
        print("  → failed")
        continue
    print(f"  → Kxc={r['Kxc']:.3g}, Ds={r['Ds']:.2e}, Ce={r['Ce']:.3g}, qe={r['qe']:.3g}")
    results.append(r)
df = pd.DataFrame(results)

print(df)






