import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from ixpy.hsdmix import HSDMIX
from ixpy.colloc import build_collocation  # optional
from pathlib import Path

data_in = "Actual_exp_Data_Provided_1.xlsx"
sheets = pd.read_excel(data_in, sheet_name=None)

# Breakthrough data sheet
brth_df = sheets["PFAS_BrTh"].copy()

# Detect PFAS columns automatically
skip_cols = {"time", "BV-C/Co"}
pfas_cols = [c for c in brth_df.columns if c not in skip_cols]

print("PFAS found:", pfas_cols)

# =========================================================
# HELPER FUNCTIONS


def rmse(yhat, y):
    yhat = np.asarray(yhat, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(yhat) & np.isfinite(y)

    if m.sum() < 5:
        return np.inf

    return float(np.sqrt(np.mean((yhat[m] - y[m]) ** 2)))


def normalize_breakthrough_ratio(y):
    y = np.asarray(y, dtype=float)
    y = np.where(np.isfinite(y), y, np.nan)
    return np.clip(y, 0.0, 2.0)


def clean_exp_series(brth_df, pfas):
    t = brth_df["time"].to_numpy(dtype=float)
    bv = brth_df["BV-C/Co"].to_numpy(dtype=float)
    y = brth_df[pfas].to_numpy(dtype=float)

    m = np.isfinite(t) & np.isfinite(bv) & np.isfinite(y)
    if m.sum() < 8:
        return None

    df = pd.DataFrame({"t": t[m], "bv": bv[m], "y": y[m]}).sort_values("t")

    # Average duplicate time points
    df = df.groupby("t", as_index=False).mean()

    t_u = df["t"].to_numpy(dtype=float)
    if t_u.size < 8 or not np.all(np.diff(t_u) > 0):
        return None

    return {
        "t_s": t_u,
        "BVx": df["bv"].to_numpy(dtype=float) / 1000.0,
        "y": normalize_breakthrough_ratio(df["y"].to_numpy(dtype=float)),
    }


def extend_influent_to_time(model, t_end_s):
    Cin_t = model.Cin_t.copy()
    t_last = float(Cin_t.index.values[-1])

    if t_end_s <= t_last + 1e-12:
        return

    last_row = Cin_t.iloc[-1].copy()
    Cin_t.loc[float(t_end_s)] = last_row
    model.Cin_t = Cin_t.sort_index()


def simulate_bestfit(model, pfas, t_eval_seconds):
    time_mult = float(model.params["time"])
    t_eval_frac = np.asarray(t_eval_seconds, dtype=float) / time_mult

    t_s, u = model.solve(t_eval=t_eval_frac, const_Cin=False, quiet=True)

    ion_index = list(model.ions.index).index(pfas)

    # Outlet concentration
    Cout = u[0, ion_index, -1, :].astype(float)

    # Influent concentration
    C0 = float(model.Cin_t.iloc[0][pfas])
    if C0 > 0:
        y_ratio = Cout / C0
    else:
        y_ratio = np.full_like(Cout, np.nan)

    # Convert model time to bed volumes
    v = float(model.params["v"])
    L = float(model.params["L"])
    BVx = (t_s * v / L) / 1000.0

    # Final outlet concentration
    Ce = float(Cout[-1])

    # Approximate final resin loading
    qe = float(np.mean(u[1:, ion_index, -1, -1]))

    return BVx, y_ratio, Ce, qe


def fit_kxc_ds_robust(pfas):
    exp = clean_exp_series(brth_df, pfas)

    if exp is None:
        print(f"[SKIP] {pfas}: insufficient clean data")
        return None

    if np.nanmax(exp["y"]) < 0.05:
        print(f"[SKIP] {pfas}: negligible breakthrough")
        return None

    t_exp = exp["t_s"]
    BVx_exp = exp["BVx"]
    y_exp = exp["y"]

    base = HSDMIX(data_in)

    if pfas not in base.ions.index:
        print(f"[SKIP] {pfas}: not found in model ions")
        return None

    extend_influent_to_time(base, float(t_exp[-1]))

    ds0 = float(base.ions.loc[pfas, "Ds"]) if "Ds" in base.ions.columns else 1e-9
    k0 = float(base.ions.loc[pfas, "Kxc"]) if "Kxc" in base.ions.columns else 1.0

    x0 = np.log10([ds0, k0])
    bounds = [(-11, -6), (-3, 5)]
    BIG = 1e6

    def objective(x):
        ds = 10 ** x[0]
        kxc = 10 ** x[1]

        try:
            m = HSDMIX(data_in)
            extend_influent_to_time(m, float(t_exp[-1]))

            m.ions.loc[pfas, "Ds"] = ds
            m.ions.loc[pfas, "Kxc"] = kxc

            _, y_model, _, _ = simulate_bestfit(m, pfas, t_exp)

            val = rmse(y_model, y_exp)
            return val if np.isfinite(val) else BIG

        except Exception as e:
            print(
                f"[WARN] objective failed for {pfas} at Ds={ds:.2e}, Kxc={kxc:.3g}: {e}"
            )
            return BIG

    res = minimize(objective, x0, method="L-BFGS-B", bounds=bounds)

    if not res.success:
        print(f"[FAIL] {pfas}: optimization failed -> {res.message}")
        return None

    ds_star = 10 ** res.x[0]
    kxc_star = 10 ** res.x[1]

    # Final best-fit simulation
    m = HSDMIX(data_in)
    extend_influent_to_time(m, float(t_exp[-1]))

    m.ions.loc[pfas, "Ds"] = ds_star
    m.ions.loc[pfas, "Kxc"] = kxc_star

    BVx_model, y_model, Ce, qe = simulate_bestfit(m, pfas, t_exp)

    return {
        "PFAS": pfas,
        "Ds_cm2_s": ds_star,
        "Kxc": kxc_star,
        "RMSE": rmse(y_model, y_exp),
        "Ce": Ce,
        "qe": qe,
        "BVx_exp": BVx_exp,
        "y_exp": y_exp,
        "BVx_model": BVx_model,
        "y_model": y_model,
    }


def plot_fit(result):
    plot_dir = Path("plots")
    plot_dir.mkdir(exist_ok=True)
    
    plt.figure(figsize=(10, 5))
    plt.plot(result["BVx_exp"], result["y_exp"], "o", alpha=0.5, label="Experimental")
    plt.plot(
        result["BVx_model"], result["y_model"], "-", linewidth=2, label="Best-fit model"
    )
    plt.axhline(1.0, linewidth=1)
    plt.xlabel("Bed Volumes")
    plt.ylabel("C/C0")
    plt.title(
        f"{result['PFAS']} | "
        f"Kxc={result['Kxc']:.3g}, "
        f"Ds={result['Ds_cm2_s']:.2e}, "
        f"Ce={result['Ce']:.4g}, "
        f"qe={result['qe']:.4g}"
    )
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    filename = plot_dir / f"{result['PFAS']}_fit.png"
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()


# =========================================================
# MAIN RUN
summary_results = []
full_results = []

for pfas in pfas_cols:
    print(f"\nRunning fit for {pfas} ...")

    result = fit_kxc_ds_robust(pfas)

    if result is None:
        continue

    print(
        f"{result['PFAS']} | "
        f"Kxc={result['Kxc']:.3g} | "
        f"Ds={result['Ds_cm2_s']:.2e} | "
        f"Ce={result['Ce']:.4g} | "
        f"qe={result['qe']:.4g} | "
        f"RMSE={result['RMSE']:.4f}"
    )

    plot_fit(result)

    full_results.append(result)
    summary_results.append(
        {
            "PFAS": result["PFAS"],
            "Kxc": result["Kxc"],
            "Ds_cm2_s": result["Ds_cm2_s"],
            "Ce": result["Ce"],
            "qe": result["qe"],
            "RMSE": result["RMSE"],
        }
    )

# Final summary table
df_results = pd.DataFrame(summary_results).sort_values("RMSE")
print("\nFinal fitted results:")
print(df_results)
