import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from pathlib import Path

try:
    from ixpy.hsdmix import HSDMIX
except ImportError as e:
    raise ImportError(
        "Could not import HSDMIX from ixpy.\n"
        "Please install the Ion Exchange Model package first.\n\n"
        "Example:\n"
        "python -m pip install /path/to/IonExchangeModel"
    ) from e


# --------------------------------------------------
# User input
# --------------------------------------------------
XLSX_PATH = Path("Actual_exp_Data_Provided_1.xlsx")
SHEET_NAME = "PFAS_BrTh"


# --------------------------------------------------
# Load data
# --------------------------------------------------
sheets = pd.read_excel(XLSX_PATH, sheet_name=None)
brth_df = sheets[SHEET_NAME].copy()

skip_cols = {"time", "BV-C/Co"}
pfas_cols = [col for col in brth_df.columns if col not in skip_cols]


# --------------------------------------------------
# Helper functions
# --------------------------------------------------
def rmse(yhat, y):
    yhat = np.asarray(yhat, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(yhat) & np.isfinite(y)
    if mask.sum() < 5:
        return np.inf

    return float(np.sqrt(np.mean((yhat[mask] - y[mask]) ** 2)))


def normalize_breakthrough_ratio(y):
    y = np.asarray(y, dtype=float)
    y = np.where(np.isfinite(y), y, np.nan)
    return np.clip(y, 0.0, 2.0)


def clean_exp_series(dataframe, pfas):
    t = dataframe["time"].to_numpy(dtype=float)
    bv = dataframe["BV-C/Co"].to_numpy(dtype=float)
    y = dataframe[pfas].to_numpy(dtype=float)

    mask = np.isfinite(t) & np.isfinite(bv) & np.isfinite(y)
    if mask.sum() < 8:
        return None

    df = pd.DataFrame(
        {
            "t": t[mask],
            "bv": bv[mask],
            "y": y[mask],
        }
    ).sort_values("t")

    df = df.groupby("t", as_index=False).mean()

    t_unique = df["t"].to_numpy(dtype=float)
    if t_unique.size < 8 or not np.all(np.diff(t_unique) > 0):
        return None

    return {
        "t_s": t_unique,
        "BVx": df["bv"].to_numpy(dtype=float) / 1000.0,
        "y": normalize_breakthrough_ratio(df["y"].to_numpy(dtype=float)),
    }


def extend_influent_to_time(model, t_end_s):
    cin_t = model.Cin_t.copy()
    t_last = float(cin_t.index.values[-1])

    if t_end_s <= t_last + 1e-9:
        return

    last_row = cin_t.iloc[-1].copy()
    cin_t.loc[float(t_end_s)] = last_row
    model.Cin_t = cin_t.sort_index()


def simulate_outlet_ratio_fixed_time(model, pfas, t_eval_seconds):
    """
    HSDMIX.solve() uses internal time scaling.
    Pass fractional time = seconds / model.params["time"].
    """
    time_mult = float(model.params["time"])
    if time_mult <= 0:
        raise ValueError("model.params['time'] must be > 0")

    t_eval_frac = np.asarray(t_eval_seconds, dtype=float) / time_mult

    t_s, u = model.solve(t_eval=t_eval_frac, const_Cin=False, quiet=True)

    ion_index = list(model.ions.index).index(pfas)

    cout = u[0, ion_index, -1, :].astype(float)
    c0 = float(model.Cin_t.iloc[0][pfas])
    y_ratio = cout / c0 if c0 > 0 else np.full_like(cout, np.nan)

    v = float(model.params["v"])
    L = float(model.params["L"])
    BVx = (t_s * v / L) / 1000.0

    return BVx, y_ratio


# --------------------------------------------------
# Fitting
# --------------------------------------------------
def fit_kxc_ds_robust(
    pfas,
    ds_bounds=(1e-11, 1e-6),
    kxc_bounds=(1e-3, 1e5),
    min_bt=0.05,
):
    exp = clean_exp_series(brth_df, pfas)
    if exp is None:
        return None

    if np.nanmax(exp["y"]) < min_bt:
        return None

    t_exp = exp["t_s"]
    BVx_exp = exp["BVx"]
    y_exp = exp["y"]

    base = HSDMIX(XLSX_PATH)
    if pfas not in base.ions.index:
        return None

    extend_influent_to_time(base, float(t_exp[-1]))

    ds0 = (
        float(base.ions.loc[pfas, "Ds"])
        if "Ds" in base.ions.columns and np.isfinite(base.ions.loc[pfas, "Ds"])
        else 1e-9
    )
    k0 = (
        float(base.ions.loc[pfas, "Kxc"])
        if "Kxc" in base.ions.columns and np.isfinite(base.ions.loc[pfas, "Kxc"])
        else 1.0
    )

    ds0 = float(np.clip(ds0, ds_bounds[0], ds_bounds[1]))
    k0 = float(np.clip(k0, kxc_bounds[0], kxc_bounds[1]))

    x0 = np.array([np.log10(ds0), np.log10(k0)], dtype=float)
    bounds = [
        (np.log10(ds_bounds[0]), np.log10(ds_bounds[1])),
        (np.log10(kxc_bounds[0]), np.log10(kxc_bounds[1])),
    ]

    big_penalty = 1e6

    def objective(x):
        ds = 10 ** x[0]
        kxc = 10 ** x[1]

        model = HSDMIX(XLSX_PATH)
        extend_influent_to_time(model, float(t_exp[-1]))

        model.ions.loc[pfas, "Ds"] = float(ds)
        model.ions.loc[pfas, "Kxc"] = float(kxc)

        try:
            _, y_model = simulate_outlet_ratio_fixed_time(model, pfas, t_exp)
            value = rmse(y_model, y_exp)
            if not np.isfinite(value):
                return big_penalty
            return value
        except Exception:
            return big_penalty

    result = minimize(objective, x0, method="L-BFGS-B", bounds=bounds)

    if not result.success:
        return None

    ds_star = 10 ** result.x[0]
    kxc_star = 10 ** result.x[1]

    final_model = HSDMIX(XLSX_PATH)
    extend_influent_to_time(final_model, float(t_exp[-1]))
    final_model.ions.loc[pfas, "Ds"] = float(ds_star)
    final_model.ions.loc[pfas, "Kxc"] = float(kxc_star)

    BVx_model, y_model = simulate_outlet_ratio_fixed_time(final_model, pfas, t_exp)
    err = rmse(y_model, y_exp)

    return {
        "PFAS": pfas,
        "Ds_cm2_s": float(ds_star),
        "Kxc": float(kxc_star),
        "RMSE": float(err),
        "BVx_exp": BVx_exp,
        "y_exp": y_exp,
        "BVx_model": BVx_model,
        "y_model": y_model,
    }


# --------------------------------------------------
# Plotting
# --------------------------------------------------
def plot_fit(result):
    pfas = result["PFAS"]

    plt.figure(figsize=(12, 6))
    plt.plot(result["BVx_exp"], result["y_exp"], "o", ms=4, alpha=0.35, label=f"{pfas} exp (C/C0)")
    plt.plot(result["BVx_model"], result["y_model"], "-", lw=2.5, label=f"{pfas} model (C/C0)")
    plt.axhline(1.0, lw=2, alpha=0.30, label="influent (C/C0 = 1)")
    plt.xlabel("Bed Volumes (×1000)")
    plt.ylabel("C/C0")
    plt.title(
        f"{pfas}: fit (Kxc, Ds) | Kxc={result['Kxc']:.3g}, "
        f"Ds={result['Ds_cm2_s']:.2e} | RMSE={result['RMSE']:.3g}"
    )
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.show()


# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    rows = []

    for pfas in pfas_cols:
        result = fit_kxc_ds_robust(pfas)

        if result is None:
            print(f"[SKIP/FAIL] {pfas}")
            continue

        print(
            f"{pfas:12s}  "
            f"Kxc={result['Kxc']:.3g}  "
            f"Ds={result['Ds_cm2_s']:.3e}  "
            f"RMSE={result['RMSE']:.4g}"
        )

        plot_fit(result)

        rows.append(
            {
                "PFAS": result["PFAS"],
                "Kxc": result["Kxc"],
                "Ds_cm2_s": result["Ds_cm2_s"],
                "RMSE": result["RMSE"],
            }
        )

    if rows:
        df = pd.DataFrame(rows).sort_values("RMSE")
        print("\n===== SUMMARY (Kxc + Ds fits; robust) =====")
        print(df.to_string(index=False))
    else:
        print("\nNo successful PFAS fits were returned.")


if __name__ == "__main__":
    main()
