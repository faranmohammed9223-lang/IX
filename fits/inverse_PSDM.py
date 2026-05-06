import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

import PSDM.PSDM as PSDM
import PSDM.PSDM_functions as PSDM_functions

EXCEL_FILE="experimental_6_PSDM.xlsx"
DATA_SHEET="data"
COLUMN_SHEET="columnSpecs"
PROPERTIES_SHEET="Properties"
COMPOUNDS_TO_FIT="ALL"

XN_FIXED=1.0
FIT_RATIO_MIN=0.001
FIT_RATIO_MAX=0.99
MIN_POINTS_TO_FIT=5

MU_CP=0.89
DEFAULT_TORTUOSITY=1.0

NR=14
NZ=19
NE=1
SOLVER="BDF"
WATER_TYPE="Organic Free"
CHEM_TYPE="PFAS"

KIX_MIN=1e-3
KIX_MAX=1e8
SPDFR_MIN=0.1
SPDFR_MAX=30.0
MAX_NFEV=250

N_KIX_START=11
N_SPDFR_START=12
SHOW_PLOTS=True

def normalize_comp_data(df):
    r={str(i).lower().strip() for i in df.index}
    c={str(i).lower().strip() for i in df.columns}
    p={"mw","molarvol","molar volume","molecular weight","molecularweight","vb"}
    return df if r&p else df.T if c&p else df

def prop(df,comp,names):
    df=normalize_comp_data(df)
    idx={str(i).lower().strip():i for i in df.index}
    for n in names:
        k=n.lower().strip()
        if k in idx:
            return float(df.loc[idx[k],comp])
    raise KeyError(f"{names} missing for {comp}")

def get_mw(df,c):return prop(df,c,["MW","MolecularWeight","Molecular Weight"])
def get_mv(df,c):return prop(df,c,["MolarVol","Molar Volume","MolarVolume","Vb"])

def dl(v):return 13.26e-5/((MU_CP**1.14)*(float(v)**0.589))
def dp(dl,t):return float(dl)/float(t)
def kmolar(K,MW):return float(K)/float(MW)/((1/MW)**XN_FIXED)

def read_inputs():
    raw,col_data,comp,_=PSDM_functions.process_input_file(EXCEL_FILE,data_sheet=DATA_SHEET,column_sheet=COLUMN_SHEET)
    comp_data=normalize_comp_data(PSDM_functions.process_input_data(EXCEL_FILE,sheet_name=PROPERTIES_SHEET))
    col=col_data[col_data.columns[0]].copy()
    col.name=col_data.columns[0]
    fit=list(comp) if COMPOUNDS_TO_FIT=="ALL" else list(COMPOUNDS_TO_FIT)
    return raw,col,comp_data,fit

def carbon(col,raw):
    if col.name in raw.columns.levels[0]:return col.name
    return [x for x in raw.columns.levels[0] if x!=col["influentID"]][0]

def series(raw,t,c):return raw[(t,c)].dropna().sort_index().astype(float)

def BV(col):return np.pi*float(col["diam"])**2*float(col["L"])/4
def BV_axis(t,col):return float(col["flrt"])*np.asarray(t)/BV(col)

def base_k(raw,col,comp_data,comps):
    return PSDM.PSDM(col,comp_data,raw,compounds=comps,xn=XN_FIXED,nr=NR,nz=NZ,ne=NE,solver=SOLVER,water_type=WATER_TYPE,chem_type=CHEM_TYPE,optimize=False).k_data.copy()

def breakthrough(comp,raw,col):
    cin=series(raw,col["influentID"],comp)
    cout=series(raw,carbon(col,raw),comp)
    t=cout.index.values.astype(float)
    cin_i=np.interp(t,cin.index.values.astype(float),cin.values.astype(float))
    cout_v=cout.values.astype(float)
    r=cout_v/np.maximum(cin_i,1e-12)
    bv=BV_axis(t,col)
    m=np.isfinite(r)&(r>=FIT_RATIO_MIN)&(r<=FIT_RATIO_MAX)
    return {"t_all":t,"bv_all":bv,"cin_all":cin_i,"ratio_all":r,
            "t_fit":t[m],"bv_fit":bv[m],"cin_fit":cin_i[m],"ratio_fit":r[m],"cout_fit":cout_v[m]}

def transport(comp,col,comp_data):
    dlv=dl(get_mv(comp_data,comp))
    t=float(col.get("tortu",DEFAULT_TORTUOSITY))
    if t<=0:t=DEFAULT_TORTUOSITY
    return dlv,dp(dlv,t),t

def run_psdm(comp,KIX,SPDFR,raw,col,comp_data,kbase):
    c=col.copy();c["psdfr"]=float(SPDFR)
    kd=kbase[[comp]].copy();kd.loc["K",comp]=KIX;kd.loc["1/n",comp]=XN_FIXED
    mw=get_mw(comp_data,comp)
    dlv,dpv,t=transport(comp,c,comp_data)

    m=pd.DataFrame(0.0,index=["kf","dp","ds"],columns=[comp])
    m.loc["dp",comp]=dpv

    model=PSDM.PSDM(c,comp_data,raw,compounds=[comp],k_data=kd,xn=XN_FIXED,
        xn_range=np.array([XN_FIXED]),test_range=np.array([KIX]),
        mass_transfer=m,nr=NR,nz=NZ,ne=NE,solver=SOLVER,
        water_type=WATER_TYPE,chem_type=CHEM_TYPE,optimize=False)

    f=model.run_psdm()[comp]

    return f,{
        "KIX":KIX,"SPDFR":SPDFR,"1/n":XN_FIXED,
        "K_molar":kmolar(KIX,mw),
        "Dl":dlv,
        "Dp":model.mass_transfer_data.loc["dp",comp],
        "Ds":model.mass_transfer_data.loc["ds",comp],
        "kf":model.mass_transfer_data.loc["kf",comp],
        "tortu":t
    }

def residual(x,comp,raw,col,comp_data,kbase):
    KIX,SPDFR=10**x[0],10**x[1]
    exp=breakthrough(comp,raw,col)
    if len(exp["t_fit"])<MIN_POINTS_TO_FIT:return np.ones(1)*1e6
    try:
        f,_=run_psdm(comp,KIX,SPDFR,raw,col,comp_data,kbase)
        m=f(exp["t_fit"])/np.maximum(exp["cin_fit"],1e-12)
        return np.nan_to_num(m-exp["ratio_fit"],nan=1e6)
    except:return np.ones(len(exp["t_fit"]))*1e6

def best_start(comp,raw,col,comp_data,kbase):
    K0=max(float(kbase.loc["K",comp]),KIX_MIN)
    kg=np.logspace(np.log10(max(KIX_MIN,K0/100)),np.log10(min(KIX_MAX,K0*1000)),N_KIX_START)
    sg=np.linspace(SPDFR_MIN,SPDFR_MAX,N_SPDFR_START)
    best=(1e9,K0,5.0)
    for k in kg:
        for s in sg:
            try:
                exp=breakthrough(comp,raw,col)
                f,_=run_psdm(comp,k,s,raw,col,comp_data,kbase)
                r=f(exp["t_fit"])/np.maximum(exp["cin_fit"],1e-12)
                rmse=np.sqrt(np.mean((r-exp["ratio_fit"])**2))
                if rmse<best[0]:best=(rmse,k,s)
            except:pass
    return best[1],best[2]

def plot_best_fit(comp,exp,f,col,meta,rmse):
    t=np.linspace(exp["t_all"].min(),exp["t_all"].max(),600)
    bv=BV_axis(t,col)
    cin=np.interp(t,exp["t_all"],exp["cin_all"])
    ratio=np.asarray(f(t),dtype=float)/np.maximum(cin,1e-12)

    plt.figure(figsize=(8.5,5.5))
    plt.scatter(exp["bv_all"],exp["ratio_all"],s=32,c="0.35",alpha=0.35,label="Experimental data")
    plt.scatter(exp["bv_fit"],exp["ratio_fit"],s=48,c="black",alpha=0.9,label="Fitting points")
    plt.plot(bv,ratio,lw=2.8,c="tab:red",label="Best PSDM fit")
    plt.axhline(FIT_RATIO_MIN,ls="--",lw=1,c="0.5",alpha=0.7)
    plt.axhline(FIT_RATIO_MAX,ls="--",lw=1,c="0.5",alpha=0.7)
    plt.xlabel("Bed volumes treated",fontsize=11)
    plt.ylabel(r"$C_t/C_0$",fontsize=11)
    plt.title(f"{comp} PSDM fit\nKIX={meta['KIX']:.2e}, SPDFR={meta['SPDFR']:.2f}, RMSE={rmse:.4f}",fontsize=12)
    plt.ylim(-0.05,1.10)
    plt.grid(True,alpha=0.25)
    plt.legend(frameon=True)
    plt.tight_layout()
    plt.show()

def fit(comp,raw,col,comp_data,kbase):
    exp=breakthrough(comp,raw,col)
    if len(exp["t_fit"])<MIN_POINTS_TO_FIT:
        return {"compound":comp,"success":False}

    K0,S0=best_start(comp,raw,col,comp_data,kbase)
    bnds=(np.log10([KIX_MIN,SPDFR_MIN]),np.log10([KIX_MAX,SPDFR_MAX]))
    res=least_squares(residual,np.log10([K0,S0]),bounds=bnds,
        args=(comp,raw,col,comp_data,kbase),max_nfev=MAX_NFEV)

    KIX,SPDFR=10**res.x[0],10**res.x[1]
    f,meta=run_psdm(comp,KIX,SPDFR,raw,col,comp_data,kbase)

    m=f(exp["t_fit"])/np.maximum(exp["cin_fit"],1e-12)
    rmse=np.sqrt(np.mean((m-exp["ratio_fit"])**2))

    if SHOW_PLOTS:
        plot_best_fit(comp,exp,f,col,meta,rmse)

    return {
        "compound":comp,
        "success":res.success,
        "RMSE":rmse,
        "KIX":meta["KIX"],
        "SPDFR":meta["SPDFR"],
        "K_molar":meta["K_molar"],
        "Dl":meta["Dl"],
        "Dp":meta["Dp"],
        "Ds":meta["Ds"],
        "kf":meta["kf"],
        "K_base":kbase.loc["K",comp],
        "q_base":kbase.loc["q",comp]
    }

def main():
    raw,col,comp_data,comps=read_inputs()
    kbase=base_k(raw,col,comp_data,comps)
    results=[fit(c,raw,col,comp_data,kbase) for c in comps]
    df=pd.DataFrame(results)
    print(df.to_string(index=False))
    return df

results_df=main()