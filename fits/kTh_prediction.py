import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_squared_error
import plotly.graph_objects as go
import os

# create output directory
out_dir = "data_out/k_pred"
os.makedirs(out_dir, exist_ok=True)

thomas_out_dir = "data_out/thomas_fit"

syms = [
    "NaCl",
    "NaCl+NOM1",
    "NaCl+SO4",
    "NaCl+NO3",
    "NaCl+NOM5",
    "NaCl+NOM1+SO4",
    "NaCl+NOM1+NO3",
    "NaCl+SO4+NO3",
    "NaCl+NOM1+SO4+NO3",
    "NaCl+SO4 (0.57 mM)",
    "NaCl+NOM1+SO4 (0.57 mM)",
    "NaCl+NOM5+SO4",
]

pfas = [
    "PFMOAA",
    "PMPA",
    "GenX",
    "NBP 1",
    "NBP 2-1",
    "NBP 2-2",
    "PFBA",
    "PFPeA",
    "PFHxA",
    "PFHpA",
    "PFOA",
    "PFNA",
    "PFDA",
    "4:2 FTS",
    "6:2 FTS",
    "8:2 FTS",
    "PFBS",
    "PFHxS",
    "PFOS",
]

colorlist = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#393b79",
    "#637939",
    "#1f77b4",
    "#2ca02c",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#393b79",
    "#637939",
]

marker_symbols = [
    "circle",
    "square",
    "diamond",
    "cross",
    "x",
    "triangle-up",
    "triangle-up",
    "pentagon",
    "hexagon",
    "star",
    "hourglass",
    "bowtie",
    "circle-open",
    "diamond",
    "diamond-open",
    "cross-open",
    "x-open",
    "triangle-up-open",
    "triangle-down-open",
    "pentagon-open",
    "hexagon-open",
    "star-open",
    "hourglass-open",
    "bowtie-open",
]

wq_style_map = {
    wq: {"symbol": marker_symbols[i], "color": colorlist[i % len(colorlist)]}
    for i, wq in enumerate(syms)
}

pfas_style_map = {
    pfas: {"symbol": marker_symbols[i], "color": colorlist[i % len(colorlist)]}
    for i, pfas in enumerate(pfas)
}

df = pd.read_excel(f"{thomas_out_dir}/thomas_params.xlsx")
df_clean = df.dropna(subset=["kTh", "qe"]).copy()

df_clean["log_q"] = np.log10(df_clean["qe"])
df_clean["log_k"] = np.log10(df_clean["kTh"])

train_df, test_df = train_test_split(df_clean, test_size=0.3, random_state=42)

model = LinearRegression()
model.fit(train_df[["log_q"]], train_df["log_k"])

test_df["log_k_pred"] = model.predict(test_df[["log_q"]])
test_df["k_pred"] = 10 ** test_df["log_k_pred"]

fig = go.Figure()

for wq, group in test_df.groupby("WQ"):

    style = wq_style_map.get(wq, {"symbol": "circle", "color": "black"})

    fig.add_trace(
        go.Scatter(
            x=group["kTh"],
            y=group["k_pred"],
            mode="markers",
            name=wq,
            marker=dict(symbol=style["symbol"], color=style["color"], size=15),
            customdata=np.stack([group["PFAS"], group["qe"]], axis=-1),
            hovertemplate=(
                "<b>WQ:</b> " + wq + "<br>"
                "<b>PFAS:</b> %{customdata[0]}<br>"
                "<b>q:</b> %{customdata[1]:.3f}<br>"
                "<b>Actual k:</b> %{x:.3e}<br>"
                "<b>Predicted k:</b> %{y:.3e}<br>"
                "<extra></extra>"
            ),
        )
    )

# 1:1 line
x_line = np.linspace(-1, 5, 400)
var = 0.3

fig.add_trace(
    go.Scatter(
        x=x_line,
        y=x_line,
        mode="lines",
        line=dict(dash="solid", color="black", width=3),
        name="1:1 line",
    )
)
fig.add_trace(
    go.Scatter(
        x=x_line,
        y=(1 + var) * x_line,
        mode="lines",
        line=dict(dash="dot", width=3, color="black"),
        showlegend=False,
    )
)

fig.add_trace(
    go.Scatter(
        x=x_line,
        y=(1 - var) * x_line,
        mode="lines",
        line=dict(dash="dot", width=3, color="black"),
        name=f"±{100*var:.0f}% variance",
    )
)

fig.update_layout(
    font={"family": "Arial", "size": 32, "color": "black"},
    plot_bgcolor="white",
    xaxis=dict(
        title="k<sub>Th</sub> (L/mg-BV)",
        type="log",
        gridcolor="lightgrey",
        linecolor="black",
        ticks="outside",
        mirror=True,
        linewidth=3,
        tickwidth=3,
        title_standoff=10,
        showgrid=False,
        range=[np.log10(0.015), np.log10(4.5)],
        tickvals=[
            0.02,
            0.03,
            0.04,
            0.05,
            0.06,
            0.07,
            0.08,
            0.09,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1,
            2,
            3,
            4,
            5,
        ],
        ticktext=[
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            0.1,
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            1,
            "",
            "",
            "",
            "",
        ],
    ),
    yaxis=dict(
        title="Predicted k<sub>Th</sub>",
        type="log",
        gridcolor="lightgrey",
        linecolor="black",
        ticks="outside",
        mirror=True,
        linewidth=3,
        tickwidth=3,
        title_standoff=10,
        showgrid=False,
        range=[np.log10(0.015), np.log10(4.5)],
        tickvals=[
            0.02,
            0.03,
            0.04,
            0.05,
            0.06,
            0.07,
            0.08,
            0.09,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1,
            2,
            3,
            4,
            5,
        ],
        ticktext=[
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            0.1,
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            1,
            "",
            "",
            "",
            "",
        ],
    ),
    legend=dict(
        title="",
        font={"size": 18},
        bgcolor="rgba(255, 255, 255, 0)",
        orientation="v",
        x=1,
        y=0,
        xanchor="left",
        yanchor="bottom",
    ),
    showlegend=True,
    width=800,
    height=600,
)

fig.write_html(f"{out_dir}/k_pred.html")
fig.show()

seeds = range(0, 1000)

pfas_results = []
wq_results = []

for seed in seeds:

    train_df, test_df = train_test_split(df_clean, test_size=0.3, random_state=seed)

    model = LinearRegression()
    model.fit(train_df[["log_q"]], train_df["log_k"])

    log_k_pred = model.predict(test_df[["log_q"]])
    k_pred = 10**log_k_pred
    k_true = test_df["kTh"].values

    test_df = test_df.copy()

    # --- log deviation (core metric) ---
    test_df["dev"] = 1 - k_pred / k_true
    test_df["abs_dev"] = np.abs(1 - k_pred / k_true)

    # --- PFAS grouping ---
    pfas_group = test_df.groupby("PFAS")[["abs_dev", "dev"]].mean().reset_index()
    pfas_group["seed"] = seed
    pfas_results.append(pfas_group)

    # --- WQ grouping ---
    wq_group = test_df.groupby("WQ")[["abs_dev", "dev"]].mean().reset_index()
    wq_group["seed"] = seed
    wq_results.append(wq_group)

pfas_df = pd.concat(pfas_results)
wq_df = pd.concat(wq_results)

pfas_summary = (
    pfas_df.groupby("PFAS")
    .agg(
        abs_bias=("abs_dev", "mean"),
        abs_uncertainty=("abs_dev", "std"),
        bias=("dev", "mean"),
        uncertainty=("dev", "std"),
    )
    .reset_index()
)

wq_summary = (
    wq_df.groupby("WQ")
    .agg(
        abs_bias=("abs_dev", "mean"),
        abs_uncertainty=("abs_dev", "std"),
        bias=("dev", "mean"),
        uncertainty=("dev", "std"),
    )
    .reset_index()
)

fig = go.Figure()

for pfas, group in pfas_summary.groupby("PFAS"):

    style = pfas_style_map.get(pfas, {"symbol": "circle", "color": "black"})

    fig.add_trace(
        go.Bar(
            x=group["PFAS"],
            y=group["bias"],
            name=pfas,
            error_y=dict(type="data", array=group["uncertainty"], visible=True),
        )
    )

fig.update_layout(
    font={"family": "Arial", "size": 14, "color": "black"},
    plot_bgcolor="white",
    xaxis={
        "gridcolor": "lightgrey",
        "linecolor": "black",
        "ticks": "outside",
        "mirror": True,
        "linewidth": 3,
        "tickwidth": 3,
        "title_standoff": 0,
        "title": "PFAS",
        "showgrid": False,
    },
    yaxis={
        "gridcolor": "lightgrey",
        "linecolor": "black",
        "ticks": "outside",
        "mirror": True,
        "linewidth": 3,
        "tickwidth": 3,
        "title_standoff": 0,
        "title": "Deviation",
        "showgrid": False,
        "tickformat": ".0%",
    },
    showlegend=False,
    height=400,
    width=1000,
)

fig.write_html(f"{out_dir}/pfas_dev.html")
fig.show()

fig = go.Figure()

for wq, group in wq_summary.groupby("WQ"):

    style = wq_style_map.get(wq, {"symbol": "circle", "color": "black"})

    fig.add_trace(
        go.Bar(
            x=group["WQ"],
            y=group["bias"],
            name=wq,
            error_y=dict(type="data", array=group["uncertainty"], visible=True),
        )
    )

fig.update_layout(
    font={"family": "Arial", "size": 14, "color": "black"},
    plot_bgcolor="white",
    xaxis={
        "gridcolor": "lightgrey",
        "linecolor": "black",
        "ticks": "outside",
        "mirror": True,
        "linewidth": 3,
        "tickwidth": 3,
        "title_standoff": 0,
        "title": "Water Matrix",
        "showgrid": False,
    },
    yaxis={
        "gridcolor": "lightgrey",
        "linecolor": "black",
        "ticks": "outside",
        "mirror": True,
        "linewidth": 3,
        "tickwidth": 3,
        "title_standoff": 0,
        "title": "Deviation",
        "showgrid": False,
        "tickformat": ".0%",
    },
    showlegend=False,
    height=500,
    width=1000,
)

fig.write_html(f"{out_dir}/wq_dev.html")
fig.show()

fig = go.Figure()

for pfas, group in pfas_summary.groupby("PFAS"):

    style = pfas_style_map.get(pfas, {"symbol": "circle", "color": "black"})

    fig.add_trace(
        go.Bar(
            x=group["PFAS"],
            y=group["abs_bias"],
            name=pfas,
            error_y=dict(type="data", array=group["abs_uncertainty"], visible=True),
        )
    )

fig.update_layout(
    font={"family": "Arial", "size": 14, "color": "black"},
    plot_bgcolor="white",
    xaxis={
        "gridcolor": "lightgrey",
        "linecolor": "black",
        "ticks": "outside",
        "mirror": True,
        "linewidth": 3,
        "tickwidth": 3,
        "title_standoff": 0,
        "title": "PFAS",
        "showgrid": False,
    },
    yaxis={
        "gridcolor": "lightgrey",
        "linecolor": "black",
        "ticks": "outside",
        "mirror": True,
        "linewidth": 3,
        "tickwidth": 3,
        "title_standoff": 0,
        "title": "Absolute Deviation",
        "showgrid": False,
        "tickformat": ".0%",
    },  #'range':[0,300000],
    #'tickvals': [0, 50000, 100000, 150000, 200000, 250000], 'ticktext': [0, '50k', '100k', '150k', '200k', '>250k']},
    showlegend=False,
    height=400,
    width=1000,
)

fig.write_html(f"{out_dir}/pfas_abs_dev.html")
fig.show()

fig = go.Figure()

for wq, group in wq_summary.groupby("WQ"):

    style = wq_style_map.get(wq, {"symbol": "circle", "color": "black"})

    fig.add_trace(
        go.Bar(
            x=group["WQ"],
            y=group["abs_bias"],
            name=wq,
            error_y=dict(type="data", array=group["abs_uncertainty"], visible=True),
        )
    )

fig.update_layout(
    font={"family": "Arial", "size": 14, "color": "black"},
    plot_bgcolor="white",
    xaxis={
        "gridcolor": "lightgrey",
        "linecolor": "black",
        "ticks": "outside",
        "mirror": True,
        "linewidth": 3,
        "tickwidth": 3,
        "title_standoff": 0,
        "title": "Water Matrix",
        "showgrid": False,
    },
    yaxis={
        "gridcolor": "lightgrey",
        "linecolor": "black",
        "ticks": "outside",
        "mirror": True,
        "linewidth": 3,
        "tickwidth": 3,
        "title_standoff": 0,
        "title": "Absolute Deviation",
        "showgrid": False,
        "tickformat": ".0%",
    },
    showlegend=False,
    height=500,
    width=1000,
)

fig.write_html(f"{out_dir}/wq_abs_dev.html")
fig.show()
