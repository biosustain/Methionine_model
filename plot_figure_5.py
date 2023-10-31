# Copyright (C) 2020 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Code for plotting DOE and answering question about Allostery"""
import os
from pathlib import Path
import tempfile

import arviz as az
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

from maud.loading_maud_inputs import load_maud_input
from maud.getting_idatas import get_idata

plt.style.use("ipynb")

MAUD_OUTPUT = Path.cwd() / "results" / "methionine"
MAUD_SIM = Path.cwd() / "results" / "methionine" / "simulate"
PLOT_OUTPUT = Path.cwd() / "figures"


def plot_regulation(regulation_df, true_df, ax):
    sns.violinplot(
        data=regulation_df,
        x="value",
        y="regulation",
        color="#005F73",
        linewidth=1,
        orient="h",
        ax=ax,
        inner="quart",
    )
    for line in ax.lines:
        line.set_marker("")
    ax.set_ylabel("")
    ax.set_xlabel(r"$\ln\frac{\text{Experiment 12}}{\text{Experiment 1}}$")
    ax.axvline(x=0, color="k")
    for i, value in enumerate(true_df):
        ax.plot(value, i, "r")
    return ax


def plot_flux(flux_df, ax):
    sns.violinplot(
        data=flux_df,
        x="flux_train",
        y="experiments",
        color="#005F73",
        linewidth=1,
        orient="h",
        ax=ax,
        inner="quart",
    )
    for line in ax.lines:
        line.set_marker("")
    ax.set_ylabel("")
    ax.set_xlabel("Flux")
    ax.set_yticklabels(["Experiment 1", "Experiment 12"])
    return ax


def get_regulation_df(infd):
    edge = "GNMT1_GNMT"
    enzyme = "GNMT1"
    exp_1 = "dataset1"
    exp_2 = "dataset12"
    regulation_df = pd.DataFrame()
    allostery_samples_1 = (
        infd.posterior.allostery_train.sel({"experiments": exp_1, "edges": edge})
        .to_dataframe()
        .reset_index()
    )
    allostery_samples_2 = (
        infd.posterior.allostery_train.sel({"experiments": exp_2, "edges": edge})
        .to_dataframe()
        .reset_index()
    )
    enz_samples_1 = (
        infd.posterior.conc_enzyme_train.sel({"experiments": exp_1, "enzymes": enzyme})
        .to_dataframe()
        .reset_index()
    )
    enz_samples_2 = (
        infd.posterior.conc_enzyme_train.sel({"experiments": exp_2, "enzymes": enzyme})
        .to_dataframe()
        .reset_index()
    )
    saturation_samples_1 = (
        infd.posterior.saturation_train.sel({"experiments": exp_1, "edges": edge})
        .to_dataframe()
        .reset_index()
    )
    saturation_samples_2 = (
        infd.posterior.saturation_train.sel({"experiments": exp_2, "edges": edge})
        .to_dataframe()
        .reset_index()
    )
    reversibility_samples_1 = (
        infd.posterior.reversibility_train.sel({"experiments": exp_1, "edges": edge})
        .to_dataframe()
        .reset_index()
    )
    reversibility_samples_2 = (
        infd.posterior.reversibility_train.sel({"experiments": exp_2, "edges": edge})
        .to_dataframe()
        .reset_index()
    )
    saturation_df = pd.DataFrame.from_dict(
        {
            "value": np.log(saturation_samples_2["saturation_train"])
            - np.log(saturation_samples_1["saturation_train"]),
            "regulation": "saturation",
        }
    )
    allostery_df = pd.DataFrame.from_dict(
        {
            "value": np.log(allostery_samples_2["allostery_train"])
            - np.log(allostery_samples_1["allostery_train"]),
            "regulation": "allostery",
        }
    )
    reversibility_df = pd.DataFrame.from_dict(
        {
            "value": np.log(reversibility_samples_2["reversibility_train"])
            - np.log(reversibility_samples_1["reversibility_train"]),
            "regulation": "reversibility",
        }
    )
    enzyme_df = pd.DataFrame.from_dict(
        {
            "value": np.log(enz_samples_2["conc_enzyme_train"])
            - np.log(enz_samples_1["conc_enzyme_train"]),
            "regulation": "enzyme",
        }
    )
    regulation_df = pd.concat(
        [
            saturation_df,
            allostery_df,
            # reversibility_df,
            enzyme_df,
        ],
    )
    return regulation_df


def get_true_regulatory_values(infd):
    edge = "GNMT1_GNMT"
    enzyme = "GNMT1"
    exp_1 = "dataset1"
    exp_2 = "dataset12"
    t_enz_reg = (
        infd.posterior.conc_enzyme_train.sel(
            {"experiments": [exp_1, exp_2], "enzymes": enzyme}
        )
        .to_dataframe()
        .reset_index()
    )
    t_allo_reg = (
        infd.posterior.allostery_train.sel(
            {"experiments": [exp_1, exp_2], "edges": edge}
        )
        .to_dataframe()
        .reset_index()
    )
    t_sat_reg = (
        infd.posterior.saturation_train.sel(
            {"experiments": [exp_1, exp_2], "edges": edge}
        )
        .to_dataframe()
        .reset_index()
    )
    true_enzyme_reg = np.log(
        t_enz_reg[t_enz_reg["experiments"] == exp_2]["conc_enzyme_train"][1]
        / t_enz_reg[t_enz_reg["experiments"] == exp_1]["conc_enzyme_train"][0]
    )
    true_allostery_reg = np.log(
        t_allo_reg[t_allo_reg["experiments"] == exp_2]["allostery_train"][1]
        / t_allo_reg[t_allo_reg["experiments"] == exp_1]["allostery_train"][0]
    )
    true_saturation_reg = np.log(
        t_sat_reg[t_sat_reg["experiments"] == exp_2]["saturation_train"][1]
        / t_sat_reg[t_sat_reg["experiments"] == exp_1]["saturation_train"][0]
    )
    return [
        true_saturation_reg,
        true_allostery_reg,
        true_enzyme_reg,
    ]


def get_flux_df(infd):
    reaction = "GNMT"
    exp_1 = "dataset1"
    exp_2 = "dataset12"
    flux_df = (
        infd.posterior.flux_train.sel(
            {"experiments": [exp_1, exp_2], "reactions": reaction}
        )
        .to_dataframe()
        .reset_index()
    )
    flux_df[flux_df["flux_train"] == 0.0] = np.nan
    return flux_df.dropna()


def plot_figure_5():
    mi = load_maud_input(data_path=MAUD_OUTPUT / "user_input")
    csvs_methionine = [
        MAUD_OUTPUT / "samples" / file
        for file in os.listdir(MAUD_OUTPUT / "samples")
        if ".csv" in file
    ]
    idata_methionine = get_idata(csvs_methionine, mi=mi, mode="train")
    regulation_df = get_regulation_df(idata_methionine)
    flux_df = get_flux_df(idata_methionine)
    true_values = az.InferenceData.from_netcdf(MAUD_SIM / "idata.nc")
    true_df = get_true_regulatory_values(true_values)

    # calling plotting scripts
    cm = 1 / 2.54  # centimeters in inches
    sns.set_context(
        "paper", rc={"font.size": 14, "axes.titlesize": 14, "axes.labelsize": 14}
    )

    fig = plt.figure(constrained_layout=True, figsize=(21.0 * cm, 10.5 * cm))
    gs = GridSpec(2, 2, figure=fig)

    ax0 = fig.add_subplot(gs[:, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])

    ax1.set_axis_off()
    ax1.imshow(mpimg.imread(PLOT_OUTPUT / "regulatory_diagram.png"), aspect="auto")
    ax0 = plot_regulation(regulation_df, true_df, ax0)
    ax2 = plot_flux(flux_df, ax2)

    plt.figure(dpi=300)
    fig.tight_layout()
    fig.savefig(PLOT_OUTPUT / "decomposition.png")


if __name__ == "__main__":
    plot_figure_5()
