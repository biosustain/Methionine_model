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
import argparse
from pathlib import Path
import tempfile
import shutil

import arviz as az
import numpy as np
import pandas as pd
import plotnine as p9
from plotnine import aes, geom_point, ggplot, labs
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

from plotting_functions import get_measurements_of_type

from maud import io

plt.style.use('ipynb')

met_13_path = Path.cwd() / "results" / "methionine_cycle_13"
met_13_true_values = Path.cwd() / "simulate"
PLOT_OUTPUT = Path.cwd() / "plots" / "figure_5"

def plot_regulation(regulation_df, true_df, ax):
    sns.violinplot(data=regulation_df, x="value", y="regulation", color="#005F73",
              linewidth=1, orient="h", ax=ax, inner="quart", marker="")
    for line in ax.lines:
        line.set_marker(None)
    ax.set_xlabel(r"ln $\frac{Regulatory\ Components_{2}}{Regulatory\ Components_{1}}$")
    ax.axvline(x=0, color="k")
    for i, value in enumerate(true_df):
        ax.plot(value, i, "r")
    return ax

def plot_flux(flux_df, ax):
    sns.violinplot(data=flux_df, x="flux", y="experiments", color="#005F73",
              linewidth=1, orient="h", ax=ax, inner="quart", marker="")
    for line in ax.lines:
        line.set_marker(None)
    return ax



def get_regulation_df(infd):
    enzyme = "GNMT1"
    exp_1 = "dataset_1"
    exp_2 = "dataset_2"
    steady_state = [not div 
                    for div 
                    in infd.sample_stats.diverging.to_dataframe().reset_index()["diverging"]]
    regulation_df = pd.DataFrame()
    allostery_samples_1 = infd.posterior.allostery.sel({"experiments": exp_1, "edges": enzyme}).to_dataframe().reset_index().loc[steady_state]
    allostery_samples_2 = infd.posterior.allostery.sel({"experiments": exp_2, "edges": enzyme}).to_dataframe().reset_index().loc[steady_state]
    enz_samples_1 = infd.posterior.conc_enzyme.sel({"experiments": exp_1, "enzymes": enzyme}).to_dataframe().reset_index().loc[steady_state]
    enz_samples_2 = infd.posterior.conc_enzyme.sel({"experiments": exp_2, "enzymes": enzyme}).to_dataframe().reset_index().loc[steady_state]
    saturation_samples_1 = infd.posterior.saturation.sel({"experiments": exp_1, "edges": enzyme}).to_dataframe().reset_index().loc[steady_state]
    saturation_samples_2 = infd.posterior.saturation.sel({"experiments": exp_2, "edges": enzyme}).to_dataframe().reset_index().loc[steady_state]
    reversibility_samples_1 = infd.posterior.reversibility.sel({"experiments": exp_1, "edges": enzyme}).to_dataframe().reset_index().loc[steady_state]
    reversibility_samples_2 = infd.posterior.reversibility.sel({"experiments": exp_2, "edges": enzyme}).to_dataframe().reset_index().loc[steady_state]
    saturation_df = pd.DataFrame.from_dict({"value": np.log(saturation_samples_2["saturation"]) - np.log(saturation_samples_1["saturation"]),
                                            "regulation": "saturation"})
    allostery_df = pd.DataFrame.from_dict({"value": np.log(allostery_samples_2["allostery"]) - np.log(allostery_samples_1["allostery"]),
                                            "regulation": "allostery"})
    reversibility_df = pd.DataFrame.from_dict({"value": np.log(reversibility_samples_2["reversibility"]) - np.log(reversibility_samples_1["reversibility"]),
                                            "regulation": "reversibility"})
    enzyme_df = pd.DataFrame.from_dict({"value": np.log(enz_samples_2["conc_enzyme"]) - np.log(enz_samples_1["conc_enzyme"]),
                                            "regulation": "enzyme"})
    regulation_df = pd.concat([
                               saturation_df, 
                               allostery_df,
                               # reversibility_df,
                               enzyme_df
                               ],
                               )
    return regulation_df

def get_true_regulatory_values(infd):
    enzyme = "GNMT1"
    exp_1 = "dataset_1"
    exp_2 = "dataset_2"

    t_enz_reg = infd.posterior.conc_enzyme.sel({"experiments":[exp_1, exp_2], "enzymes": enzyme}).to_dataframe().reset_index()
    t_allo_reg = infd.posterior.allostery.sel({"experiments":[exp_1, exp_2], "edges": enzyme}).to_dataframe().reset_index()
    t_sat_reg = infd.posterior.saturation.sel({"experiments":[exp_1, exp_2], "edges": enzyme}).to_dataframe().reset_index()
    true_enzyme_reg = np.log(t_enz_reg[t_enz_reg["experiments"]==exp_2]["conc_enzyme"][1] / t_enz_reg[t_enz_reg["experiments"]==exp_1]["conc_enzyme"][0])
    true_allostery_reg = np.log(t_allo_reg[t_allo_reg["experiments"]==exp_2]["allostery"][1] / t_allo_reg[t_allo_reg["experiments"]==exp_1]["allostery"][0])
    true_saturation_reg = np.log(t_sat_reg[t_sat_reg["experiments"]==exp_2]["saturation"][1] / t_sat_reg[t_sat_reg["experiments"]==exp_1]["saturation"][0])

    return [
        true_saturation_reg,
        true_allostery_reg,
        true_enzyme_reg,  
    ]



def get_flux_df(infd):
    reaction = "GNMT"
    exp_1 = "dataset_1"
    exp_2 = "dataset_2"

    flux_df = infd.posterior.flux \
    .sel(
        {
        "experiments": [exp_1, exp_2], 
        "reactions": reaction
        }
        ) \
    .to_dataframe() \
    .reset_index()
    flux_df[flux_df["flux"]==0.0] = np.nan

    return flux_df.dropna()

def plot_figure_5():
    met_13_infd = az.InferenceData.from_netcdf(met_13_path / "infd.nc")
    regulation_df = get_regulation_df(met_13_infd)
    flux_df = get_flux_df(met_13_infd)
    true_values = az.InferenceData.from_netcdf(met_13_true_values / "infd.nc")
    true_df = get_true_regulatory_values(true_values)


    # calling plotting scripts
    cm = 1/2.54  # centimeters in inches
    sns.set_context("paper", rc={"font.size":14,"axes.titlesize":14,"axes.labelsize":14})

    fig = plt.figure(constrained_layout=True, figsize = (21.0*cm, 10.5*cm))
    gs = GridSpec(2, 2, figure=fig)

    ax0 = fig.add_subplot(gs[:, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])

    ax1.set_axis_off()
    ax1.imshow(mpimg.imread('plots/regulatory_diagram.png'), aspect="auto")
    ax0 = plot_regulation(regulation_df, true_df, ax0)
    ax2 = plot_flux(flux_df, ax2)

    plt.figure(dpi = 300)
    fig.tight_layout()
    fig.savefig("plots/figure_5.png")

if __name__ == "__main__":
    plot_figure_5()
