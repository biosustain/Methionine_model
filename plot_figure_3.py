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

"""Code for plotting posterior distribution."""
import os
from pathlib import Path
import arviz as az
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
import matplotlib.transforms as transforms
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from maud.loading_maud_inputs import load_maud_input
from maud.getting_idatas import get_idata

plt.style.use('ipynb')


HELP_MSG = """
This script plots boxplots of enzyme concentrations,
metabolite concentrations, and fluxes of out-of-sample
predictions. Measurements are also included in these
plots.
"""

MAUD_OUTPUT_MISSING_AHCYS = Path.cwd() / "results" / "methionine_missing_ahcys"
MAUD_OUTPUT = Path.cwd() / "results" / "methionine"
PLOT_OUTPUT = Path.cwd() / "figures" / "figure_3.png"

colour_scheme = {
    "complete measurement set": "#001219",
    "missing ahcys measurement": "#005F73"
}

OFFSET = 0.05
cm = 1/2.54  # centimeters in inches

def plot_residual_state_variables(residual_df, ax):
    offset = lambda p: transforms.ScaledTranslation(p/72.,0, plt.gcf().dpi_scale_trans)
    for id, sample in residual_df.loc[residual_df["model"]=="missing ahcys measurement"].reset_index().iterrows():
        ax.vlines(
            x=id,
            ymin=sample[0.05],
            ymax=sample[0.95],
            linewidth=2.5,
            color=colour_scheme["missing ahcys measurement"],
            transform=ax.transData+offset(+1),
            alpha=0.7,
        )
    for id, sample in residual_df.loc[residual_df["model"]=="complete measurement set"].reset_index().iterrows():
        ax.vlines(
            x=id,
            ymin=sample[0.05],
            ymax=sample[0.95],
            linewidth=2.5,
            color=colour_scheme["complete measurement set"],
            transform=ax.transData+offset(-1),
            alpha=0.7
        )
    ax.axhline(y=0, xmin=0, xmax=len(residual_df)/2, color="red")
    ax.set_ylabel("log residual plot of balanced metabolites")
    ax.get_xaxis().set_visible(False)
    return ax

def plot_box(df, ax):
    ax.vlines(
        x="complete measurement set",
        ymin=df.loc[df["model"]=="complete measurement set"]["lower_5_lp"],
        ymax=df.loc[df["model"]=="complete measurement set"]["upper_95_lp"],
        linewidth=25,
        color=colour_scheme["complete measurement set"],
        alpha=0.7,
    )
    ax.vlines(
        x="missing ahcys_c measurement",
        ymin=df.loc[df["model"]=="missing ahcys measurement"]["lower_5_lp"],
        ymax=df.loc[df["model"]=="missing ahcys measurement"]["upper_95_lp"],
        linewidth=25,
        color=colour_scheme["missing ahcys measurement"],
        alpha=0.7,
    )
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylabel("log likelihood of flux measurements")
    ax.get_xaxis().set_visible(False)
    return ax

def ll_normal(value, mean, std):
    return -np.log(2*np.pi)/2-np.log(std**2)/2-(value-mean)**2/(2*std**2)

def plot_parameter_hist(parameter_df, ax):
    sns.histplot(parameter_df, x="ln_km", bins=30, hue="model", stat="probability", ax=ax, legend=False, alpha=0.7)
    ax.set_xlabel("ln $K_m^{ACH1, ahcys_c}$")
    ax.axvline(x=-10.6713582793, color="red")
    return ax

def plot_figure():
    mi = load_maud_input(data_path=MAUD_OUTPUT / "user_input")
    mi_missing = load_maud_input(data_path=MAUD_OUTPUT_MISSING_AHCYS / "user_input")
    csvs_methionine = [MAUD_OUTPUT/ "samples" / file for file in os.listdir(MAUD_OUTPUT/ "samples") if ".csv" in file]
    infd_methionine = get_idata(csvs_methionine, mi=mi, mode="train")
    csvs_methionine_missing = [MAUD_OUTPUT_MISSING_AHCYS/ "samples" / file for file in os.listdir(MAUD_OUTPUT_MISSING_AHCYS/ "samples") if ".csv" in file]
    infd_methionine_missing = get_idata(csvs_methionine_missing, mi=mi_missing, mode="train")

    flux_measurements = [{
        "flux": meas.value, 
        "uncertainty": meas.error_scale,
        "reaction": meas.reaction, 
        "experiment": meas.experiment
    } 
        for exp in mi.experiments 
        for meas in exp.measurements 
        if (meas.target_type=="flux")
    ]
    flux_meas_df = pd.DataFrame.from_records(flux_measurements)
    flux_meas_df["id"] = flux_meas_df["experiment"]+flux_meas_df["reaction"]
    flux_meas_df = flux_meas_df.set_index("id")

    conc_measurements = [{
        "conc": meas.value, 
        "uncertainty": meas.error_scale,
        "metabolite": meas.metabolite, 
        "compartment": meas.compartment,
        "experiment": meas.experiment
    } 
        for exp in mi.experiments 
        for meas in exp.measurements 
        if (meas.target_type=="mic")
    ]
    conc_meas_df = pd.DataFrame.from_records(conc_measurements)
    conc_meas_df["id"] = conc_meas_df["experiment"]+conc_meas_df["metabolite"]+"_"+conc_meas_df["compartment"]
    conc_meas_df = conc_meas_df.set_index("id")
        
    flux_lp_df = pd.DataFrame(columns=["lower_5_lp", "upper_95_lp", "model"])
    flux_methionine = infd_methionine.posterior.flux_train.to_dataframe().dropna().droplevel(1).reset_index()
    flux_methionine_missing = infd_methionine_missing.posterior.flux_train.to_dataframe().dropna().droplevel(1).reset_index()

    lp_flux_methionine = []
    for _, row in flux_methionine.iterrows():
        id = row["experiments"]+row["reactions"]
        if id in flux_meas_df.index:
            lp_flux_methionine.append(ll_normal(row["flux_train"],flux_meas_df.loc[id]["flux"],flux_meas_df.loc[id]["uncertainty"]))
    flux_lp_df.loc[0,"lower_5_lp"] = np.quantile(lp_flux_methionine, 0.05)
    flux_lp_df.loc[0,"upper_95_lp"] = np.quantile(lp_flux_methionine, 0.95)
    flux_lp_df.loc[0, "model"] = "complete measurement set"
    lp_flux_methionine_missing = []
    for _, row in flux_methionine_missing.iterrows():
        id = row["experiments"]+row["reactions"]
        if id in flux_meas_df.index:
            lp_flux_methionine_missing.append(ll_normal(row["flux_train"],flux_meas_df.loc[id]["flux"],flux_meas_df.loc[id]["uncertainty"]))

    flux_lp_df.loc[1,"lower_5_lp"] = np.quantile(lp_flux_methionine_missing, 0.05)
    flux_lp_df.loc[1,"upper_95_lp"] = np.quantile(lp_flux_methionine_missing, 0.95)
    flux_lp_df.loc[1, "model"] = "missing ahcys measurement"

    balanced_mics = [mic.id for mic in mi.kinetic_model.mics if mic.balanced == True]
    conc_methionine = infd_methionine.posterior.conc_train.sel(mics=balanced_mics).to_dataframe().dropna().droplevel(1).reset_index()
    conc_methionine["model"] = "complete measurement set"
    conc_methionine_missing = infd_methionine_missing.posterior.conc_train.sel(mics=balanced_mics).to_dataframe().dropna().droplevel(1).reset_index()
    conc_methionine_missing["model"] = "missing ahcys measurement"

    met_AHC1_ahcys = infd_methionine.posterior.km.sel(kms=["AHC1_ahcys_c"]).to_dataframe().reset_index()
    met_AHC1_ahcys["model"] = "complete measurement set"
    met_missing_ahcys_AHC1_ahcys = infd_methionine_missing.posterior.km.sel(kms=["AHC1_ahcys_c"]).to_dataframe().reset_index()
    met_missing_ahcys_AHC1_ahcys["model"] = "missing ahcys measurement"
    AHC1_ahcys_df = pd.concat([met_AHC1_ahcys, met_missing_ahcys_AHC1_ahcys]).reset_index()
    AHC1_ahcys_df["ln_km"] = AHC1_ahcys_df["km"].apply(lambda x: np.log(x))

    for i, row in conc_methionine.iterrows():
        id = row["experiments"]+row["mics"]
        conc_methionine.loc[i, "log_residual"] = np.log(row["conc_train"]) - np.log(conc_meas_df.loc[id]["conc"])
    conc_methionine = conc_methionine.loc[conc_methionine["experiments"]!="validation"]
   
    for i, row in conc_methionine_missing.iterrows():
        id = row["experiments"]+row["mics"]
        conc_methionine_missing.loc[i, "log_residual"] = np.log(row["conc_train"]) - np.log(conc_meas_df.loc[id]["conc"])
    conc_methionine_missing = conc_methionine_missing.loc[conc_methionine_missing["experiments"]!="validation"]

    conc_validation = pd.concat([conc_methionine, conc_methionine_missing], axis=0)
    conc_validation = conc_validation.groupby(["experiments", "mics","model"])\
                                     .quantile([0.05, 0.95])\
                                     .reset_index()\
                                     .pivot(index=["mics", "experiments", "model"], columns="level_3", values="log_residual")\
                                     .reset_index()

    # calling plotting scripts
    cm = 1/2.54  # centimeters in inches
    sns.set_context("paper", rc={"font.size":14,"axes.titlesize":14,"axes.labelsize":14})

    fig = plt.figure(constrained_layout=True, figsize = (21.0*cm, 21.0*cm))
    gs = GridSpec(2, 2, figure=fig)

    samples_legend = mpatches.Patch(color=colour_scheme["complete measurement set"], label="complete measurement set")
    measurements_legend = mpatches.Patch(color=colour_scheme["missing ahcys measurement"], label="missing ahcys measurement")
    true_value_legend = Line2D([0], [0], color='r', lw=3, label='true value')


    fig.legend(handles=[samples_legend, measurements_legend, true_value_legend], 
               loc='upper center',
               ncol=3, 
               borderaxespad=0.,
               fontsize=12,)
    fig.suptitle(' ')

    ax0 = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])

    ax0 = plot_residual_state_variables(conc_validation, ax0)
    ax1 = plot_parameter_hist(AHC1_ahcys_df, ax1)
    ax2 = plot_box(flux_lp_df, ax2)
    fig.savefig(PLOT_OUTPUT)

if __name__ == "__main__":
    plot_figure()
