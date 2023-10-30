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
import argparse
from pathlib import Path
import tempfile
import shutil
import seaborn as sns

import arviz as az
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
import matplotlib.transforms as transforms

from maud.loading_maud_inputs import load_maud_input
from maud.getting_idatas import get_idata


plt.style.use('ipynb')


HELP_MSG = """
This script plots boxplots of enzyme concentrations,
metabolite concentrations, and fluxes of out-of-sample
predictions. Measurements are also included in these
plots.
"""

MAUD_OUTPUT_TRAINING = Path.cwd() / "results" / "example_ode"
MAUD_OUTPUT_VALIDATION = MAUD_OUTPUT_TRAINING / "validation"
MAUD_OUTPUT_LAPLACE_TRAINING = Path.cwd() / "results" / "example_ode_laplace"
MAUD_OUTPUT_VALIDATION = MAUD_OUTPUT_LAPLACE_TRAINING / "validation"
PLOT_OUTPUT = Path.cwd() / "figures" / "figure_4.png"

colour_scheme = {
    "measurement": "0.2",
    "laplace": "#005F73",
    "posterior": "#001219"
}

OFFSET = 0.05
cm = 1/2.54  # centimeters in inches

def plot_lp_testing(axes, lp):
    color_map = {
        "posterior": "#001219",
        "laplace_approximation": "#005F73"
    }
    sns.histplot(
    data=lp, x="lp", hue="sample",
    log_scale=False, element="step", fill=False,
    cumulative=True, stat="probability", marker="",
    ax=axes[0,0],
    palette=color_map,
    )
    axes[0,0].get_legend().remove()
    axes[0,0] = axes[0,0].twinx()
    sns.histplot(
    data=lp, x="lp", hue="sample",
    cumulative=False,
    ax=axes[0,0],
    palette=color_map,
    )
    axes[0,0].set(xlim=(-500, -0))
    axes[0,0].get_legend().remove()
    return axes

def lpdf_normal(flux, meas, scale):
    """calculates the log probability density of a normal distribution"""
    return np.log((1.0/(scale*np.sqrt(np.pi*2.0)))*np.exp(-0.5*((flux-meas)/scale)**2.0))


def plot_pair_corr(df, par_1, par_2, axes):
    """Return pair plot axes of correlated parameters."""
    color_map = {
        "posterior": "#001219",
        "laplace_approximation": "#005F73"
    }
    plt.scatter(
        x=df[par_1], 
        y=df[par_2],
        c=df["sample"].apply(lambda x: color_map[x]),
        alpha = 0.7,
        axes=axes[1,1]
        )
    axes[1,1].set(xlabel=r"$Km_{r1}^A$")
    axes[1,1].set(ylabel=r"$kcat_{r1}$")
    return axes

def plot_llik_flux_validation(df, axes):
    """Return forest plot axes of fluxes."""
    count = 0
    for rxn in df.reaction.unique():
        for exp in df.experiment.unique():
            axes[1,0].vlines(
                x=count+0.1,
                ymin= df.loc[
                    (df["sample"]=="posterior") &
                    (df["reaction"]==rxn) &
                    (df["experiment"]==exp)]["lower_5_lp"],
                ymax=df.loc[
                    (df["sample"]=="posterior") &
                    (df["reaction"]==rxn) &
                    (df["experiment"]==exp)]["upper_95_lp"],
                linewidth=2.5,
                color=colour_scheme["posterior"],
                alpha=0.7,
            )
            axes[1,0].vlines(
                x=count- 0.1,
                ymin=df.loc[
                    (df["sample"]=="laplace_approximation") &
                    (df["reaction"]==rxn) &
                    (df["experiment"]==exp)]["lower_5_lp"],
                ymax=df.loc[
                    (df["sample"]=="laplace_approximation") &
                    (df["reaction"]==rxn) &
                    (df["experiment"]==exp)]["upper_95_lp"],
                linewidth=2.5,
                color=colour_scheme["laplace"],
                alpha=0.7,
            )
            count += 1
    axes[1,0].set_ylabel("Log likelihood of predicted flux values")
    axes[1,0].get_xaxis().set_visible(False)
    return axes

def plot_llik_flux_training(log_lik_df, axes):
    for _, sample in log_lik_df.loc[log_lik_df["sample"]=="posterior"].iterrows():
        axes[0,1].vlines(
            x=sample["id"]+0.2,
            ymin=sample["lower_5"],
            ymax=sample["upper_95"],
            linewidth=2.5,
            color=colour_scheme["posterior"],
            alpha=0.7,
        )
    for id, sample in log_lik_df.loc[log_lik_df["sample"]=="laplace_approximation"].iterrows():
        axes[0,1].vlines(
            x=sample["id"]-0.2,
            ymin=sample["lower_5"],
            ymax=sample["upper_95"],
            linewidth=2.5,
            color=colour_scheme["laplace"],
            alpha=0.7,
        )
    axes[0,1].set_ylabel("Log likelihood of sampled fluxes")
    axes[0,1].get_xaxis().set_visible(False)
    return axes


def plot_figure():
    
    mi = load_maud_input(data_path=MAUD_OUTPUT_TRAINING / "user_input")
    csvs_training = [MAUD_OUTPUT_TRAINING / "samples" / file for file in os.listdir(MAUD_OUTPUT_TRAINING / "samples") if ".csv" in file]
    csvs_validation= [MAUD_OUTPUT_TRAINING / "validation" / file for file in os.listdir(MAUD_OUTPUT_TRAINING / "validation") if ".csv" in file]
    csvs_laplace_training = [MAUD_OUTPUT_LAPLACE_TRAINING /"samples" / file for file in os.listdir(MAUD_OUTPUT_LAPLACE_TRAINING / "samples") if ".csv" in file]   
    csvs_laplace_validation = [MAUD_OUTPUT_LAPLACE_TRAINING / "validation" / file for file in os.listdir(MAUD_OUTPUT_LAPLACE_TRAINING / "validation") if ".csv" in file]       
    idata_train = get_idata(csvs_training, mi=mi, mode="sample")
    idata_validation = get_idata(csvs_validation, mi=mi, mode="test")
    idata_laplace_train= get_idata(csvs_laplace_training, mi=mi, mode="sample")
    idata_laplace_validation = get_idata(csvs_laplace_validation, mi=mi, mode="test")
    lp = [{
        "sample": sample,
        "lp": lp}
        for idata, sample, term in [
            [idata_train, "posterior", "lp"], 
            [idata_laplace_train,"laplace_approximation", "log_p"]]
        for lp in idata.sample_stats[term].to_series().values
    ]
    lp_df = pd.DataFrame.from_records(lp)

    llik_flux = [{
        "sample": sample,
        "id": flux,
        "lower_5": idata.posterior.llik_flux_train.sel(llik_flux_train_dim_0=[flux]).to_series().quantile(0.05),
        "upper_95": idata.posterior.llik_flux_train.sel(llik_flux_train_dim_0=[flux]).to_series().quantile(0.95)
    }
        for idata, sample in [
            [idata_train, "posterior"], 
            [idata_laplace_train,"laplace_approximation"]]
        for flux in idata.posterior.coords["llik_flux_train_dim_0"].values
    ]
    llik_flux_df = pd.DataFrame.from_records(llik_flux)

    measurements = [{
        "flux": meas.value, 
        "uncertainty": meas.error_scale,
        "reaction": meas.reaction, 
        "experiment": meas.experiment
    } 
        for exp in mi.experiments 
        for meas in exp.measurements 
        if meas.target_type=="flux"
    ]
    meas_df = pd.DataFrame.from_records(measurements)
    meas_df["id"] = meas_df["experiment"]+meas_df["reaction"]
    meas_df = meas_df.set_index("id")

    
    flux_interval = [{
        "sample": sample,
        "reaction": rxn,
        "experiment": exp,
        "lower_5_lp": 
            idata.posterior.flux_test.sel(reactions=[rxn], experiments=[exp]).to_series().apply(
                lambda x: lpdf_normal(x, meas_df.loc[exp+rxn].flux, meas_df.loc[exp+rxn].uncertainty)).quantile(0.05),
        "upper_95_lp": 
            idata.posterior.flux_test.sel(reactions=[rxn], experiments=[exp]).to_series().apply(
                lambda x: lpdf_normal(x, meas_df.loc[exp+rxn].flux, meas_df.loc[exp+rxn].uncertainty)).quantile(0.95)
    }
        for idata, sample in [
            [idata_validation, "posterior"], 
            [idata_laplace_validation,"laplace_approximation"]]
        for rxn in idata.posterior.coords["reactions"].values
        for exp in idata.posterior.coords["experiments"].values
        if exp+rxn in meas_df.index
    ]
    flux_interval_df = pd.DataFrame.from_records(flux_interval)

    print(scipy.stats.ks_2samp(
          lp_df[lp_df["sample"]=="posterior"]["lp"],
          lp_df[lp_df["sample"]=="laplace_approximation"]["lp"],
          alternative="two-sided"))

    par_1 = "r1_A_c"
    par_2 = "r1"
    km_posterior = idata_train.posterior.km.sel(kms=[par_1]).to_dataframe().unstack().km.droplevel(1).reset_index()
    kcat_posterior = idata_train.posterior.kcat.sel(enzymes=[par_2]).to_dataframe().unstack().kcat.droplevel(1).reset_index()    
    posterior_pars = pd.concat([km_posterior, kcat_posterior], axis=1)
    posterior_pars["sample"]="posterior"
    km_laplace= idata_laplace_train.posterior.km.sel(kms=[par_1]).to_dataframe().unstack().km.droplevel(1).reset_index()
    kcat_laplace= idata_laplace_train.posterior.kcat.sel(enzymes=[par_2]).to_dataframe().unstack().kcat.droplevel(1).reset_index()    
    laplace_pars = pd.concat([km_laplace, kcat_laplace], axis=1)
    laplace_pars["sample"]="laplace_approximation"
    par_df = pd.concat([laplace_pars, posterior_pars])
    # calling plotting scripts

    cm = 1/2.54  # centimeters in inches
    fig, axes = plt.subplots(2, 2, figsize = (21.0*cm, 21.0*cm))
    axes = plot_pair_corr(par_df, par_1, par_2, axes)
    axes = plot_lp_testing(axes, lp_df)
    axes = plot_llik_flux_validation(flux_interval_df, axes)
    axes = plot_llik_flux_training(llik_flux_df, axes)
    true_posterior_legend = mpatches.Patch(color=colour_scheme["posterior"],
                                           alpha=0.7,
                                           label="posterior")
    laplace_legend = mpatches.Patch(color=colour_scheme["laplace"],
                                alpha=0.7,
                                label="laplace_approximation")
    fig.legend(handles=[true_posterior_legend, laplace_legend], 
               loc='upper center',
               ncol=2, 
               bbox_to_anchor=(0.5,1.0),
               fontsize=12,)
    plt.figure(dpi = 300)
    fig.tight_layout()
    fig.subplots_adjust(top=0.93)
    fig.savefig(PLOT_OUTPUT)

if __name__ == "__main__":
    plot_figure()
