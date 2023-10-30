import builtins
from pathlib import Path
import tempfile
import shutil
import os

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import arviz as az

from maud.loading_maud_inputs import load_maud_input
from maud.getting_idatas import get_idata

from plotting_functions import get_measurements_of_type, get_csvs

MAUD_OUTPUT = Path.cwd() / "methionine_model"/ "results" / "methionine"
PLOT_OUTPUT = Path.cwd() / "methionine_model"/ "plots" / "figure_2"

plt.style.use('ipynb')

colour_scheme = {
    "measurement": "0.6",
    "sample": "#8F2727"
}

def plot_flux_forest(mi, infd, axes):
    """Return forest plot axes of fluxes."""
    flux_measurements = get_measurements_of_type("flux", "yflux", mi.measurements)
    flux_measurements["id"] = flux_measurements["experiments"] + '-' + flux_measurements["flux"]
    flux_samples = infd.posterior["flux"].to_dataframe().reset_index()
    flux_samples["id"] = flux_samples["experiments"] + '-' + flux_samples["reactions"]
    flux_samples = flux_samples.groupby("id") \
                               .quantile([0.05, 0.95]) \
                               .reset_index() \
                               .pivot(index="id", columns="level_1", values="flux")
    sorted_fluxes = list(flux_measurements.sort_values("measurement")["id"])
    for id, sample in flux_samples.reindex(sorted_fluxes).iterrows():
        axes[0,0].vlines(
            x=id,
            ymin=sample[0.05],
            ymax=sample[0.95],
            linewidth=2.5,
            color=colour_scheme["sample"],
        )
    for _, meas in flux_measurements.sort_values("measurement").iterrows():
        axes[0,0].vlines(
            x=meas["id"],
            ymin=meas["measurement"]-1.65*meas["error_scale"],
            ymax=meas["measurement"]+1.65*meas["error_scale"],
            linewidth=2,
            color=colour_scheme["measurement"],
            alpha=0.7,
        )
    axes[0,0].set_xticklabels(axes[0,0].get_xticklabels(),rotation = 45, ha="right")
    axes[0,0].set_ylabel("Flux Measurement [mM/s]")
    # axes[0,0].set_xlabel("flux ID")
    axes[0,0].set(yscale="log")
    return axes

def plot_conc_forest(mi, infd, axes):
    """Return forest plot axes of concentrations."""
    conc_measurements = get_measurements_of_type("conc", "yconc", mi.measurements)
    conc_measurements["id"] = conc_measurements["experiments"] + '-' + conc_measurements["conc"]
    conc_samples = infd.posterior["conc"].to_dataframe().reset_index()
    conc_samples["id"] = conc_samples["experiments"] + '-' + conc_samples["mics"]
    conc_samples = conc_samples.groupby("id") \
                               .quantile([0.05, 0.95]) \
                               .reset_index() \
                               .pivot(index="id", columns="level_1", values="conc")
    sorted_concs = list(conc_measurements.sort_values("measurement")["id"])
    for id, sample in conc_samples.reindex(sorted_concs).iterrows():
        axes[0,1].vlines(
            x=id,
            ymin=sample[0.05],
            ymax=sample[0.95],
            linewidth=2.5,
            color=colour_scheme["sample"],
        )
    for _, meas in conc_measurements.sort_values("measurement").iterrows():
        axes[0,1].vlines(
            x=meas["id"],
            ymin=np.exp(np.log(meas["measurement"])-1.65*meas["error_scale"]),
            ymax=np.exp(np.log(meas["measurement"])+1.65*meas["error_scale"]),
            linewidth=2.5,
            color=colour_scheme["measurement"],
            alpha=0.5,
        )
    axes[0,1].set_xticklabels(axes[0,1].get_xticklabels(),rotation = 45, ha="right")
    axes[0,1].set_ylabel("Concentration Measurement [mM]")
    axes[0,1].set(yscale="log")
    return axes


def plot_pair_corr(mi, infd):
    """Return pair plot axes of correlated parameters."""

    par_1 = "AHC1-hcys-L_c"
    par_2 = "AHC1-ahcys_c"
    par_1_name = "km - "+par_1
    par_2_name = "km - "+par_2
    km_df = infd.posterior["km"].to_dataframe().reset_index()
    par_1_df = km_df.loc[km_df["kms"]==par_1].reset_index(drop=True)
    par_2_df = km_df.loc[km_df["kms"]==par_2].reset_index(drop=True)
    pairplot_df = pd.DataFrame()
    pairplot_df[par_1_name] = par_1_df["km"]
    pairplot_df[par_2_name] = par_2_df["km"]
    pairplot_df = np.log(pairplot_df)

    g0 = sns.jointplot(
        data=pairplot_df, 
        x=par_1_name, 
        y=par_2_name,
        )
    return g0


def plot_pair_uncorr(mi, infd):
    """Return pair plot axes of correlated parameters."""

    par_1 = "MAT3-atp_c"
    par_2 = "MAT3"
    par_1_name = "km - "+par_1
    par_2_name = "kcat - "+par_2
    km_df = infd.posterior["km"].to_dataframe().reset_index()
    par_1_df = km_df.loc[km_df["kms"]==par_1].reset_index(drop=True)
    kcat_df = infd.posterior["kcat"].to_dataframe().reset_index()
    par_2_df = kcat_df.loc[kcat_df["enzymes"]==par_2].reset_index(drop=True)
    pairplot_df = pd.DataFrame()
    pairplot_df[par_1_name] = par_1_df["km"]
    pairplot_df[par_2_name] = par_2_df["kcat"]
    pairplot_df = np.log(pairplot_df)


    g0 = sns.jointplot(
        data=pairplot_df, 
        x=par_1_name, 
        y=par_2_name,
        )
    return g0



def plot_figure():

    # Loading Data
    mi = load_maud_input(data_path=MAUD_OUTPUT / "user_input")
    idata = get_idata([MAUD_OUTPUT / "samples" / "model-20231017125941_1.csv"], mi=mi, mode="sample")

    # calling plotting scripts
    cm = 1/2.54  # centimeters in inches
    sns.set_context("paper", rc={"font.size":14,"axes.titlesize":14,"axes.labelsize":14})

    # Running scripts that generate their own figures
    tempdir = tempfile.mkdtemp(prefix="temp_figures-")

    g0 = plot_pair_corr(mi, infd)
    g0_file = os.path.join(tempdir, "g0.png")
    g0.savefig(g0_file)
    plt.close(g0.fig)

    g1 = plot_pair_uncorr(mi, infd)
    g1_file = os.path.join(tempdir, "g1.png")
    g1.savefig(g1_file)
    plt.close(g1.fig)

    sns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":8})
    samples_legend = mpatches.Patch(color=colour_scheme["sample"], label="5 - 95 % sampling quantile")
    measurements_legend = mpatches.Patch(color=colour_scheme["measurement"], label="true value $\pm$ $45\%$ error")

    fig, axes = plt.subplots(2, 2, figsize = (21.0*cm, 21.0*cm))
    axes = plot_flux_forest(mi, infd, axes)
    axes = plot_conc_forest(mi, infd, axes)
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles=[samples_legend, measurements_legend], loc='upper center',
                      ncol=2, borderaxespad=0.)
    axes[1,0].set_axis_off()
    axes[1,0].imshow(mpimg.imread(g0_file), aspect="auto")
    axes[1,1].set_axis_off()
    axes[1,1].imshow(mpimg.imread(g1_file), aspect="auto")

    plt.figure(dpi = 300)
    fig.tight_layout()
    fig.subplots_adjust(top=0.97)
    fig.savefig("plots/figure_3.png")
    shutil.rmtree(tempdir)

    

if __name__ == "__main__":
    plot_figure()


    
