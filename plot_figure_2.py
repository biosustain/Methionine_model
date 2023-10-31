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

from maud.loading_maud_inputs import load_maud_input
from maud.getting_idatas import get_idata


MAUD_OUTPUT = Path.cwd() / "results" / "methionine"
PLOT_OUTPUT = Path.cwd() / "figures" / "posterior.png"

plt.style.use("ipynb")

colour_scheme = {"measurement": "#8F2727", "sample": "black"}


def plot_flux_forest(mi, idata_methionine, axes, flux_meas_df):
    """Return forest plot axes of fluxes."""
    flux_samples = idata_methionine.posterior["flux_train"].to_dataframe().reset_index()
    flux_samples["id"] = flux_samples["experiments"] + flux_samples["reactions"]
    flux_samples = (
        flux_samples[["id", "flux_train"]]
        .groupby("id")
        .quantile([0.05, 0.95])
        .reset_index()
        .pivot(index="id", columns="level_1", values="flux_train")
    )
    sorted_fluxes = list(flux_meas_df.sort_values("flux").index)
    for id, sample in flux_samples.reindex(sorted_fluxes).iterrows():
        axes[0, 0].vlines(
            x=id,
            ymin=sample[0.05],
            ymax=sample[0.95],
            linewidth=2.5,
            color=colour_scheme["sample"],
            zorder=0,
        )
    f = flux_meas_df.reindex(sorted_fluxes).reset_index()
    axes[0, 0].scatter(
        f.index,
        f["flux"],
        marker="o",
        color=colour_scheme["measurement"],
    )
    axes[0, 0].set_xticklabels([])
    axes[0, 0].set_ylabel("Flux [mM/s]")
    # axes[0,0].set_xlabel("flux ID")
    axes[0, 0].set(yscale="log")
    return axes


def plot_conc_forest(mi, idata_methionine, axes, conc_meas_df):
    """Return forest plot axes of concentrations."""
    conc_samples = idata_methionine.posterior["conc_train"].to_dataframe().reset_index()
    conc_samples["id"] = conc_samples["experiments"] + conc_samples["mics"]
    conc_samples = (
        conc_samples[["id", "conc_train"]]
        .groupby("id")
        .quantile([0.05, 0.95])
        .reset_index()
        .pivot(index="id", columns="level_1", values="conc_train")
    )
    sorted_concs = list(conc_meas_df.sort_values("conc").index)
    for id, sample in conc_samples.reindex(sorted_concs).iterrows():
        axes[0, 1].vlines(
            x=id,
            ymin=sample[0.05],
            ymax=sample[0.95],
            linewidth=2.5,
            color=colour_scheme["sample"],
            zorder=0,
        )

    c = conc_meas_df.reindex(sorted_concs).reset_index()
    axes[0, 1].scatter(
        c.index,
        c["conc"],
        marker="o",
        color=colour_scheme["measurement"],
    )
    axes[0, 1].set_xticklabels([])
    axes[0, 1].set_ylabel("Concentration [mM]")
    axes[0, 1].set(yscale="log")
    return axes


def plot_pair_corr(mi, infd):
    """Return pair plot axes of correlated parameters."""

    par_1 = "AHC1_hcys-L_c"
    par_2 = "AHC1_ahcys_c"
    par_1_name = "km - " + par_1
    par_2_name = "km - " + par_2
    km_df = infd.posterior["km"].to_dataframe().reset_index()
    par_1_df = km_df.loc[km_df["kms"] == par_1].reset_index(drop=True)
    par_2_df = km_df.loc[km_df["kms"] == par_2].reset_index(drop=True)
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

    par_1 = "MAT3_atp_c"
    par_2 = "MAT3"
    par_1_name = "km - " + par_1
    par_2_name = "kcat - " + par_2
    km_df = infd.posterior["km"].to_dataframe().reset_index()
    par_1_df = km_df.loc[km_df["kms"] == par_1].reset_index(drop=True)
    kcat_df = infd.posterior["kcat"].to_dataframe().reset_index()
    par_2_df = kcat_df.loc[kcat_df["enzymes"] == par_2].reset_index(drop=True)
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
    csvs_methionine = [
        MAUD_OUTPUT / "samples" / file
        for file in os.listdir(MAUD_OUTPUT / "samples")
        if ".csv" in file
    ]
    idata_methionine = get_idata(csvs_methionine, mi=mi, mode="train")

    flux_measurements = [
        {
            "flux": meas.value,
            "uncertainty": meas.error_scale,
            "reaction": meas.reaction,
            "experiment": meas.experiment,
        }
        for exp in mi.experiments
        for meas in exp.measurements
        if (meas.target_type == "flux") & (exp.is_train == True)
    ]
    flux_meas_df = pd.DataFrame.from_records(flux_measurements)
    flux_meas_df["id"] = flux_meas_df["experiment"] + flux_meas_df["reaction"]
    flux_meas_df = flux_meas_df.set_index("id")

    conc_measurements = [
        {
            "conc": meas.value,
            "uncertainty": meas.error_scale,
            "metabolite": meas.metabolite,
            "compartment": meas.compartment,
            "experiment": meas.experiment,
        }
        for exp in mi.experiments
        for meas in exp.measurements
        if (meas.target_type == "mic") & (exp.is_train == True)
    ]
    conc_meas_df = pd.DataFrame.from_records(conc_measurements)
    conc_meas_df["id"] = (
        conc_meas_df["experiment"]
        + conc_meas_df["metabolite"]
        + "_"
        + conc_meas_df["compartment"]
    )
    conc_meas_df = conc_meas_df.set_index("id")

    # calling plotting scripts
    cm = 1 / 2.54  # centimeters in inches
    sns.set_context(
        "paper", rc={"font.size": 14, "axes.titlesize": 14, "axes.labelsize": 14}
    )

    # Running scripts that generate their own figures
    tempdir = tempfile.mkdtemp(prefix="temp_figures-")

    g0 = plot_pair_corr(mi, idata_methionine)
    g0_file = os.path.join(tempdir, "g0.png")
    g0.savefig(g0_file)
    plt.close(g0.fig)

    g1 = plot_pair_uncorr(mi, idata_methionine)
    g1_file = os.path.join(tempdir, "g1.png")
    g1.savefig(g1_file)
    plt.close(g1.fig)

    sns.set_context(
        "paper", rc={"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8}
    )
    samples_legend = mpatches.Patch(
        color=colour_scheme["sample"], label="5% - 95% sampling quantile"
    )
    measurements_legend = mpatches.Patch(
        color=colour_scheme["measurement"], label="Measured value"
    )

    fig, axes = plt.subplots(2, 2, figsize=(21.0 * cm, 21.0 * cm))
    axes = plot_flux_forest(mi, idata_methionine, axes, flux_meas_df)
    axes = plot_conc_forest(mi, idata_methionine, axes, conc_meas_df)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles=[samples_legend, measurements_legend],
        loc="upper center",
        ncol=2,
        borderaxespad=0.0,
    )
    axes[1, 0].set_axis_off()
    axes[1, 0].imshow(mpimg.imread(g0_file), aspect="auto")
    axes[1, 1].set_axis_off()
    axes[1, 1].imshow(mpimg.imread(g1_file), aspect="auto")

    plt.figure(dpi=300)
    fig.tight_layout()
    fig.subplots_adjust(top=0.97)
    fig.savefig(PLOT_OUTPUT)
    shutil.rmtree(tempdir)


if __name__ == "__main__":
    plot_figure()
