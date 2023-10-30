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
from pathlib import Path
import tempfile

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

from maud import io

plt.style.use('ipynb')


HELP_MSG = """
This script plots boxplots of enzyme concentrations,
metabolite concentrations, and fluxes of out-of-sample
predictions. Measurements are also included in these
plots.
"""

MAUD_OUTPUT_MISSING_AHCYS = Path.cwd() / "results" / "methionine_cycle_6_missing_ahcys_validation"
MAUD_OUTPUT = Path.cwd() / "results" / "methionine_cycle_6_validation"
PLOT_OUTPUT = Path.cwd() / "plots" / "figure_4"

colour_scheme = {
    "complete_measurement_set": "#001219",
    "missing_ahcys_measurement": "#005F73"
}

OFFSET = 0.05
cm = 1/2.54  # centimeters in inches

def plot_residual_state_variables(residual_df, ax):

    offset = lambda p: transforms.ScaledTranslation(p/72.,0, plt.gcf().dpi_scale_trans)
    for id, sample in residual_df.loc[residual_df["model"]=="missing ahcys_c measurement"].reset_index().iterrows():
        ax.vlines(
            x=id,
            ymin=sample[0.05],
            ymax=sample[0.95],
            linewidth=2.5,
            color=colour_scheme["missing_ahcys_measurement"],
            transform=ax.transData+offset(+1),
            alpha=0.7,
        )
    for id, sample in residual_df.loc[residual_df["model"]=="complete measurement set"].reset_index().iterrows():
        ax.vlines(
            x=id,
            ymin=sample[0.05],
            ymax=sample[0.95],
            linewidth=2.5,
            color=colour_scheme["complete_measurement_set"],
            transform=ax.transData+offset(-1),
            alpha=0.7
        )
    ax.axhline(y=0, color="red")
    ax.set_ylabel("log residual plot of balanced metabolites")
    ax.get_xaxis().set_visible(False)

    return ax

def plot_box(df, ax):
    ax.vlines(
        x="complete measurement set",
        ymin=df.loc["complete measurement set"]["q_0.05"],
        ymax=df.loc["complete measurement set"]["q_0.95"],
        linewidth=25,
        color=colour_scheme["complete_measurement_set"],
        alpha=0.7,
    )
    ax.vlines(
        x="missing ahcys_c measurement",
        ymin=df.loc["missing ahcys_c measurement"]["q_0.05"],
        ymax=df.loc["missing ahcys_c measurement"]["q_0.95"],
        linewidth=25,
        color=colour_scheme["missing_ahcys_measurement"],
        alpha=0.7,
    )
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylabel("log likelihood of flux measurements")
    ax.get_xaxis().set_visible(False)
    return ax

def ll_normal(value, mean, std):
    return -np.log(2*np.pi)/2-np.log(std**2)/2-(value-mean)**2/(2*std**2)

def extract_log_likelihood_flux(df, mi):
    measurements = mi.measurements \
                     .yflux \
                     .reset_index() \
                     .rename(columns={"target_id": "reactions", "experiment_id": "experiments"})
    df = df.merge(measurements, left_on=["experiments", "reactions"], right_on=["experiments", "reactions"])[["model", "experiments", "reactions", "flux", "measurement", "error_scale"]]
    df["ll"] = df.apply(lambda x: ll_normal(x["flux"], x["measurement"], x["error_scale"]), axis=1)

    return df


def filter_dataframe(df):
    df_filtered = pd.DataFrame()
    df_quantiles = df.groupby(["model"])["ll"] \
                     .quantile([0.05, 0.95]) \
                     .unstack() \
                     .add_prefix("q_") \
                     .reset_index() \
                     .set_index("model")

    for model in df_quantiles.index:
        df_filtered = pd.concat([df_filtered, df[(df["model"] == model) \
                                & (df["ll"] > df_quantiles.loc[model]["q_0.05"]) \
                                & (df["ll"] < df_quantiles.loc[model]["q_0.95"])]])
    return df_quantiles

def plot_parameter_hist(parameter_df, ax):
    sns.histplot(parameter_df, x="ln_km", bins=30, hue="model", stat="probability", ax=ax, legend=False, alpha=0.7)
    ax.set_xlabel("ln $K_m^{ACH1, ahcys_c}$")
    ax.axvline(x=-10.6713582793, color="red")
    return ax

def plot_figure():
    mi_6 = io.load_maud_input(data_path=MAUD_OUTPUT / "user_input", mode="predict")
    mi_6_missing_ahcys = io.load_maud_input(data_path=MAUD_OUTPUT_MISSING_AHCYS / "user_input", mode="predict")
    met_6_infd = az.InferenceData.from_netcdf(MAUD_OUTPUT / "trained_samples" / "infd.nc")
    met_6_missing_ahcys_infd = az.InferenceData.from_netcdf(MAUD_OUTPUT_MISSING_AHCYS / "trained_samples" / "infd.nc")
    met_6_predict = pd.read_csv(MAUD_OUTPUT / "oos_samples" / "flux.csv")
    met_6_predict = met_6_predict[met_6_predict["flux"] != 0]
    met_6_predict["model"] = "complete measurement set"
    met_6_missing_ahcys_predict = pd.read_csv(MAUD_OUTPUT_MISSING_AHCYS / "oos_samples" / "flux.csv")
    met_6_missing_ahcys_predict = met_6_missing_ahcys_predict[met_6_missing_ahcys_predict["flux"] != 0]
    met_6_missing_ahcys_predict["model"] = "missing ahcys_c measurement"
    model_ll = pd.concat([met_6_predict, met_6_missing_ahcys_predict])

    target_reaction = "AHC"

    balanced_mics = [mic.id for mic in mi_6.kinetic_model.mics if mic.balanced == True]
    concs_validation = pd.DataFrame()
    concs_validation_6 = pd.read_csv(MAUD_OUTPUT / "oos_samples" / "conc.csv")
    concs_validation_6 = concs_validation_6.loc[concs_validation_6["conc"] != 0]
    concs_validation_6 = concs_validation_6.loc[concs_validation_6["mics"].isin(balanced_mics)][["experiments", "mics", "conc"]]
    concs_validation_6["model"] = "complete measurement set"
    concs_validation_6_missing_ahcys = pd.read_csv(MAUD_OUTPUT_MISSING_AHCYS / "oos_samples" / "conc.csv")
    concs_validation_6_missing_ahcys = concs_validation_6_missing_ahcys.loc[concs_validation_6_missing_ahcys["conc"] != 0]
    concs_validation_6_missing_ahcys = concs_validation_6_missing_ahcys.loc[concs_validation_6_missing_ahcys["mics"].isin(balanced_mics)][["experiments", "mics", "conc"]]
    concs_validation_6_missing_ahcys["model"] = "missing ahcys_c measurement"
    concs_validation = pd.concat([concs_validation_6, concs_validation_6_missing_ahcys])
    concs_validation = concs_validation.groupby(["mics", "experiments", "model"]) \
                                       .quantile([0.05, 0.95]) \
                                       .reset_index() \
                                       .pivot(index=["mics", "model", "experiments"], columns="level_3", values="conc") \
                                       .reset_index()

    met_6_flux_samples = met_6_infd.posterior \
                                  .flux \
                                  .sel(reactions=[target_reaction]) \
                                  .to_dataframe() \
                                  .reset_index()[["experiments", "reactions", "flux"]] \
                                  .groupby(["experiments", "reactions"]) \
                                  .quantile([0.05, 0.95]) \
                                  .unstack()["flux"] \
                                  .add_prefix("q_") \
                                  .reset_index()
    met_6_missing_ahcys_flux_samples = met_6_missing_ahcys_infd.posterior \
                                                              .flux \
                                                              .sel(reactions=[target_reaction]) \
                                                              .to_dataframe() \
                                                              .reset_index()[["experiments", "reactions", "flux"]] \
                                                              .groupby(["experiments", "reactions"]) \
                                                              .quantile([0.05, 0.95]) \
                                                              .unstack()["flux"] \
                                                              .add_prefix("q_") \
                                                              .reset_index()

    df_ll = extract_log_likelihood_flux(model_ll, mi_6)
    df_ll = filter_dataframe(df_ll)
    met_6_AHC1_ahcys = met_6_infd.posterior.km.sel(kms=["AHC1-ahcys_c"]).to_dataframe().reset_index()
    met_6_AHC1_ahcys["model"] = "complete measurement set"
    met_6_missing_ahcys_AHC1_ahcys = met_6_missing_ahcys_infd.posterior.km.sel(kms=["AHC1-ahcys_c"]).to_dataframe().reset_index()
    met_6_missing_ahcys_AHC1_ahcys["model"] = "missing ahcys_c measurement"
    AHC1_ahcys_df = pd.concat([met_6_AHC1_ahcys, met_6_missing_ahcys_AHC1_ahcys]).reset_index()
    AHC1_ahcys_df["ln_km"] = AHC1_ahcys_df["km"].apply(lambda x: np.log(x))


    for i, row in concs_validation.iterrows():
        concs_validation.loc[i, "measurement"] = mi_6.measurements.yconc.loc[row["experiments"], row["mics"]].measurement

    concs_validation[0.05] = concs_validation.apply(lambda x: np.log(x[0.05]) - np.log(x["measurement"]), axis=1)
    concs_validation[0.95] = concs_validation.apply(lambda x: np.log(x[0.95]) - np.log(x["measurement"]), axis=1)
    # calling plotting scripts
    cm = 1/2.54  # centimeters in inches
    sns.set_context("paper", rc={"font.size":14,"axes.titlesize":14,"axes.labelsize":14})

    fig = plt.figure(constrained_layout=True, figsize = (21.0*cm, 21.0*cm))
    gs = GridSpec(2, 2, figure=fig)

    samples_legend = mpatches.Patch(color=colour_scheme["complete_measurement_set"], label="complete measurement set")
    measurements_legend = mpatches.Patch(color=colour_scheme["missing_ahcys_measurement"], label="missing ahcys measurement")
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

    ax0 = plot_residual_state_variables(concs_validation, ax0)
    ax1 = plot_parameter_hist(AHC1_ahcys_df, ax1)
    ax2 = plot_box(df_ll, ax2)
    fig.savefig("plots/figure_3_5.png")

if __name__ == "__main__":
    plot_figure()
