""""""

import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import LinearRegression
from pandas import DataFrame, read_csv, merge
from statsmodels.stats.proportion import proportion_confint


# Snakemake parameters
SITES = snakemake.input[0]  # type: ignore
PVALS = snakemake.input[1]  # type: ignore
WINDOW = snakemake.params.window  # type: ignore
THRESH = snakemake.params.threshold  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def read_sites(filepath: str) -> DataFrame:
    """Reads ix sites fast"""
    return read_csv(
        filepath,
        sep="\t",
        dtype={"score": int, "summary": str},
        engine="c",
        names=["score", "summary"],
    )


def read_pvals(filepath: str) -> DataFrame:
    """Reads the pvals file"""
    return read_csv(filepath, sep="\s+", header=None, names=["score", "pval", "perc"])


def main() -> None:
    """Main function"""
    # Read in the sites file
    sites = read_sites(SITES)
    pvals = read_pvals(PVALS)

    # Merge with pvals
    sites = merge(sites, pvals, on="score", how="left")

    # Read score and summary and make dict of (un)bound:count
    summary_dict = {
        score: dict([i.split(":") for i in summary.split(",")])
        for score, summary in zip(sites["score"], sites["summary"])
    }

    # Add '0' or '1' to each key if not present
    for key in summary_dict.keys():
        if "0" not in summary_dict[key].keys():
            summary_dict[key]["0"] = 0
        if "1" not in summary_dict[key].keys():
            summary_dict[key]["1"] = 0

    # Add n_unbound and n_bound to sites
    sites["n_unbound"] = [int(summary_dict[s]["0"]) for s in sites["score"]]
    sites["n_bound"] = [int(summary_dict[s]["1"]) for s in sites["score"]]

    # Total motifs
    sites["n_motifs"] = sites["n_bound"] + sites["n_unbound"]

    # Proportion bound - Activity
    sites["activity"] = sites["n_bound"] / sites["n_motifs"]

    # Add standard error on proportion - wilson interval to help with small proportions
    sites["ci95_lbound"], sites["ci95_rbound"] = proportion_confint(
        sites["n_bound"], sites["n_motifs"], method="wilson"
    )

    # Width of CI
    sites["ci_width"] = sites["ci95_rbound"] - sites["ci95_lbound"]

    # Rolling motif count
    sites["ci_width_roll"] = sites["ci_width"].rolling(window=WINDOW).mean()

    # This may be driven by small proportions of 1 - weightign by CI width, Fill with min - these are missed in the rolling
    weights = sites["ci_width_roll"].fillna(sites["ci_width_roll"].min())

    # Train the linear regression model with weights
    X = np.array(sites["score"]).reshape(-1, 1)
    y = sites["activity"]
    model = LinearRegression().fit(X, y, sample_weight=weights)

    # Get the predicted values
    y_pred = model.predict(X)

    # Calculate the Mean Squared Error
    mse = mean_squared_error(y, y_pred)

    # Write to file
    sites.to_csv(OUTPUT, sep="\t", index=False)

    # Save auc
    with open(OUTPUT.replace(".tsv", ".mse"), "w") as f:
        f.write(str(mse))


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
