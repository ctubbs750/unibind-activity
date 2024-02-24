""""""

import re
from sklearn import metrics
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

    # Try summary
    try:
        # Get count summary info
        sites[[0, 1, 2, 3]] = (
            sites["summary"]
            .str.replace(",", ":")
            .str.split(":", expand=True)
            .fillna(0)
            .astype(int)
        )
        # Whenever value in col 0 is 1, then switch cols 0,1 with 2,3
        sites.loc[sites[0] == 1, [0, 1, 2, 3]] = sites.loc[
            sites[0] == 1, [2, 3, 0, 1]
        ].values

        # Drop and rename
        sites = sites.drop(columns=[0, 2]).rename(
            columns={1: "n_unbound", 3: "n_bound"}
        )
    # This is when there is no bound - not exactly true below but quick fix for now #NOTE
    except ValueError:
        sites["n_unbound"] = 0
        sites["n_bound"] = 0

    # Total motifs
    sites["n_motifs"] = sites["n_bound"] + sites["n_unbound"]

    # Proportion bound - Activity
    sites["activity"] = sites["n_bound"] / sites["n_motifs"]

    # Add standard error on proportion
    sites["ci95_lbound"], sites["ci95_rbound"] = proportion_confint(
        sites["n_bound"], sites["n_motifs"]
    )

    # Rolling motif count
    sites["nmotif_roll"] = sites["n_motifs"].rolling(window=WINDOW).mean()

    # Try, exept index error means motif count is good across the board, so can skip the process
    try:
        last_index = sites[sites["nmotif_roll"] < THRESH].index[0] - WINDOW

        # handle negatives
        if last_index < 0:
            last_index = 0

        # Update to values at the cap index
        columns_to_update = [
            "n_unbound",
            "n_bound",
            "n_motifs",
            "activity",
            "ci95_lbound",
            "ci95_rbound",
            "nmotif_roll",
        ]

        # Truncate
        for column in columns_to_update:
            sites.loc[last_index:, column] = sites.loc[last_index, column]

        # Flat that you adjusted
        sites["adjusted"] = [1 if x > last_index else 0 for x in sites.index]

        # Bound to 0 and 1
        x = sites["perc"] / 10
        y = sites["activity"]

        # Get AUC
        auc = metrics.auc(x, y)

        # Write to file
        sites.to_csv(OUTPUT, sep="\t", index=False)

        # Save auc
        with open(OUTPUT.replace(".tsv", ".auc"), "w") as f:
            f.write(str(auc))

    except IndexError:
        # Flag not adjusted
        sites["adjusted"] = 0

        # Bound to 0 and 1
        x = sites["perc"] / 10
        y = sites["activity"]

        # Get AUC
        auc = metrics.auc(x, y)

        # Write to file
        sites.to_csv(OUTPUT, sep="\t", index=False)

        # Save auc
        with open(OUTPUT.replace(".tsv", ".auc"), "w") as f:
            f.write(str(auc))


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
