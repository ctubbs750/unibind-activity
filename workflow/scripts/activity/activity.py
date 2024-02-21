"""d"""

from sklearn import metrics
from pandas import DataFrame, read_csv, merge
from statsmodels.stats.proportion import proportion_confint

# Snakemake
SITES = snakemake.input[0]  # type: ignore
WINDOW = snakemake.params.window  # type: ignore
THRESH = snakemake.params.threshold  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_ix_sites(filepath: str) -> DataFrame:
    """Reads ix sites fast"""
    return read_csv(
        filepath,
        sep="\t",
        dtype={"score": int, "bound": int, "pval": float, "perc": float},
        compression="gzip",
        engine="c",
    )


def main() -> None:
    """D"""
    # Read in
    ix_sites = read_ix_sites(SITES)

    # Convert bound to binary
    ix_sites["bound"] = ix_sites["bound"].apply(lambda x: 1 if x > 0 else 0)

    # Summary - Nans are no bound or all bound
    activity = (
        ix_sites.groupby("score")["bound"]
        .value_counts(normalize=False)
        .unstack()
        .fillna(0)
        .rename(columns={0: "n_unbound", 1: "n_bound"})
    ).reset_index()
    
    # Total motifs
    activity["n_motifs"] = activity["n_bound"] + activity["n_unbound"]

    # Proportion bound - Activity
    activity["activity"] = activity["n_bound"] / activity["n_motifs"]

    # Add standard error on proportion
    activity["ci95_lbound"], activity["ci95_rbound"] = proportion_confint(
        activity["n_bound"], activity["n_motifs"]
    )

    # Rolling average of CI width
    activity["ci_width_roll"] = (
        (activity["ci95_rbound"] - activity["ci95_lbound"])
        .rolling(window=WINDOW)
        .mean()
    )
    
    # Reset index for adjustments
    activity.reset_index(inplace=True)

    # Find first row where CI width is less than 0.1
    last_index = activity[activity["ci_width_roll"] < THRESH].index[0] - WINDOW

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
        "ci_width_roll",
    ]

    # Truncate
    for column in columns_to_update:
        activity.loc[last_index:, column] = activity.loc[last_index, column]

    # Flat that you adjusted
    activity["adjusted"] = [1 if x > last_index else 0 for x in activity.index]

    # Dict where score are keys and values are perc
    score_perc = dict(zip(ix_sites["score"], ix_sites["perc"]))

    # Update with perc
    activity['perc'] = activity['score'].map(score_perc)

    # Bound to 0 and 1
    x = activity["perc"] / 100
    y = activity["activity"]

    # Get AUC
    auc = metrics.auc(x, y)

    # Write to file
    activity.to_csv(OUTPUT, sep="\t", index=False)

    # Save auc
    with open(OUTPUT.replace(".tsv", ".auc"), "w") as f:
        f.write(str(auc))


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
