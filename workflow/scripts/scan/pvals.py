""""""

from pandas import read_csv, DataFrame, merge


# Snakemake parameters
IP_SITES = snakemake.input[0]  # type: ignore
IP_PVALS = snakemake.input[1]  # type: ignore
OP_SITES = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def read_sites(filepath: str) -> DataFrame:
    """Reads the sites file"""
    return read_csv(
        filepath, sep="\s+", header=None, names=["score", "bound"], compression="gzip"
    )


def read_pvals(filepath: str) -> DataFrame:
    """Reads the pvals file"""
    return read_csv(filepath, sep="\s+", header=None, names=["score", "pval", "perc"])


def main() -> None:
    """Main function"""
    # Read in the sites file
    sites = read_sites(IP_SITES)
    pvals = read_pvals(IP_PVALS)

    # Merge the dataframes
    merged = merge(sites, pvals, on="score", how="left")

    # Write the merged dataframe
    merged.to_csv(OP_SITES, sep="\t", index=False, compression="gzip")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
