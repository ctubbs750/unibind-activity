from pandas import DataFrame
from Bio import motifs

# Snakemake parameters
IP_PFM = snakemake.input[0]  # type: ignore
OP_PFM = snakemake.output[0]  # type: ignore
PROFILE = snakemake.params.profile  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def save_pfm(ip_jaspar: str, output: str) -> None:
    """Main program"""
    # Some motifs from UniBind aren't in JASPAR download, ok to toss
    # Connect to file
    fh = open(ip_jaspar)

    # Loop over motifs
    for m in motifs.parse(fh, "jaspar"):
        # Get counts
        counts = m.counts

        # Format and save counts as PFM
        pfm = DataFrame.from_dict(counts).T

        # Save PFM
        pfm.to_csv(output, sep="\t", header=False, index=False)


def main() -> None:
    """Main program"""
    # Save pfm
    save_pfm(IP_PFM, OP_PFM)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
