import requests
import tempfile
from pathlib import Path
from pandas import DataFrame
from Bio import motifs

# Snakemake parameters
IP_DIR = snakemake.input[0]  # type: ignore
OP_DIR = snakemake.output[0]  # type: ignore
PFM_URL = snakemake.params[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def get_unique_profiles(ip_dir: str) -> list:
    """Get unique profiles from input directory."""
    # Get unique profiles
    profiles = set()
    for file in Path(ip_dir).glob("*.pfm"):
        profile = ".".join(file.name.split(".")[-4:-2])
        profiles.add(profile)
    return list(profiles)


def make_output_dir(op_dir: str) -> None:
    """Make output directory if it doesn't exist."""
    Path(op_dir).mkdir(parents=True, exist_ok=True)


def main() -> None:
    """Main program"""
    # Make output directory
    make_output_dir(OP_DIR)

    # Get unique profiles
    profiles = get_unique_profiles(IP_DIR)

    # Save URL to temp file
    with tempfile.NamedTemporaryFile(
        suffix=".jaspar", delete=False, mode="w+"
    ) as temp_file:
        temp_file.write(requests.get(PFM_URL).text)

    # Connect to file
    fh = open(temp_file.name)

    # Loop over motifs
    for m in motifs.parse(fh, "jaspar"):
        # Check if motif is in unibind profiles
        if m.matrix_id not in profiles:
            pass

        # Get output file name
        output_file = f"{OP_DIR}/{m.matrix_id}.pfm"

        # Get counts
        counts = m.counts

        # Format and save counts as PFM
        pfm = DataFrame.from_dict(counts).T

        # Save PFM
        pfm.to_csv(output_file, sep="\t", header=False, index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
