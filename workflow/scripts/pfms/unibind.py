from pandas import DataFrame
from Bio import SeqIO, motifs

# Snakemake parameters
IP_FA = snakemake.input[0]  # type: ignore
OP_PFM = snakemake.output[0]  # type: ignore
PROFILE = snakemake.params.profile  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def save_pfm(fasta: str = IP_FA, output: str = OP_PFM) -> None:
    """Main program"""
    # Report binding sequences
    instances = []
    for record in SeqIO.parse(fasta, "fasta"):
        # Check if seq has masked nucleotides, if so skip
        if "N" in record.seq or "n" in record.seq:
            continue
        instances.append(record.seq.upper())

    # Create motif from instances
    m = motifs.create(instances, alphabet="ACGT")

    # Get counts
    counts = m.counts

    # Format and save counts as PFM
    pfm = DataFrame.from_dict(counts).T
    pfm.to_csv(output, sep="\t", header=False, index=False)


def main() -> None:
    """Main program"""
    # Save pfm
    save_pfm(IP_FA, OP_PFM)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
