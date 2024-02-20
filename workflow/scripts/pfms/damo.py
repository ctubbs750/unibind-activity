from pathlib import Path
from pandas import DataFrame
from Bio import SeqIO, motifs
from multiprocessing import Pool

# Snakemake parameters
IP_DIR = snakemake.input[0]  # type: ignore
OP_DIR = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def make_output_dir(op_dir: str) -> None:
    """Make output directory if it doesn't exist."""
    Path(op_dir).mkdir(parents=True, exist_ok=True)


def get_fasta_files(ip_dir: str) -> list:
    """Return list of fasta files from input directory."""
    return [i for i in Path(ip_dir).glob("*.fa")]


def save_pfm(fasta: str, output: str) -> None:
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


# Function to save PFM for each FASTA
def save_pfm_for_fasta(fasta):
    # Set output file, calculate and save
    output_pfm = f"{OP_DIR}/{'.'.join(fasta.name.split('.')[:-1])}.pfm"
    save_pfm(fasta, output_pfm)


def main() -> None:
    """Main program"""
    # Make output directory
    make_output_dir(OP_DIR)

    # Fasta files
    fasta_files = get_fasta_files(IP_DIR)

    # Number of processes
    num_processes = snakemake.threads  # type: ignore

    # Parallelize
    with Pool(processes=num_processes) as p:
        p.map(save_pfm_for_fasta, fasta_files)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
