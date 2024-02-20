from pathlib import Path
from Bio import motifs
from pandas import DataFrame
from Bio.motifs.jaspar import calculate_pseudocounts
from multiprocessing import Pool
import logomaker
from matplotlib import pyplot as plt

# Snakemake parameters
IP_DIR = snakemake.input[0]  # type: ignore
OP_DIR = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def make_output_dir(op_dir: str) -> None:
    """Make output directory if it doesn't exist."""
    Path(op_dir).mkdir(parents=True, exist_ok=True)


def get_pfms(ip_dir: str) -> list:
    """Return list of fasta files from input directory."""
    return [i for i in Path(ip_dir).glob("*.pfm")]


def read_pfm(filepath: str):
    """Read PFM from from file and return transposed matrix for Biopython"""
    with open(filepath) as handle:
        m = motifs.read(handle, "pfm")
    return m


def save_logo(pssm_data: DataFrame, outpath: str) -> None:
    """Save logo plot"""

    pssm_data.columns = ["A", "C", "G", "T"]
    # Mask all values below 0
    pssm_data = pssm_data.mask(pssm_data < 0, 0)

    # # Get index in 1-base
    pssm_data.index = [i + 1 for i in pssm_data.index]

    # Create Logo object
    pwm_logo = logomaker.Logo(pssm_data, color_scheme="classic", ax=None)

    # Style using Logo methods
    pwm_logo.style_spines(visible=False)
    pwm_logo.style_spines(spines=["left", "bottom"], visible=True)
    pwm_logo.style_xticks(rotation=0, fmt="%d", anchor=0)

    # Style using Axes methods
    pwm_logo.ax.set_ylabel("Bits", labelpad=-1)
    pwm_logo.ax.xaxis.set_ticks_position("none")
    pwm_logo.ax.xaxis.set_tick_params(pad=-1)
    pwm_logo.ax.set_title("d", color="r")

    # # Save figure
    plt.savefig(outpath, dpi=100)


def save_pssm(pfm: str, outdir: str) -> None:
    """Main program"""
    # Read input PFM
    pfm = read_pfm(pfm)

    # Calculate pseudocounts
    pfm.pseudocounts = calculate_pseudocounts(pfm)

    # Caclulate PSSM
    pssm = list(map(list, zip(*[pfm.pssm[nt] for nt in "ACGT"])))  # type: ignore

    # Grab IC and GC content
    ic = pfm.pssm.mean()
    gc = pfm.counts.gc_content

    # Make intlogodds PSSM
    intlogodds = (DataFrame.from_dict(pssm) * 100).round().astype(int)

    # Format ic and gc as stats file
    stats_df = DataFrame.from_dict({"IC": [ic], "GC": [gc]})

    # Make dir if doesn't exist
    make_output_dir(outdir)

    # Write outputs
    save_logo((DataFrame.from_dict(pssm)), f"{outdir}/logo.png")
    (DataFrame.from_dict(pssm)).to_csv(
        f"{outdir}/pssm", sep="\t", header=False, index=False
    )
    intlogodds.to_csv(f"{outdir}/pssm.intLogOdds", sep="\t", header=False, index=False)
    stats_df.to_csv(f"{outdir}/pssm.stats", sep="\t", header=True, index=False)


# Function to save PFM for each FASTA
def save_pssm_for_pfm(pfm):
    # Set output file, calculate and save
    output_pssm_dir = f"{OP_DIR}/{'.'.join(Path(pfm).name.split('.')[:-1])}"
    save_pssm(pfm, output_pssm_dir)


def main() -> None:
    """Main program"""
    # Make initial output directory
    make_output_dir(OP_DIR)

    # PFM files
    pfm_files = get_pfms(IP_DIR)

    # Number of processes
    num_processes = snakemake.threads  # type: ignore

    # Parallelize
    with Pool(processes=num_processes) as p:
        p.map(save_pssm_for_pfm, pfm_files)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
